"""Check the contents of a clean, public release checkout without printing matches."""

import argparse
from collections import Counter
import csv
import io
import json
from pathlib import Path
import re
import stat


ROOT_FILES = {
    ".Renviron.example",
    ".gitattributes",
    ".gitignore",
    "CITATION.cff",
    "LICENSE",
    "README.md",
}
SOURCE_TYPES = {".py", ".r", ".sh"}
SUMMARY_TYPES = {".csv", ".json"}
SIZE_LIMIT = 50_000_000
JWT = re.compile(r"\b" + "ey" + r"J[A-Za-z0-9_-]{6,}")
HOST_PATH = re.compile(
    r"(?i)(?<![\w])(?:[a-z]:[/\\]|/(?:root|home|mnt/[a-z]|Users|tmp)/|"
    r"\\\\[^\\\s]+\\[^\\\s]+)"
)
OTHER_CREDENTIALS = re.compile(
    r"""(?ix)
    \b(?:gh[pousr]_[a-z0-9]{20,}|github_pat_[a-z0-9_]{30,}|AKIA[A-Z0-9]{16})\b
    | -----BEGIN\ (?:[A-Z]+\ )?PRIVATE\ KEY-----
    | \bBearer[ \t]+[A-Za-z0-9._/+~=-]{8,}
    | ://[^/\s"'<>]+:[^/\s"'<>]+@
    | \b(?:https?|socks5h?)://(?:localhost|127\.0\.0\.1|0\.0\.0\.0|10\.\d+\.\d+\.\d+|192\.168\.\d+\.\d+):\d+
    """
)
SECRET_ASSIGNMENT = re.compile(
    r"""(?imx)\b(?:[A-Za-z0-9_]*_(?:jwt|token|password|secret|api_key)|"""
    r"""jwt|token|password|secret|api_key)\s*(?:<-|=|:)\s*["']([^"'\n]{8,})["']"""
)
EMAIL = re.compile(r"(?i)\b[A-Za-z0-9_.%+-]+@[A-Za-z0-9.-]+\.[A-Z]{2,}\b")
CONTACT = "luzhongxin" + "@zxhospital.com"
SAMPLE_ID = re.compile(
    r"(?i)\b(?:TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}(?:-[A-Z0-9]+)*|"
    r"(?:GSM|SRR|SRS|SAMN)\d{5,}|[ACGT]{12,}-\d+)\b"
)
INDIVIDUAL_COLUMNS = {
    "barcode", "cellid", "dateofbirth", "dob", "firstname", "lastname",
    "origident", "patient", "patientid", "participantid", "runid",
    "sampleid", "subjectid",
}
PRIVATE_SUMMARIES = {
    "gmsig4_ici_subgroup_exactscores.csv",
    "scrna_tnbc_patient_pseudobulk_and_props.csv",
    "scrna_tnbc_patientlevel_cellcomm.csv",
    "scrna_tnbc_umap_metadata_keygenes.csv",
}
CHECKS = {
    "renviron": "Real .Renviron files",
    "jwt": "JWT-like strings",
    "host_path": "Hard-coded machine paths",
    "oversize": "Files larger than 50 MB",
    "manuscript": "Manuscript files or unapproved Markdown",
    "credential": "Other credentials or private proxy endpoints",
    "email": "Unapproved email addresses",
    "individual": "Result files with participant or cell identifiers",
    "allowlist": "Non-allowlisted files, directories, or links",
    "format": "Invalid text, metadata, or result files",
}


def scan(root):
    counts = Counter()
    findings = []
    inventory = []

    def record(check, path, matches=1):
        if matches:
            counts[check] += matches
            findings.append((CHECKS[check], path, matches))

    def visit(directory):
        for entry in sorted(directory.iterdir(), key=lambda path: path.name.casefold()):
            rel = entry.relative_to(root).as_posix()
            if entry.is_symlink() or entry.is_junction():
                record("allowlist", rel)
                continue
            attributes = entry.lstat()
            if getattr(attributes, "st_file_attributes", 0) & stat.FILE_ATTRIBUTE_REPARSE_POINT:
                record("allowlist", rel)
                continue
            if entry.is_dir():
                if rel == ".git":
                    continue
                if rel not in {"00_scripts", "02_results"}:
                    record("allowlist", rel)
                visit(entry)
                continue
            if not entry.is_file():
                record("allowlist", rel)
                continue

            size = attributes.st_size
            inventory.append((rel, size))
            suffix = entry.suffix.lower()
            name = entry.name.lower()
            record("renviron", rel, int(name == ".renviron"))
            record("manuscript", rel, int(suffix in {".docx", ".docm"} or (
                suffix == ".md" and rel != "README.md"
            )))
            allowed = rel in ROOT_FILES or (
                entry.parent == root / "00_scripts" and suffix in SOURCE_TYPES
            ) or (
                entry.parent == root / "02_results" and suffix in SUMMARY_TYPES
            )
            record("allowlist", rel, int(not allowed))
            if size > SIZE_LIMIT:
                record("oversize", rel)
                continue
            try:
                content = entry.read_bytes().decode("utf-8-sig")
                if "\0" in content:
                    raise ValueError("NUL in text")
            except (UnicodeError, ValueError):
                record("format", rel)
                continue

            record("jwt", rel, len(JWT.findall(content)))
            record("host_path", rel, len(HOST_PATH.findall(content)))
            record("credential", rel, len(OTHER_CREDENTIALS.findall(content)))
            record("credential", rel, len([
                value for value in SECRET_ASSIGNMENT.findall(content)
                if value != "your_opengwas_jwt_token_here"
            ]))
            record("email", rel, len([
                email for email in EMAIL.findall(content)
                if email != CONTACT or rel not in {"CITATION.cff", "README.md"}
            ]))

            if entry.parent == root / "02_results":
                inspect_result(content, name, rel, record)
            if rel == ".Renviron.example":
                expected = ["OPENGWAS_JWT=your_opengwas_jwt_token_here", "PROJECT_ROOT=."]
                record("format", rel, int(content.splitlines() != expected))

    visit(root)
    present = {rel for rel, _ in inventory}
    for missing in sorted(ROOT_FILES - present):
        record("format", missing)
    for directory in ("00_scripts", "02_results"):
        if not any(rel.startswith(directory + "/") for rel in present):
            record("format", directory)
    return counts, inventory, findings


def inspect_result(content, name, relative_path, record):
    try:
        if name.endswith(".json"):
            json.loads(content)
        elif name.endswith(".csv"):
            reader = csv.reader(io.StringIO(content, newline=""), strict=True)
            fields = next(reader)
            normalized = {re.sub(r"[^a-z0-9]", "", field.lower()) for field in fields}
            if normalized & INDIVIDUAL_COLUMNS or name in PRIVATE_SUMMARIES or name.startswith("gmsig4_patientscores_"):
                record("individual", relative_path)
            if any(len(row) != len(fields) for row in reader):
                record("format", relative_path)
    except (csv.Error, json.JSONDecodeError, StopIteration):
        record("format", relative_path)
    if SAMPLE_ID.search(content):
        record("individual", relative_path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "root", nargs="?", type=Path, default=Path(__file__).resolve().parents[1],
        help="Clean release checkout. The current checkout is used by default.",
    )
    args = parser.parse_args()
    counts, inventory, findings = scan(args.root.resolve(strict=True))
    total = sum(size for _, size in inventory)
    scripts = sum(rel.startswith("00_scripts/") for rel, _ in inventory)
    results = sum(rel.startswith("02_results/") for rel, _ in inventory)
    print(f"Release files: {len(inventory)} ({scripts} scripts, {results} result tables, "
          f"{len(inventory) - scripts - results} root files)")
    print(f"Release size: {total:,} bytes ({total / 1_000_000:.2f} MB)")
    print(f"Largest file: {max((size for _, size in inventory), default=0):,} bytes")
    for check, label in CHECKS.items():
        print(f"{label}: {counts[check]}")
    for label, relative_path, count in findings:
        print(f"{label}: {relative_path} ({count})")
    print("Security scan: " + ("FAIL" if findings else "PASS"))
    return 1 if findings else 0


if __name__ == "__main__":
    raise SystemExit(main())
