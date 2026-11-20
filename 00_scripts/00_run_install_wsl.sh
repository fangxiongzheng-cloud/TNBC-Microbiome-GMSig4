#!/usr/bin/env bash
set -e
project_root="${PROJECT_ROOT:-.}"
exec Rscript "$project_root/00_scripts/00_install_r_packages.R"
