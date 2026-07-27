#!/usr/bin/env bash
# Run the whole project from a clean clone, in order.
#
#   ./run_all.sh          use cached simulation grids where available
#   REFRESH=TRUE ./run_all.sh   re-simulate everything (~90 min on 11 cores)
#
# set -e matters: the phases depend on each other, and an earlier failure that
# lets later steps run produces a report full of stale numbers. This project
# already got bitten once by a `;`-chained pipeline that carried on past a crash.
set -euo pipefail
cd "$(dirname "$0")"

step () { printf '\n\033[1m== %s\033[0m\n' "$1"; }

step "Phase 0: simulation (no real data touched)"
Rscript sim/01_ideal.R
Rscript sim/02_null_calibration.R
Rscript sim/03_degrade.R
Rscript sim/04_power.R
Rscript sim/05_measurement_error.R
Rscript sim/07_likert_artifact.R
Rscript sim/06_summary.R | tee output/tables/06_summary.txt

step "Phase 1: download + reproduce"
Rscript analysis/00_download.R
Rscript analysis/01_clean.R
Rscript analysis/02_reproduce.R
Rscript analysis/03_artifact_test.R

step "Phase 2: specification curve (preregistered)"
Rscript analysis/04_speccurve.R

step "Phase 3: design recommendations (preregistered)"
Rscript analysis/05_design.R

step "Report"
quarto render report/report.qmd --to html

printf '\n\033[1mDone.\033[0m report/report.html\n'
