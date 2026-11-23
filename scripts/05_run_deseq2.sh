#!/usr/bin/env bash
# Step 5: differential expression, pathway analysis, run QC, and dashboard data (runs deseq2.R).
set -uo pipefail
# shellcheck source=../config.sh
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/config.sh"
need_cmds Rscript

[ -s "$COUNTS_FILE" ] || { echo "ERROR: counts table not found: $COUNTS_FILE (run step 04)." >&2; exit 1; }
cd "$PROJECT_DIR" || exit 1
export COUNTS_FILE BAM_DIR
Rscript scripts/deseq2.R
