#!/usr/bin/env bash
# Run every step in order. Each step skips work that is already done, so it is safe to re-run.
# Usage:  bash scripts/run_all.sh            # all steps
#         bash scripts/run_all.sh 4          # start from step 4 (counts) onward
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
START="${1:-1}"
steps=("01_download_sra.sh" "02_trim_fastp.sh" "03_align_hisat2.sh" "04_count_featurecounts.sh" "05_run_deseq2.sh")
for i in "${!steps[@]}"; do
    n=$((i + 1))
    [ "$n" -lt "$START" ] && continue
    echo; echo "######## Step $n: ${steps[$i]} ########"
    bash "$HERE/${steps[$i]}" || { echo "Step $n failed; stopping." >&2; exit 1; }
done
echo; echo "Done. Launch the dashboard with:  Rscript -e 'shiny::runApp(\"app\")'"
