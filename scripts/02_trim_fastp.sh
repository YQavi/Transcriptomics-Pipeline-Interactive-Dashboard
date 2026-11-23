#!/usr/bin/env bash
# Step 2: adapter detection, quality filtering and per-sample QC reports with fastp.
# Skips samples whose cleaned FASTQ files already exist.
set -uo pipefail
# shellcheck source=../config.sh
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/config.sh"
need_cmds fastp

mkdir -p "$CLEAN_DIR" "$REPORT_DIR"
FAILED=()

for sample in "${SAMPLES[@]}"; do
    echo "=== fastp $sample at $(date) ==="
    in1="$RAW_DIR/${sample}_1.fastq.gz";  out1="$CLEAN_DIR/${sample}_1.clean.fastq.gz"
    in2="$RAW_DIR/${sample}_2.fastq.gz";  out2="$CLEAN_DIR/${sample}_2.clean.fastq.gz"

    if [ "$PAIRED" = true ]; then
        if [ -s "$out1" ] && [ -s "$out2" ]; then echo "Already trimmed, skipping."; continue; fi
        if [ ! -f "$in1" ] || [ ! -f "$in2" ]; then
            echo "ERROR: missing input for $sample ($in1 / $in2)" >&2; FAILED+=("$sample (missing fastq)"); continue
        fi
        # FASTP_EXTRA is intentionally unquoted so its options split into separate words
        # shellcheck disable=SC2086
        if ! fastp -i "$in1" -I "$in2" -o "$out1" -O "$out2" \
                -h "$REPORT_DIR/${sample}_fastp.html" -j "$REPORT_DIR/${sample}_fastp.json" \
                $FASTP_EXTRA --detect_adapter_for_pe -w "$THREADS_FASTP"; then
            rm -f "$out1" "$out2"; FAILED+=("$sample (fastp)")
        fi
    else
        if [ -s "$out1" ]; then echo "Already trimmed, skipping."; continue; fi
        if [ ! -f "$in1" ]; then
            echo "ERROR: missing input for $sample ($in1)" >&2; FAILED+=("$sample (missing fastq)"); continue
        fi
        # shellcheck disable=SC2086
        if ! fastp -i "$in1" -o "$out1" \
                -h "$REPORT_DIR/${sample}_fastp.html" -j "$REPORT_DIR/${sample}_fastp.json" \
                $FASTP_EXTRA -w "$THREADS_FASTP"; then
            rm -f "$out1"; FAILED+=("$sample (fastp)")
        fi
    fi
done

report_failures "${FAILED[@]+"${FAILED[@]}"}"
