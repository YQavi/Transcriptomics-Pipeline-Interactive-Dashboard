#!/usr/bin/env bash
# Step 1: download SRA runs and convert to gzipped FASTQ.
# Skips samples whose FASTQ files already exist.
set -uo pipefail
# shellcheck source=../config.sh
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/config.sh"
need_cmds prefetch fasterq-dump gzip

mkdir -p "$RAW_DIR" "$SCRATCH_DIR"
if [ "$PAIRED" = true ]; then SUFFIXES=(1 2); else SUFFIXES=(1); fi
FAILED=()

for i in "${!SAMPLES[@]}"; do
    sample="${SAMPLES[$i]}"; acc="${ACCESSIONS[$i]}"
    echo "=== $sample ($acc) at $(date) ==="

    have_all=true
    for s in "${SUFFIXES[@]}"; do [ -s "$RAW_DIR/${sample}_${s}.fastq.gz" ] || have_all=false; done
    if [ "$have_all" = true ]; then echo "FASTQ already present, skipping."; continue; fi

    work="$SCRATCH_DIR/sra_${sample}"
    rm -rf "$work"; mkdir -p "$work"

    if ! prefetch "$acc" -O "$work"; then
        FAILED+=("$sample (prefetch)"); continue
    fi
    # --temp on internal disk avoids the slow external drive for small random writes
    if ! fasterq-dump --split-files --threads "$THREADS_SRA" --temp "$work" -O "$work" "$work/$acc"; then
        FAILED+=("$sample (fasterq-dump)"); continue
    fi

    # single end runs are written as <acc>.fastq (no _1 suffix)
    if [ -f "$work/${acc}.fastq" ]; then mv "$work/${acc}.fastq" "$work/${acc}_1.fastq"; fi

    ok=true
    for s in "${SUFFIXES[@]}"; do
        if [ ! -f "$work/${acc}_${s}.fastq" ]; then echo "ERROR: missing ${acc}_${s}.fastq" >&2; ok=false; break; fi
        # write to a temp name, then rename, so a half-written file is never mistaken for done
        gzip -c "$work/${acc}_${s}.fastq" > "$RAW_DIR/${sample}_${s}.fastq.gz.part" \
            && mv "$RAW_DIR/${sample}_${s}.fastq.gz.part" "$RAW_DIR/${sample}_${s}.fastq.gz" || ok=false
    done
    [ "$ok" = true ] || FAILED+=("$sample (convert)")
    rm -rf "$work"
done

report_failures "${FAILED[@]+"${FAILED[@]}"}"
