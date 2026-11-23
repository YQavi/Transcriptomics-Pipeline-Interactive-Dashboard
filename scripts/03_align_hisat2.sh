#!/usr/bin/env bash
# Step 3: splice-aware alignment with HISAT2, piped straight into samtools sort (no SAM on disk).
# Skips samples whose sorted BAM and index already exist.
# Not using `set -e`: we log a failure per sample and keep going.
set -uo pipefail
# shellcheck source=../config.sh
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/config.sh"
need_cmds hisat2 samtools

# Fail early and clearly if the index is missing (e.g. external drive not mounted)
if ! compgen -G "${INDEX_PREFIX}.*.ht2*" > /dev/null; then
    echo "ERROR: no HISAT2 index files found for prefix '$INDEX_PREFIX' (expected ${INDEX_PREFIX}.1.ht2 ...)." >&2
    echo "       Check INDEX_PREFIX in config.sh and that the drive is mounted." >&2
    exit 1
fi

mkdir -p "$BAM_DIR" "$SCRATCH_DIR"
AVAIL_GB=$(free_gb "$BAM_DIR")
if [ "$AVAIL_GB" -lt 20 ]; then
    echo "WARNING: only ${AVAIL_GB}GB free where BAMs are written. Sorting large BAMs may fail mid-run." >&2
fi

# HISAT2 strandedness flag matching featureCounts -s (0 none, 1 stranded, 2 reversely stranded)
STRAND_OPT=""
if [ "$PAIRED" = true ]; then fwd="FR"; rev="RF"; else fwd="F"; rev="R"; fi
[ "$STRANDEDNESS" = "1" ] && STRAND_OPT="--rna-strandness $fwd"
[ "$STRANDEDNESS" = "2" ] && STRAND_OPT="--rna-strandness $rev"

FAILED=()
for sample in "${SAMPLES[@]}"; do
    echo "=========================================="
    echo "HISAT2 alignment for $sample at $(date)"
    echo "=========================================="

    bam="$BAM_DIR/${sample}.sorted.bam"
    if [ -s "$bam" ] && [ -s "${bam}.bai" ]; then echo "BAM and index exist, skipping."; continue; fi

    R1="$CLEAN_DIR/${sample}_1.clean.fastq.gz"; R2="$CLEAN_DIR/${sample}_2.clean.fastq.gz"
    if [ "$PAIRED" = true ]; then need_files=("$R1" "$R2"); else need_files=("$R1"); fi
    missing=false
    for f in "${need_files[@]}"; do [ -f "$f" ] || { echo "ERROR: missing $f" >&2; missing=true; }; done
    if [ "$missing" = true ]; then FAILED+=("$sample (missing fastq)"); continue; fi

    tmp="$SCRATCH_DIR/sort_${sample}"; rm -rf "$tmp"; mkdir -p "$tmp"

    # STRAND_OPT is unquoted on purpose (it holds two words); all file paths are quoted.
    if [ "$PAIRED" = true ]; then
        # shellcheck disable=SC2086
        keep_awake hisat2 -p "$THREADS_HISAT2" -x "$INDEX_PREFIX" -1 "$R1" -2 "$R2" $STRAND_OPT \
            --summary-file "$BAM_DIR/${sample}_summary.txt" \
          | samtools sort -@ "$THREADS_SAMTOOLS" -m "$SAMTOOLS_MEM_PER_THREAD" \
              -T "$tmp/tmp" -o "${bam}.part" -
    else
        # shellcheck disable=SC2086
        keep_awake hisat2 -p "$THREADS_HISAT2" -x "$INDEX_PREFIX" -U "$R1" $STRAND_OPT \
            --summary-file "$BAM_DIR/${sample}_summary.txt" \
          | samtools sort -@ "$THREADS_SAMTOOLS" -m "$SAMTOOLS_MEM_PER_THREAD" \
              -T "$tmp/tmp" -o "${bam}.part" -
    fi

    # With pipefail, check BOTH sides of the pipe. Copy PIPESTATUS immediately and give each
    # element a default: bash 3.2 on macOS can throw "unbound variable" under `set -u`.
    PS_COPY=("${PIPESTATUS[@]}")
    hisat2_rc="${PS_COPY[0]:-1}"; sort_rc="${PS_COPY[1]:-1}"
    rm -rf "$tmp"
    if [ "$hisat2_rc" -ne 0 ] || [ "$sort_rc" -ne 0 ]; then
        echo "ERROR: hisat2/samtools pipeline failed for $sample." >&2
        rm -f "${bam}.part"; FAILED+=("$sample (align/sort)"); continue
    fi
    mv "${bam}.part" "$bam"

    echo "Indexing ${sample}.sorted.bam..."
    if ! samtools index -@ "$THREADS_SAMTOOLS" "$bam"; then
        FAILED+=("$sample (index)")
    fi
    echo "Finished $sample at $(date)"
done

rm -rf "$SCRATCH_DIR"
report_failures "${FAILED[@]+"${FAILED[@]}"}"
