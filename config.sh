#!/usr/bin/env bash
# shellcheck disable=SC2034,SC2207  # variables are used by the scripts that source this file; word splitting of IDs is intended
# Settings shared by every shell step. Edit this file only; the scripts in scripts/ never need changes.
# Sourced by the scripts, so it must not contain spaces around "=".

# Repo root (works from any directory, including paths that contain spaces)
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ---- Where data lives ------------------------------------------------------
# Layout under DATA_DIR:  <acc>_1.fastq.gz  clean_fastq/  aligned_bam/  reports/  hisat2_index/  gene_counts.txt
# This matches the folder layout you already have (the repo sits inside the project folder).
# For a fresh clone use instead:  DATA_DIR="$PROJECT_DIR/data"
DATA_DIR="$(cd "$PROJECT_DIR/.." && pwd)"

SAMPLE_SHEET="$PROJECT_DIR/config/samples.tsv"   # columns: sample, condition, accession
RAW_DIR="$DATA_DIR"                              # raw <sample>_1.fastq.gz / _2.fastq.gz
CLEAN_DIR="$DATA_DIR/clean_fastq"
BAM_DIR="$DATA_DIR/aligned_bam"
REPORT_DIR="$DATA_DIR/reports"
COUNTS_FILE="$DATA_DIR/gene_counts.txt"
SCRATCH_DIR="/tmp/rnaseq_scratch"                # fast internal disk for SRA dumps and sort files

# ---- Reference -------------------------------------------------------------
INDEX_PREFIX="$DATA_DIR/hisat2_index/hg38"       # needs hg38.1.ht2 ... hg38.8.ht2
GTF="$DATA_DIR/annotation/genes.gtf"             # only needed by step 04; must match the genome build

# ---- Library type (check the GEO sample pages / methods) --------------------
PAIRED=true            # true = paired end, false = single end
STRANDEDNESS=0         # featureCounts -s: 0 unstranded, 1 stranded, 2 reversely stranded
                       # Verify with the assignment rate printed at the end of step 04.

# ---- Resources --------------------------------------------------------------
THREADS_SRA=4
THREADS_FASTP=4        # fastp allows at most 16
THREADS_HISAT2=5
THREADS_SAMTOOLS=4
THREADS_FEATURECOUNTS=4
SAMTOOLS_MEM_PER_THREAD="1500M"   # peak sort memory = THREADS_SAMTOOLS x this
FASTP_EXTRA="-q 20"               # --detect_adapter_for_pe is added automatically for paired end

# ---- Samples come from the sample sheet (single source of truth) ------------
if [ ! -f "$SAMPLE_SHEET" ]; then
    echo "ERROR: sample sheet not found: $SAMPLE_SHEET" >&2
    exit 1
fi
SAMPLES=($(tail -n +2 "$SAMPLE_SHEET" | cut -f1 | tr -d '\r'))
# accession column (3rd); falls back to the sample name when empty
ACCESSIONS=($(tail -n +2 "$SAMPLE_SHEET" | tr -d '\r' | awk -F'\t' '{print ($3 == "" ? $1 : $3)}'))
if [ "${#SAMPLES[@]}" -eq 0 ]; then
    echo "ERROR: no samples found in $SAMPLE_SHEET" >&2
    exit 1
fi

# ---- Helpers -----------------------------------------------------------------
need_cmds() {
    for cmd in "$@"; do
        if ! command -v "$cmd" >/dev/null 2>&1; then
            echo "ERROR: '$cmd' not found on PATH. Activate your conda environment (e.g. 'conda activate txnomics')." >&2
            exit 1
        fi
    done
}
# Free space in GB on the filesystem holding $1 (portable: works on macOS and Linux)
free_gb() { df -Pk "$1" | awk 'NR==2 {printf "%d", $4/1024/1024}'; }
# Prefix a command with caffeinate on macOS so the laptop does not sleep mid-run
keep_awake() { if command -v caffeinate >/dev/null 2>&1; then caffeinate -i -s "$@"; else "$@"; fi; }
report_failures() {
    # usage: report_failures "${FAILED[@]+"${FAILED[@]}"}"   (safe for empty arrays in bash 3.2)
    echo "=========================================="
    if [ "$#" -eq 0 ]; then
        echo "All samples completed successfully."
    else
        echo "Completed with failures:"
        printf '  - %s\n' "$@"
        return 1
    fi
}
