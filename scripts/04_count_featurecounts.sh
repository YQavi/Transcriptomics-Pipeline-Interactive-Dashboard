#!/usr/bin/env bash
# Step 4: gene-level counts for all samples with featureCounts (one table, one column per BAM).
# Skips if the counts table already exists; delete it (or set FORCE=1) to recount.
set -uo pipefail
# shellcheck source=../config.sh
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/config.sh"
need_cmds featureCounts

if [ -s "$COUNTS_FILE" ] && [ "${FORCE:-0}" != "1" ]; then
    echo "Counts table exists ($COUNTS_FILE), skipping. Use FORCE=1 to recount."
    exit 0
fi
if [ ! -f "$GTF" ]; then
    echo "ERROR: annotation not found: $GTF" >&2
    echo "       Set GTF in config.sh to a GTF on the SAME genome build and chromosome naming as the HISAT2 index." >&2
    exit 1
fi

BAMS=()
for sample in "${SAMPLES[@]}"; do
    b="$BAM_DIR/${sample}.sorted.bam"
    [ -s "$b" ] || { echo "ERROR: missing BAM $b (run step 03)." >&2; exit 1; }
    BAMS+=("$b")
done

# Paired-end: Subread >= 2.0.2 needs --countReadPairs alongside -p; older versions do not know it.
PAIR_OPT=""
if [ "$PAIRED" = true ]; then
    usage="$(featureCounts 2>&1 || true)"
    if grep -q -- '--countReadPairs' <<< "$usage"; then PAIR_OPT="-p --countReadPairs"; else PAIR_OPT="-p"; fi
fi

mkdir -p "$(dirname "$COUNTS_FILE")"
echo "featureCounts on ${#BAMS[@]} BAMs (strandedness -s $STRANDEDNESS) at $(date)"
# PAIR_OPT is unquoted on purpose so its flags split into separate words
# shellcheck disable=SC2086
if ! featureCounts -T "$THREADS_FEATURECOUNTS" -a "$GTF" -t exon -g gene_id \
        -s "$STRANDEDNESS" $PAIR_OPT -o "${COUNTS_FILE}.part" "${BAMS[@]}"; then
    echo "ERROR: featureCounts failed." >&2; rm -f "${COUNTS_FILE}.part" "${COUNTS_FILE}.part.summary"; exit 1
fi
mv "${COUNTS_FILE}.part" "$COUNTS_FILE"
mv "${COUNTS_FILE}.part.summary" "${COUNTS_FILE}.summary"

# Assignment rate per sample: a low value usually means the wrong -s, or a GTF/genome mismatch
echo "Assigned-read percentage per sample (aim for well above 50%):"
awk -F'\t' 'NR==1 {n=NF; for (i=2;i<=NF;i++) name[i]=$i; next}
            {for (i=2;i<=NF;i++) {tot[i]+=$i; if ($1=="Assigned") asg[i]=$i}}
            END {for (i=2;i<=n;i++) if (tot[i]>0) printf "  %s  %.1f%%\n", name[i], 100*asg[i]/tot[i]}' \
    "${COUNTS_FILE}.summary"
