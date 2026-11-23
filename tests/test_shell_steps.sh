#!/usr/bin/env bash
# Stubbed end-to-end test of the shell steps (01, 02, 03, 04, 05).
# Real bioinformatics tools are replaced by tiny stubs so this runs in seconds anywhere. It checks:
#   - paths containing spaces work, - finished work is skipped on rerun,
#   - a failure on either side of the HISAT2 | samtools pipe is caught,
#   - strandedness and paired/single-end options reach the tools, - step 5 receives the right paths.
set -uo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
PASS=0; FAIL=0
check() { local name="$1"; shift; if "$@"; then echo "ok   - $name"; PASS=$((PASS + 1)); else echo "FAIL - $name"; FAIL=$((FAIL + 1)); fi; }

DATA="$WORK/My Project"            # deliberately contains a space
mkdir -p "$DATA" "$WORK/bin"
cp -R "$REPO" "$DATA/repo"
R="$DATA/repo"
mkdir -p "$DATA/hisat2_index" "$DATA/annotation"
for i in 1 2 3 4 5 6 7 8; do touch "$DATA/hisat2_index/hg38.$i.ht2"; done
echo x > "$DATA/annotation/genes.gtf"
for s in SRR40723926 SRR40723927 SRR40723928; do for r in 1 2; do echo data | gzip > "$DATA/${s}_${r}.fastq.gz"; done; done

# ---- stub tools; every call is logged to $WORK/calls.log ----
cat > "$WORK/bin/fastp" <<'STUB'
#!/usr/bin/env bash
echo "fastp $*" >> "$(dirname "$0")/../calls.log"
while [ $# -gt 0 ]; do case "$1" in -o|-O|-h|-j) echo stub > "$2"; shift;; esac; shift; done
STUB
cat > "$WORK/bin/hisat2" <<'STUB'
#!/usr/bin/env bash
echo "hisat2 $*" >> "$(dirname "$0")/../calls.log"
args=("$@")
for ((i = 0; i < ${#args[@]}; i++)); do
  if [ "${args[$i]}" = "--summary-file" ]; then
    printf '100 reads; of these:\n  100 (100.00%%) were paired\n95.50%% overall alignment rate\n' > "${args[$((i + 1))]}"
  fi
done
echo "@HD	VN:1.0"
[ "${FAIL_HISAT:-0}" = 1 ] && exit 3
exit 0
STUB
cat > "$WORK/bin/samtools" <<'STUB'
#!/usr/bin/env bash
if [ "$1" = sort ]; then
  out=""; while [ $# -gt 0 ]; do [ "$1" = -o ] && out="$2"; shift; done
  cat > "$out"
  [ "${FAIL_SORT:-0}" = 1 ] && exit 4
elif [ "$1" = index ]; then
  echo idx > "${@: -1}.bai"
fi
exit 0
STUB
cat > "$WORK/bin/featureCounts" <<'STUB'
#!/usr/bin/env bash
if [ $# -eq 0 ]; then echo "usage: featureCounts ... --countReadPairs"; exit 1; fi
echo "featureCounts $*" >> "$(dirname "$0")/../calls.log"
out=""; bams=()
while [ $# -gt 0 ]; do
  case "$1" in -o) out="$2"; shift;; -T|-a|-t|-g|-s) shift;; -p|--countReadPairs) ;; *) bams+=("$1");; esac
  shift
done
{ echo "# Program:featureCounts"; printf 'Geneid\tChr\tStart\tEnd\tStrand\tLength'; for b in "${bams[@]}"; do printf '\t%s' "$b"; done; echo
  printf 'G1\tchr1\t1\t2\t+\t2'; for b in "${bams[@]}"; do printf '\t10'; done; echo; } > "$out"
{ printf 'Status'; for b in "${bams[@]}"; do printf '\t%s' "$b"; done; echo
  printf 'Assigned'; for b in "${bams[@]}"; do printf '\t80'; done; echo
  printf 'Unassigned_NoFeatures'; for b in "${bams[@]}"; do printf '\t20'; done; echo; } > "$out.summary"
STUB
for t in prefetch fasterq-dump; do printf '#!/usr/bin/env bash\nexit 0\n' > "$WORK/bin/$t"; done
cat > "$WORK/bin/Rscript" <<'STUB'
#!/usr/bin/env bash
echo "Rscript $* COUNTS_FILE=$COUNTS_FILE BAM_DIR=$BAM_DIR" >> "$(dirname "$0")/../calls.log"
STUB
chmod +x "$WORK/bin/"*
export PATH="$WORK/bin:$PATH"
LOG="$WORK/calls.log"; : > "$LOG"
cd "$R" || exit 1

# ---------------- config parsing ----------------
n_samples=$(bash -c 'source config.sh; echo "${#SAMPLES[@]}"')
check "config.sh reads 3 samples from the sample sheet" test "$n_samples" = 3
data_dir=$(bash -c 'source config.sh; echo "$DATA_DIR"')
check "DATA_DIR resolves to the folder above the repo (path with a space)" test "$data_dir" = "$DATA"

# ---------------- step 1 ----------------
out1=$(bash scripts/01_download_sra.sh 2>&1)
check "step 1 skips samples whose FASTQ already exist" test "$(grep -c 'already present' <<< "$out1")" = 3

# ---------------- step 2 ----------------
bash scripts/02_trim_fastp.sh > /dev/null 2>&1
check "step 2 writes paired clean FASTQ and reports" test -s "$DATA/clean_fastq/SRR40723926_2.clean.fastq.gz" -a -s "$DATA/reports/SRR40723926_fastp.html"
out2=$(bash scripts/02_trim_fastp.sh 2>&1)
check "step 2 skips on rerun" test "$(grep -c 'Already trimmed' <<< "$out2")" = 3
check "step 2 passes --detect_adapter_for_pe for paired data" grep -q -- '--detect_adapter_for_pe' "$LOG"

# ---------------- step 3 ----------------
bash scripts/03_align_hisat2.sh > /dev/null 2>&1
check "step 3 writes sorted BAM, index and summary" test -s "$DATA/aligned_bam/SRR40723927.sorted.bam" -a -s "$DATA/aligned_bam/SRR40723927.sorted.bam.bai" -a -s "$DATA/aligned_bam/SRR40723927_summary.txt"
out3=$(bash scripts/03_align_hisat2.sh 2>&1)
check "step 3 skips finished samples on rerun" test "$(grep -c 'skipping' <<< "$out3")" = 3

rm -rf "$DATA/aligned_bam"
FAIL_HISAT=1 bash scripts/03_align_hisat2.sh > /dev/null 2>&1; rc=$?
check "step 3 exits non-zero when hisat2 fails" test "$rc" -ne 0
check "step 3 leaves no BAM behind after a hisat2 failure" test "$(find "$DATA/aligned_bam" -name '*.sorted.bam' | wc -l)" -eq 0
rm -rf "$DATA/aligned_bam"
FAIL_SORT=1 bash scripts/03_align_hisat2.sh > /dev/null 2>&1; rc=$?
check "step 3 exits non-zero when samtools sort fails" test "$rc" -ne 0

# strandedness and single-end plumbing (edit the copy of config.sh used by this test, then restore it)
cp config.sh "$WORK/config.sh.orig"
rm -rf "$DATA/aligned_bam"; : > "$LOG"
sed -i.bak 's/^STRANDEDNESS=0/STRANDEDNESS=2/' config.sh
bash scripts/03_align_hisat2.sh > /dev/null 2>&1
check "STRANDEDNESS=2 becomes --rna-strandness RF for paired data" grep -q -- '--rna-strandness RF' "$LOG"
cp "$WORK/config.sh.orig" config.sh

rm -rf "$DATA/aligned_bam" "$DATA/clean_fastq"; : > "$LOG"
sed -i.bak 's/^PAIRED=true/PAIRED=false/' config.sh
bash scripts/02_trim_fastp.sh > /dev/null 2>&1
bash scripts/03_align_hisat2.sh > /dev/null 2>&1
check "single-end: step 2 writes only read 1" test -s "$DATA/clean_fastq/SRR40723926_1.clean.fastq.gz" -a ! -e "$DATA/clean_fastq/SRR40723926_2.clean.fastq.gz"
check "single-end: step 3 uses -U" grep -q 'hisat2 .* -U ' "$LOG"
cp "$WORK/config.sh.orig" config.sh; rm -f config.sh.bak

# rebuild the paired-end inputs for the remaining steps
rm -rf "$DATA/aligned_bam" "$DATA/clean_fastq"
bash scripts/02_trim_fastp.sh > /dev/null 2>&1
bash scripts/03_align_hisat2.sh > /dev/null 2>&1

# ---------------- step 4 ----------------
: > "$LOG"
out4=$(bash scripts/04_count_featurecounts.sh 2>&1)
check "step 4 writes the counts table and its summary" test -s "$DATA/gene_counts.txt" -a -s "$DATA/gene_counts.txt.summary"
check "step 4 prints a per-sample assigned percentage" test "$(grep -c '80.0%' <<< "$out4")" = 3
check "step 4 adds --countReadPairs for paired data" grep -q -- '--countReadPairs' "$LOG"
out4b=$(bash scripts/04_count_featurecounts.sh 2>&1)
check "step 4 skips when the counts table exists" grep -q 'skipping' <<< "$out4b"

# ---------------- step 5 ----------------
: > "$LOG"
bash scripts/05_run_deseq2.sh > /dev/null 2>&1
check "step 5 runs deseq2.R with COUNTS_FILE and BAM_DIR set" grep -q "Rscript scripts/deseq2.R COUNTS_FILE=$DATA/gene_counts.txt BAM_DIR=$DATA/aligned_bam" "$LOG"
rm -f "$DATA/gene_counts.txt"
bash scripts/05_run_deseq2.sh > /dev/null 2>&1; rc=$?
check "step 5 refuses to run without a counts table" test "$rc" -ne 0

echo
echo "passed: $PASS   failed: $FAIL"
[ "$FAIL" -eq 0 ]
