#!/usr/bin/env bash
set -Eeuo pipefail

CONFIG="${1:-}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck disable=SC1090
source "$SCRIPT_DIR/lib.sh" "$CONFIG"

progress "1/9" "checking programs, configuration, and input files"
for program in python Rscript fastp bwa samtools bcftools bgzip tabix; do
  require_command "$program"
done
if [[ "$RUN_FREEBAYES" == 1 ]]; then
  require_command freebayes
fi

[[ -r "$SAMPLES_TSV" ]] || fail "sample sheet is not readable: $SAMPLES_TSV" \
  "create a tab-separated file with header sample_id<TAB>R1<TAB>R2"
[[ -s "$REFERENCE" ]] || fail "reference FASTA is missing or empty: $REFERENCE" \
  "set REFERENCE to the Plasmodium FASTA in config/pipeline.env"
[[ -s "$TARGETS" ]] || fail "target file is missing or empty: $TARGETS" \
  "provide the frozen REF/ALT target file used for the population panel"
[[ -s "$POPULATION_PANEL" ]] || fail "population panel is missing or empty: $POPULATION_PANEL" \
  "provide a target-aligned panel containing a Global column"

if [[ "$TARGETS" == *.gz && ! -s "$TARGETS.tbi" && ! -s "$TARGETS.csi" ]]; then
  fail "compressed target file has no tabix/CSI index: $TARGETS" \
    "run tabix -s 1 -b 2 -e 2 $TARGETS or point TARGETS at an indexed file"
fi

header="$(head -n 1 "$SAMPLES_TSV" | tr -d '\r')"
[[ "$header" == $'sample_id\tR1\tR2' ]] || fail \
  "sample sheet header must be exactly: sample_id<TAB>R1<TAB>R2" \
  "fix the first line of $SAMPLES_TSV"

awk -F '\t' '
  NR > 1 && $1 != "" {
    if (NF != 3) { print "line " NR ": expected 3 tab-separated fields" > "/dev/stderr"; bad=1 }
    if (++seen[$1] > 1) { print "duplicate sample_id: " $1 > "/dev/stderr"; bad=1 }
    n++
  }
  END { if (n == 0) { print "sample sheet has no samples" > "/dev/stderr"; bad=1 } exit bad }
' "$SAMPLES_TSV" || fail "invalid sample sheet: $SAMPLES_TSV" \
  "remove blank/duplicate rows and keep sample_id, R1, R2 columns"

sample_total=0
while IFS=$'\t' read -r sample r1 r2; do
  [[ -z "${sample:-}" ]] && continue
  [[ "$sample" =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*$ ]] || fail \
    "unsafe sample_id '$sample'" "use letters, numbers, '.', '_' or '-' only"
  r1_path="$(sample_path "$r1")"
  r2_path="$(sample_path "$r2")"
  [[ -s "$r1_path" ]] || fail "R1 FASTQ missing for $sample: $r1_path" \
    "fix the R1 path in $SAMPLES_TSV or RAW_DIR in config/pipeline.env"
  [[ -s "$r2_path" ]] || fail "R2 FASTQ missing for $sample: $r2_path" \
    "fix the R2 path in $SAMPLES_TSV or RAW_DIR in config/pipeline.env"
  sample_total=$((sample_total + 1))
done < <(sample_rows)

R_SCRIPT="$(tool_path Rscript)"
if ! "$R_SCRIPT" -e 'ok <- requireNamespace("moimix", quietly=TRUE) || requireNamespace("flexmix", quietly=TRUE); quit(status=if (ok) 0 else 1)' >/dev/null 2>&1; then
  fail "neither the R package moimix nor flexmix is installed" \
    "run ./setup.sh again; flexmix is an explicit fallback when moimix cannot build"
fi

echo "[1/9] OK: $sample_total paired samples; reference, targets, panel, and MOI package found"
