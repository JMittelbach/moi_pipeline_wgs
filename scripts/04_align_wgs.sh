#!/usr/bin/env bash
set -Eeuo pipefail

CONFIG="${1:-}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck disable=SC1090
source "$SCRIPT_DIR/lib.sh" "$CONFIG"

progress "4/9" "aligning reads to the Plasmodium reference with BWA-MEM"
require_command bwa
require_command samtools
mkdir -p "$PROCESSED_DIR/bam"

sample_total=0
while IFS=$'\t' read -r sample _r1 _r2; do
  [[ -z "${sample:-}" ]] && continue
  input1="$(trim_r1 "$sample")"
  input2="$(trim_r2 "$sample")"
  output="$(name_bam "$sample")"
  log="$PROCESSED_DIR/logs/04_align_${sample}.log"
  sample_total=$((sample_total + 1))
  if [[ "$RESUME" == 1 && -s "$output" ]]; then
    echo "[4/9] [$sample_total] $sample: alignment already exists (resume)"
    continue
  fi
  echo "[4/9] [$sample_total] $sample: BWA-MEM + name sort"
  rg="@RG\\tID:${sample}\\tSM:${sample}\\tPL:ILLUMINA"
  mkdir -p "$(dirname "$output")"
  if ! { bwa mem -t "$THREADS" -R "$rg" "$REFERENCE_LINK" "$input1" "$input2" \
      | samtools sort -n -@ "$THREADS" -o "$output" -; } >"$log" 2>&1; then
    echo "[ERROR] alignment failed for $sample. Last log lines:" >&2
    tail -n 40 "$log" >&2 || true
    fail "alignment failed for $sample" "inspect $log and verify FASTQs/reference"
  fi
  [[ -s "$output" ]] || fail "alignment output is empty for $sample" \
    "inspect $log and check available disk space"
done < <(sample_rows)
echo "[4/9] OK: alignments completed for $sample_total samples"
