#!/usr/bin/env bash
set -Eeuo pipefail

CONFIG="${1:-}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck disable=SC1090
source "$SCRIPT_DIR/lib.sh" "$CONFIG"

progress "3/9" "trimming paired reads with fastp"
require_command fastp
mkdir -p "$OUTPUT_DIR/trimmed" "$OUTPUT_DIR/qc/trim"

sample_total=0
while IFS=$'\t' read -r sample r1 r2; do
  [[ -z "${sample:-}" ]] && continue
  in1="$(sample_path "$r1")"
  in2="$(sample_path "$r2")"
  out1="$(trim_r1 "$sample")"
  out2="$(trim_r2 "$sample")"
  json="$OUTPUT_DIR/qc/trim/${sample}.json"
  html="$OUTPUT_DIR/qc/trim/${sample}.html"
  log="$OUTPUT_DIR/logs/03_trim_${sample}.log"
  sample_total=$((sample_total + 1))
  if [[ "$RESUME" == 1 && -s "$out1" && -s "$out2" ]]; then
    echo "[3/9] [$sample_total] $sample: trimmed files already exist (resume)"
    continue
  fi
  echo "[3/9] [$sample_total] $sample: trimming"
  fastp_command=(fastp --in1 "$in1" --in2 "$in2" --out1 "$out1" --out2 "$out2"
    --thread "$THREADS" --qualified_quality_phred "$TRIM_QUALITY"
    --length_required "$MIN_READ_LENGTH" --detect_adapter_for_pe
    --json "$json" --html "$html")
  [[ -n "${ADAPTER_R1:-}" ]] && fastp_command+=(--adapter_sequence "$ADAPTER_R1")
  [[ -n "${ADAPTER_R2:-}" ]] && fastp_command+=(--adapter_sequence_r2 "$ADAPTER_R2")
  run_logged "$log" "fastp for $sample" "${fastp_command[@]}"
  [[ -s "$out1" && -s "$out2" ]] || fail "fastp produced no reads for $sample" \
    "inspect $log; check input compression and read pairing"
done < <(sample_rows)
echo "[3/9] OK: trimming completed for $sample_total samples"
