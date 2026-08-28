#!/usr/bin/env bash
set -Eeuo pipefail

CONFIG="${1:-}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck disable=SC1090
source "$SCRIPT_DIR/lib.sh" "$CONFIG"

progress "5/9" "fixing mate tags, coordinate-sorting, and removing duplicates"
require_command samtools
mkdir -p "$OUTPUT_DIR/bam" "$OUTPUT_DIR/qc/duplicates"

sample_total=0
while IFS=$'\t' read -r sample _r1 _r2; do
  [[ -z "${sample:-}" ]] && continue
  name_sorted="$(name_bam "$sample")"
  fixmate="$OUTPUT_DIR/bam/$sample.fixmate.bam"
  coordinate="$OUTPUT_DIR/bam/$sample.coordinate.bam"
  output="$(dedup_bam "$sample")"
  bai="$output.bai"
  metrics="$OUTPUT_DIR/qc/duplicates/${sample}.flagstat.txt"
  log="$OUTPUT_DIR/logs/05_duplicates_${sample}.log"
  sample_total=$((sample_total + 1))
  if [[ "$RESUME" == 1 && -s "$output" && -s "$bai" ]]; then
    echo "[5/9] [$sample_total] $sample: duplicate-marked BAM already exists (resume)"
    continue
  fi
  echo "[5/9] [$sample_total] $sample: mark duplicates"
  if ! {
    samtools fixmate -m "$name_sorted" "$fixmate"
    samtools sort -@ "$THREADS" -o "$coordinate" "$fixmate"
    samtools markdup -r -@ "$THREADS" "$coordinate" "$output"
    samtools index -@ "$THREADS" "$output"
    samtools flagstat "$output" > "$metrics"
  } >"$log" 2>&1; then
    echo "[ERROR] duplicate handling failed for $sample. Last log lines:" >&2
    tail -n 40 "$log" >&2 || true
    fail "duplicate handling failed for $sample" \
      "inspect $log; input must be a valid name-sorted paired-end BAM"
  fi
  rm -f "$fixmate" "$coordinate"
  [[ -s "$output" && -s "$bai" ]] || fail "duplicate-marked BAM is incomplete for $sample" \
    "inspect $log and check available disk space"
done < <(sample_rows)
echo "[5/9] OK: duplicate removal and BAM indexing completed for $sample_total samples"
