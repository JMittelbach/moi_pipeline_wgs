#!/usr/bin/env bash
set -Eeuo pipefail

CONFIG="${1:-}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck disable=SC1090
source "$SCRIPT_DIR/lib.sh" "$CONFIG"

if [[ "$RUN_FREEBAYES" != 1 ]]; then
  echo "[6a/9] FreeBayes variant calling disabled (set RUN_FREEBAYES=1 to enable)"
  exit 0
fi

progress "6a/9" "calling variants with FreeBayes"
require_command freebayes
require_command bgzip
require_command tabix
mkdir -p "$OUTPUT_DIR/variants"

sample_total=0
while IFS=$'\t' read -r sample _r1 _r2; do
  [[ -z "${sample:-}" ]] && continue
  bam="$(dedup_bam "$sample")"
  output="$(freebayes_vcf "$sample")"
  index="$output.tbi"
  log="$PROCESSED_DIR/logs/06_freebayes_${sample}.log"
  sample_total=$((sample_total + 1))
  if [[ "$RESUME" == 1 && -s "$output" && -s "$index" ]]; then
    echo "[6a/9] [$sample_total] $sample: FreeBayes VCF already exists (resume)"
    continue
  fi
  [[ -s "$bam" ]] || fail "deduplicated BAM is missing for $sample: $bam" \
    "complete step 05 before running FreeBayes"
  echo "[6a/9] [$sample_total] $sample: FreeBayes"
  temporary_dir="$(mktemp -d "$PROCESSED_DIR/variants.XXXXXX")"
  temporary_vcf="$temporary_dir/$sample.freebayes.vcf.gz"
  cleanup() { rm -rf "$temporary_dir"; }
  if {
    set -o pipefail
    freebayes \
      --fasta-reference "$REFERENCE_LINK" \
      --ploidy "$FREEBAYES_PLOIDY" \
      --min-mapping-quality "$MIN_MAPQ" \
      --min-base-quality "$MIN_BASEQ" \
      --min-alternate-count "$FREEBAYES_MIN_ALT_COUNT" \
      --min-alternate-fraction "$FREEBAYES_MIN_ALT_FRACTION" \
      "$bam" | bgzip -c > "$temporary_vcf"
    tabix -f -p vcf "$temporary_vcf"
    mkdir -p "$(dirname "$output")"
    mv "$temporary_vcf" "$output"
    mv "$temporary_vcf.tbi" "$index"
  } >"$log" 2>&1; then
    cleanup
  else
    cleanup
    echo "[ERROR] FreeBayes variant calling failed for $sample. Last log lines:" >&2
    tail -n 40 "$log" >&2 || true
    fail "FreeBayes variant calling failed for $sample" \
      "inspect $log and verify the reference/BAM and FreeBayes thresholds"
  fi
  [[ -s "$output" && -s "$index" ]] || fail "FreeBayes output is incomplete for $sample" \
    "inspect $log and check available disk space"
done < <(sample_rows)
echo "[6a/9] OK: FreeBayes VCFs completed for $sample_total samples"
