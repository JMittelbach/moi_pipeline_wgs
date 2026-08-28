#!/usr/bin/env bash
set -Eeuo pipefail

CONFIG="${1:-}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck disable=SC1090
source "$SCRIPT_DIR/lib.sh" "$CONFIG"

if [[ "$RUN_GATK4" != 1 ]]; then
  echo "[6c/9] GATK4 variant calling disabled (set RUN_GATK4=1 to enable)"
  exit 0
fi

progress "6c/9" "calling variants with GATK4 HaplotypeCaller/GenotypeGVCFs"
require_command gatk
mkdir -p "$OUTPUT_DIR/variants" "$PROCESSED_DIR/gatk4"

sample_total=0
while IFS=$'\t' read -r sample _r1 _r2; do
  [[ -z "${sample:-}" ]] && continue
  bam="$(dedup_bam "$sample")"
  gvcf="$(gatk4_gvcf "$sample")"
  gvcf_index="$gvcf.tbi"
  output="$(gatk4_vcf "$sample")"
  index="$output.tbi"
  log="$PROCESSED_DIR/logs/06_gatk4_${sample}.log"
  sample_total=$((sample_total + 1))
  if [[ "$RESUME" == 1 && -s "$output" && -s "$index" ]]; then
    echo "[6c/9] [$sample_total] $sample: GATK4 VCF already exists (resume)"
    continue
  fi
  [[ -s "$bam" ]] || fail "deduplicated BAM is missing for $sample: $bam" \
    "complete step 05 before running GATK4"
  echo "[6c/9] [$sample_total] $sample: HaplotypeCaller + GenotypeGVCFs"
  temporary_dir="$(mktemp -d "$PROCESSED_DIR/gatk4/.tmp_${sample}.XXXXXX")"
  temporary_gvcf="$temporary_dir/$sample.g.vcf.gz"
  temporary_vcf="$temporary_dir/$sample.vcf.gz"
  cleanup() { rm -rf "$temporary_dir"; }
  run_gatk_for_sample() {
    haplotype_command=(
      gatk --java-options "$GATK_JAVA_OPTIONS" HaplotypeCaller
      -R "$REFERENCE_LINK"
      -I "$bam"
      -O "$temporary_gvcf"
      -ERC GVCF
      --sample-ploidy "$GATK_PLOIDY"
      --heterozygosity "$GATK_HETEROZYGOSITY"
      --indel-heterozygosity "$GATK_INDEL_HETEROZYGOSITY"
      --min-assembly-region-size "$GATK_MIN_ASSEMBLY_REGION_SIZE"
      --min-base-quality-score "$GATK_MIN_BASE_QUALITY_SCORE"
      --base-quality-score-threshold "$GATK_BASE_QUALITY_SCORE_THRESHOLD"
    )
    if [[ "$GATK_DISABLE_MAPPING_QUALITY_FILTER" == 1 ]]; then
      haplotype_command+=(--disable-read-filter MappingQualityReadFilter)
    fi
    "${haplotype_command[@]}" || return 1
    gatk --java-options "$GATK_JAVA_OPTIONS" IndexFeatureFile -I "$temporary_gvcf" || return 1
    genotype_command=(
      gatk --java-options "$GATK_JAVA_OPTIONS" GenotypeGVCFs
      -R "$REFERENCE_LINK"
      -V "$temporary_gvcf"
      -O "$temporary_vcf"
      --sample-ploidy "$GATK_PLOIDY"
      --stand-call-conf "$GATK_STAND_CALL_CONF"
    )
    "${genotype_command[@]}" || return 1
    gatk --java-options "$GATK_JAVA_OPTIONS" IndexFeatureFile -I "$temporary_vcf" || return 1
    mkdir -p "$(dirname "$gvcf")" "$(dirname "$output")" || return 1
    mv "$temporary_gvcf" "$gvcf" || return 1
    mv "$temporary_gvcf.tbi" "$gvcf_index" || return 1
    mv "$temporary_vcf" "$output" || return 1
    mv "$temporary_vcf.tbi" "$index" || return 1
  }
  if run_gatk_for_sample >"$log" 2>&1; then
    cleanup
  else
    cleanup
    echo "[ERROR] GATK4 variant calling failed for $sample. Last log lines:" >&2
    tail -n 40 "$log" >&2 || true
    fail "GATK4 variant calling failed for $sample" \
      "inspect $log and verify the reference dictionary, BAM, and GATK parameters"
  fi
  [[ -s "$output" && -s "$index" ]] || fail "GATK4 output is incomplete for $sample" \
    "inspect $log and check available disk space"
done < <(sample_rows)
echo "[6c/9] OK: GATK4 VCFs completed for $sample_total samples"
