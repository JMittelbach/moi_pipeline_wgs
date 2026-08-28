#!/usr/bin/env bash
set -Eeuo pipefail

CONFIG="${1:-}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck disable=SC1090
source "$SCRIPT_DIR/lib.sh" "$CONFIG"

progress "6/9" "extracting REF/ALT allele counts at frozen targets"
require_command bcftools
require_command tabix
mkdir -p "$OUTPUT_DIR/counts"

sample_total=0
while IFS=$'\t' read -r sample _r1 _r2; do
  [[ -z "${sample:-}" ]] && continue
  bam="$(dedup_bam "$sample")"
  output="$(counts_tsv "$sample")"
  log="$OUTPUT_DIR/logs/06_counts_${sample}.log"
  sample_total=$((sample_total + 1))
  if [[ "$RESUME" == 1 && -s "$output" ]]; then
    echo "[6/9] [$sample_total] $sample: counts already exist (resume)"
    continue
  fi
  echo "[6/9] [$sample_total] $sample: bcftools mpileup/call"
  temporary_dir="$(mktemp -d "$OUTPUT_DIR/counts/.tmp_${sample}.XXXXXX")"
  temporary_bcf="$temporary_dir/counts.bcf"
  temporary_tsv="$temporary_dir/counts.tsv"
  cleanup() { rm -rf "$temporary_dir"; }
  if {
    set -o pipefail
    bcftools mpileup -Ou -f "$REFERENCE_LINK" -T "$TARGETS" \
      -q "$MIN_MAPQ" -Q "$MIN_BASEQ" -d "$MAX_DEPTH" -I \
      -a FORMAT/AD,FORMAT/DP "$bam" \
      | bcftools call -m -A -i --ploidy 1 -C alleles -T "$TARGETS" \
          -Ob -o "$temporary_bcf"
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD\t%DP]\n' \
      "$temporary_bcf" > "$temporary_tsv"
    mv "$temporary_tsv" "$output"
  } >"$log" 2>&1; then
    cleanup
  else
    cleanup
    echo "[ERROR] allele-count extraction failed for $sample. Last log lines:" >&2
    tail -n 40 "$log" >&2 || true
    fail "allele-count extraction failed for $sample" \
      "inspect $log; verify TARGETS is indexed and matches the reference contig names"
  fi
  [[ -s "$output" ]] || fail "empty allele-count table for $sample" \
    "inspect $log and verify BAM coverage at the frozen targets"
done < <(sample_rows)
echo "[6/9] OK: allele-count tables completed for $sample_total samples"
