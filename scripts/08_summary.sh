#!/usr/bin/env bash
set -Eeuo pipefail

CONFIG="${1:-}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck disable=SC1090
source "$SCRIPT_DIR/lib.sh" "$CONFIG"

progress "8/9" "writing one compact MOI/Fws summary"
summary="$OUTPUT_DIR/moi_fws_summary.tsv"
per_sample="$OUTPUT_DIR/moi_per_sample.tsv"
first=1
for metrics in "$OUTPUT_DIR"/metrics/*.moi_fws.tsv; do
  [[ -e "$metrics" ]] || continue
  if [[ "$first" == 1 ]]; then
    cat "$metrics" > "$summary"
    first=0
  else
    tail -n +2 "$metrics" >> "$summary"
  fi
done
[[ "$first" == 0 && -s "$summary" ]] || fail "no metrics files found" "complete step 07 first"
rows=$(( $(wc -l < "$summary") - 1 ))

# Human-readable one-line result per sample.  BIC weight is model-selection
# support, not a calibrated statistical confidence interval (moimix does not
# provide a formal CI in this interface).
temporary="$per_sample.tmp"
printf 'sample_id\tmoi\tmoi_status\tbic\tbic_delta\tbic_weight\tconfidence\tfws\n' > "$temporary"
awk -F '\t' '
  NR == 1 { next }
  $2 == "binommix_moi" { moi[$1]=$5; status[$1]=$3; bic[$1]=$9; delta[$1]=$10; weight[$1]=$11; confidence[$1]=$12 }
  $2 == "fws_moimix_compatible" { fws[$1]=$5 }
  END { for (sample in moi) print sample "\t" moi[sample] "\t" status[sample] "\t" bic[sample] "\t" delta[sample] "\t" weight[sample] "\t" confidence[sample] "\t" fws[sample] }
' "$summary" | sort >> "$temporary"
mv "$temporary" "$per_sample"
sample_rows=$(( $(wc -l < "$per_sample") - 1 ))
echo "[8/9] OK: $rows metric rows for $sample_rows samples"
echo "        long table: $summary"
echo "        per-sample table: $per_sample"
