#!/usr/bin/env bash
set -Eeuo pipefail

PIPELINE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

usage() {
  cat <<'EOF'
Usage: ./scripts/00_run_pipeline.sh [CONFIG]
       ./scripts/00_run_pipeline.sh --help

Run the nine-step paired-end Plasmodium WGS MOI/Fws pipeline sequentially.
CONFIG defaults to config/pipeline.env. Completed files are reused when
RESUME=1. Detailed command logs are written below PROCESSED_DIR/logs; the
terminal shows step and sample progress only.

Options:
  -h, --help    show this help and exit without running a step
EOF
}

case "${1:-}" in
  -h|--help ) usage; exit 0 ;;
  "" ) CONFIG="$PIPELINE_ROOT/config/pipeline.env" ;;
  * ) CONFIG="$1" ;;
esac

if [[ $# -gt 1 ]]; then
  echo "[ERROR] too many arguments" >&2
  usage >&2
  exit 2
fi

if [[ "$CONFIG" != /* ]]; then
  CONFIG="$PIPELINE_ROOT/$CONFIG"
fi

[[ -f "$CONFIG" ]] || {
  echo "[ERROR] configuration file not found: $CONFIG" >&2
  echo "        Fix: edit config/pipeline.env or pass its path as the first argument" >&2
  exit 1
}

echo "Plasmodium WGS MOI/Fws pipeline"
echo "Config: $CONFIG"
echo "The nine steps run sequentially; completed files are reused when RESUME=1."

steps=(
  01_validate_inputs.sh
  02_prepare_reference.sh
  03_trim_reads.sh
  04_align_wgs.sh
  05_mark_duplicates.sh
  06_extract_counts.sh
  07_estimate_moi_fws.py
  08_summary.sh
  09_check_main_table_and_plots.py
)
for step in "${steps[@]}"; do
  script="$PIPELINE_ROOT/scripts/$step"
  [[ -x "$script" ]] || {
    echo "[ERROR] pipeline step is missing or not executable: $script" >&2
    echo "        Fix: run chmod +x scripts/*.sh scripts/*.py" >&2
    exit 1
  }
  "$script" "$CONFIG"
done

echo
echo "DONE: intermediates/QC/logs are in PROCESSED_DIR; final tables and plots are in OUTPUT_DIR."
