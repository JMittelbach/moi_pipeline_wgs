#!/usr/bin/env bash
set -Eeuo pipefail

PIPELINE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

usage() {
  cat <<'EOF'
Usage: ./run_pipeline.sh build-metadata
       ./run_pipeline.sh [CONFIG]
       ./run_pipeline.sh --help

build-metadata scans data/ for Illumina FASTQ/FASTQ.GZ files, writes the
reviewable data/metadata.txt sample/lane/read manifest, and then stops.
Inspect or correct that file manually before running ./run_pipeline.sh again.
With no argument, the normal pipeline uses config/pipeline.env.
EOF
}

if [[ "${1:-}" == "--help" || "${1:-}" == "-h" ]]; then
  usage
  exit 0
fi

if [[ "${1:-}" == "build-metadata" ]]; then
  [[ $# -eq 1 ]] || { echo "[ERROR] build-metadata does not accept extra arguments" >&2; usage >&2; exit 2; }
  "$PIPELINE_ROOT/scripts/build_metadata.py" "$PIPELINE_ROOT/data" "$PIPELINE_ROOT/data/metadata.txt"
  echo
  echo "Bitte data/metadata.txt jetzt öffnen und Sample-, Lane- und Read-Zuordnung prüfen."
  if [[ -t 0 && -t 1 ]]; then
    read -r -p "Hast du die Datei geprüft? [j/N] " answer || answer=""
    case "$answer" in
      j|J|y|Y) echo "Danke. Die eigentliche Pipeline startet erst mit dem nächsten Aufruf: ./run_pipeline.sh" ;;
      *) echo "Bitte korrigieren und danach erneut ./run_pipeline.sh aufrufen." ;;
    esac
  else
    echo "Nicht-interaktive Sitzung: prüfen/korrigieren und danach erneut ./run_pipeline.sh aufrufen."
  fi
  exit 0
fi

if [[ $# -gt 1 ]]; then
  echo "[ERROR] too many arguments" >&2
  usage >&2
  exit 2
fi

CONFIG="${1:-$PIPELINE_ROOT/config/pipeline.env}"
if [[ "$CONFIG" != /* ]]; then
  CONFIG="$PIPELINE_ROOT/$CONFIG"
fi
[[ -f "$CONFIG" ]] || { echo "[ERROR] configuration file not found: $CONFIG" >&2; exit 1; }

"$PIPELINE_ROOT/scripts/metadata_to_samples.py" \
  "$PIPELINE_ROOT/data/metadata.txt" \
  "$PIPELINE_ROOT/data" \
  "$PIPELINE_ROOT/data/samples.tsv"
metadata_state="$PIPELINE_ROOT/data/.metadata.sha256"
current_metadata_hash="$(<"$metadata_state")"
previous_metadata_hash=""
if [[ -f "$metadata_state.previous" ]]; then
  previous_metadata_hash="$(<"$metadata_state.previous")"
fi
force_recompute=0
if [[ "$previous_metadata_hash" != "$current_metadata_hash" ]]; then
  force_recompute=1
fi
cp "$metadata_state" "$metadata_state.previous"
if [[ "$force_recompute" == 1 ]]; then
  echo "Metadata changed or was prepared for the first time; recomputing pipeline outputs (RESUME=0)."
fi
FORCE_RECOMPUTE="$force_recompute" exec "$PIPELINE_ROOT/scripts/00_run_pipeline.sh" "$CONFIG"
