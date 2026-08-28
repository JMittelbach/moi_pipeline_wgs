#!/usr/bin/env bash
set -Eeuo pipefail

# Shared, deliberately small helpers for the numbered pipeline steps.
PIPELINE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONFIG_FILE="${1:-$PIPELINE_ROOT/config/pipeline.env}"

fail() {
  local message="$1"
  local fix="${2:-}"
  echo >&2
  echo "[ERROR] $message" >&2
  [[ -n "$fix" ]] && echo "        Fix: $fix" >&2
  if [[ -n "${PROCESSED_DIR:-}" ]]; then
    echo "        Logs: ${PROCESSED_DIR}/logs" >&2
  fi
  exit 1
}

[[ -f "$CONFIG_FILE" ]] || fail \
  "configuration file not found: $CONFIG_FILE" \
  "copy config/pipeline.env and edit its paths"

# The configuration is a local, user-controlled shell file.  Values containing
# spaces must be quoted in it.
# shellcheck disable=SC1090
source "$CONFIG_FILE"

resolve_path() {
  local value="$1"
  if [[ "$value" == /* ]]; then
    printf '%s\n' "$value"
  else
    printf '%s\n' "$PIPELINE_ROOT/$value"
  fi
}

required_setting() {
  local name="$1"
  [[ -n "${!name:-}" ]] || fail \
    "required setting $name is empty" \
    "set $name in $CONFIG_FILE"
}

for setting in SAMPLES_TSV RAW_DIR REFERENCE TARGETS POPULATION_PANEL PROCESSED_DIR OUTPUT_DIR; do
  required_setting "$setting"
done

SAMPLES_TSV="$(resolve_path "$SAMPLES_TSV")"
RAW_DIR="$(resolve_path "$RAW_DIR")"
REFERENCE="$(resolve_path "$REFERENCE")"
TARGETS="$(resolve_path "$TARGETS")"
POPULATION_PANEL="$(resolve_path "$POPULATION_PANEL")"
PROCESSED_DIR="$(resolve_path "$PROCESSED_DIR")"
OUTPUT_DIR="$(resolve_path "$OUTPUT_DIR")"
PLOTS_DIR="${PLOTS_DIR:-plots}"
if [[ "$PLOTS_DIR" == /* ]]; then
  PLOTS_DIR="$PLOTS_DIR"
else
  PLOTS_DIR="$OUTPUT_DIR/$PLOTS_DIR"
fi
REFERENCE_LINK="$PROCESSED_DIR/reference/reference.fasta"
BWA_PREFIX="$PROCESSED_DIR/reference/pf3d7"

THREADS="${THREADS:-8}"
TRIM_QUALITY="${TRIM_QUALITY:-30}"
MIN_READ_LENGTH="${MIN_READ_LENGTH:-75}"
MIN_MAPQ="${MIN_MAPQ:-30}"
MIN_BASEQ="${MIN_BASEQ:-20}"
MAX_DEPTH="${MAX_DEPTH:-100000}"
POPULATION="${POPULATION:-Global}"
FWS_MIN_DEPTH="${FWS_MIN_DEPTH:-50}"
MAX_UNMODELLED_FRACTION="${MAX_UNMODELLED_FRACTION:-0.02}"
K_VALUES="${K_VALUES:-1,2,3,4,5}"
MOI_COVERAGE_THRESHOLD="${MOI_COVERAGE_THRESHOLD:-10}"
MOI_NITER="${MOI_NITER:-1000}"
RESUME="${RESUME:-1}"

[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || fail "THREADS must be a positive integer" "set THREADS in $CONFIG_FILE"
[[ "$RESUME" =~ ^[01]$ ]] || fail "RESUME must be 0 or 1" "set RESUME=1 or RESUME=0"

mkdir -p "$PROCESSED_DIR/logs"

progress() {
  echo "[$1] ${*:2}"
}

require_command() {
  local command_name="$1"
  command -v "$command_name" >/dev/null 2>&1 || fail \
    "required program not found: $command_name" \
    "run ./setup.sh, then activate the environment with: conda activate moi_pipeline"
}

run_logged() {
  local logfile="$1"
  local description="$2"
  shift 2
  mkdir -p "$(dirname "$logfile")"
  if ! "$@" >"$logfile" 2>&1; then
    echo "[ERROR] $description failed. Last log lines:" >&2
    tail -n 40 "$logfile" >&2 || true
    fail "$description failed" "inspect $logfile and fix the reported input/tool problem"
  fi
}

sample_rows() {
  tail -n +2 "$SAMPLES_TSV" | tr -d '\r'
}

sample_path() {
  local value="$1"
  if [[ "$value" == /* ]]; then
    printf '%s\n' "$value"
  else
    printf '%s\n' "$RAW_DIR/$value"
  fi
}

trim_r1() { printf '%s\n' "$PROCESSED_DIR/trimmed/$1_R1.fastq.gz"; }
trim_r2() { printf '%s\n' "$PROCESSED_DIR/trimmed/$1_R2.fastq.gz"; }
name_bam() { printf '%s\n' "$PROCESSED_DIR/bam/$1.name.bam"; }
dedup_bam() { printf '%s\n' "$PROCESSED_DIR/bam/$1.dedup.bam"; }
counts_tsv() { printf '%s\n' "$PROCESSED_DIR/counts/$1.tsv"; }
metrics_tsv() { printf '%s\n' "$OUTPUT_DIR/metrics/$1.moi_fws.tsv"; }
