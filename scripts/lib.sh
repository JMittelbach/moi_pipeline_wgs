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

# On macOS, shell activation may append Conda after Homebrew in PATH. Put the
# active environment first so all pipeline tools (not only Rscript) agree.
if [[ -n "${CONDA_PREFIX:-}" && -d "$CONDA_PREFIX/bin" ]]; then
  PATH="$CONDA_PREFIX/bin:$PATH"
  export PATH
fi

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
REFERENCE_DICT="$PROCESSED_DIR/reference/reference.dict"
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
  tool_path "$command_name" >/dev/null || fail \
    "required program not found: $command_name" \
    "run ./setup.sh, then activate the environment with: conda activate moi_pipeline"
}

# Conda activation can append its bin directory after Homebrew's bin on macOS.
# Prefer the active environment's executable explicitly so host R/Python tools
# cannot accidentally be mixed with Conda libraries.
tool_path() {
  local command_name="$1"
  if [[ -n "${CONDA_PREFIX:-}" && -x "$CONDA_PREFIX/bin/$command_name" ]]; then
    printf '%s\n' "$CONDA_PREFIX/bin/$command_name"
  else
    command -v "$command_name" 2>/dev/null || return 1
  fi
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
freebayes_vcf() { printf '%s\n' "$OUTPUT_DIR/variants/$1.freebayes.vcf.gz"; }
gatk4_gvcf() { printf '%s\n' "$PROCESSED_DIR/gatk4/$1.g.vcf.gz"; }
gatk4_vcf() { printf '%s\n' "$OUTPUT_DIR/variants/$1.gatk4.vcf.gz"; }

# FreeBayes is an optional whole-genome variant-calling branch.  The defaults
# keep the historical fixed-target MOI/Fws path unchanged when it is disabled.
RUN_FREEBAYES="${RUN_FREEBAYES:-0}"
FREEBAYES_PLOIDY="${FREEBAYES_PLOIDY:-1}"
FREEBAYES_MIN_ALT_COUNT="${FREEBAYES_MIN_ALT_COUNT:-2}"
FREEBAYES_MIN_ALT_FRACTION="${FREEBAYES_MIN_ALT_FRACTION:-0.2}"
[[ "$RUN_FREEBAYES" =~ ^[01]$ ]] || fail "RUN_FREEBAYES must be 0 or 1" "set RUN_FREEBAYES=0 or RUN_FREEBAYES=1 in $CONFIG_FILE"
[[ "$FREEBAYES_PLOIDY" =~ ^[1-9][0-9]*$ ]] || fail "FREEBAYES_PLOIDY must be a positive integer" "set FREEBAYES_PLOIDY in $CONFIG_FILE"
[[ "$FREEBAYES_MIN_ALT_COUNT" =~ ^[1-9][0-9]*$ ]] || fail "FREEBAYES_MIN_ALT_COUNT must be a positive integer" "set FREEBAYES_MIN_ALT_COUNT in $CONFIG_FILE"
awk -v fraction="$FREEBAYES_MIN_ALT_FRACTION" 'BEGIN { exit !(fraction >= 0 && fraction <= 1) }' || fail \
  "FREEBAYES_MIN_ALT_FRACTION must be between 0 and 1" "set FREEBAYES_MIN_ALT_FRACTION in $CONFIG_FILE"

# GATK4 is an optional HaplotypeCaller -> gVCF -> GenotypeGVCFs branch.  The
# defaults mirror the optimized Pf WGS settings described in the linked paper.
RUN_GATK4="${RUN_GATK4:-0}"
[[ "$RUN_GATK4" =~ ^[01]$ ]] || fail "RUN_GATK4 must be 0 or 1" "set RUN_GATK4=0 or RUN_GATK4=1 in $CONFIG_FILE"
GATK_PLOIDY="${GATK_PLOIDY:-2}"
GATK_HETEROZYGOSITY="${GATK_HETEROZYGOSITY:-0.0029}"
GATK_INDEL_HETEROZYGOSITY="${GATK_INDEL_HETEROZYGOSITY:-0.0017}"
GATK_MIN_ASSEMBLY_REGION_SIZE="${GATK_MIN_ASSEMBLY_REGION_SIZE:-100}"
GATK_MIN_BASE_QUALITY_SCORE="${GATK_MIN_BASE_QUALITY_SCORE:-5}"
GATK_BASE_QUALITY_SCORE_THRESHOLD="${GATK_BASE_QUALITY_SCORE_THRESHOLD:-12}"
GATK_STAND_CALL_CONF="${GATK_STAND_CALL_CONF:-30}"
GATK_JAVA_OPTIONS="${GATK_JAVA_OPTIONS:--Xmx4g}"
GATK_DISABLE_MAPPING_QUALITY_FILTER="${GATK_DISABLE_MAPPING_QUALITY_FILTER:-1}"
[[ "$GATK_PLOIDY" =~ ^[1-9][0-9]*$ ]] || fail "GATK_PLOIDY must be a positive integer" "set GATK_PLOIDY in $CONFIG_FILE"
[[ "$GATK_MIN_ASSEMBLY_REGION_SIZE" =~ ^[1-9][0-9]*$ ]] || fail "GATK_MIN_ASSEMBLY_REGION_SIZE must be a positive integer" "set GATK_MIN_ASSEMBLY_REGION_SIZE in $CONFIG_FILE"
[[ "$GATK_MIN_BASE_QUALITY_SCORE" =~ ^[0-9]+$ ]] || fail "GATK_MIN_BASE_QUALITY_SCORE must be a non-negative integer" "set GATK_MIN_BASE_QUALITY_SCORE in $CONFIG_FILE"
[[ "$GATK_BASE_QUALITY_SCORE_THRESHOLD" =~ ^[0-9]+$ ]] || fail "GATK_BASE_QUALITY_SCORE_THRESHOLD must be a non-negative integer" "set GATK_BASE_QUALITY_SCORE_THRESHOLD in $CONFIG_FILE"
[[ "$GATK_DISABLE_MAPPING_QUALITY_FILTER" =~ ^[01]$ ]] || fail "GATK_DISABLE_MAPPING_QUALITY_FILTER must be 0 or 1" "set GATK_DISABLE_MAPPING_QUALITY_FILTER=0 or 1 in $CONFIG_FILE"
awk -v value="$GATK_HETEROZYGOSITY" 'BEGIN { exit !(value > 0 && value < 1) }' || fail \
  "GATK_HETEROZYGOSITY must be between 0 and 1" "set GATK_HETEROZYGOSITY in $CONFIG_FILE"
awk -v value="$GATK_INDEL_HETEROZYGOSITY" 'BEGIN { exit !(value > 0 && value < 1) }' || fail \
  "GATK_INDEL_HETEROZYGOSITY must be between 0 and 1" "set GATK_INDEL_HETEROZYGOSITY in $CONFIG_FILE"
awk -v value="$GATK_STAND_CALL_CONF" 'BEGIN { exit !(value >= 0) }' || fail \
  "GATK_STAND_CALL_CONF must be non-negative" "set GATK_STAND_CALL_CONF in $CONFIG_FILE"
