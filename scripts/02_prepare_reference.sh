#!/usr/bin/env bash
set -Eeuo pipefail

CONFIG="${1:-}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck disable=SC1090
source "$SCRIPT_DIR/lib.sh" "$CONFIG"

progress "2/9" "preparing reference indexes"
require_command samtools
require_command bwa
if [[ "$RUN_GATK4" == 1 ]]; then
  require_command gatk
fi
mkdir -p "$PROCESSED_DIR/reference"

if [[ -e "$REFERENCE_LINK" && ! -L "$REFERENCE_LINK" ]]; then
  fail "pipeline reference path is occupied by a real file: $REFERENCE_LINK" \
    "move it aside or choose another PROCESSED_DIR"
fi
if [[ ! -e "$REFERENCE_LINK" && ! -L "$REFERENCE_LINK" ]]; then
  ln -s "$REFERENCE" "$REFERENCE_LINK"
fi

if [[ ! -s "$REFERENCE_LINK.fai" ]]; then
  run_logged "$PROCESSED_DIR/logs/02_faidx.log" "samtools faidx" \
    samtools faidx "$REFERENCE_LINK"
fi

if [[ ! -s "$BWA_PREFIX.bwt" ]]; then
  run_logged "$PROCESSED_DIR/logs/02_bwa_index.log" "BWA index" \
    bwa index -p "$BWA_PREFIX" "$REFERENCE_LINK"
fi

if [[ "$RUN_GATK4" == 1 && ! -s "$REFERENCE_DICT" ]]; then
  run_logged "$PROCESSED_DIR/logs/02_gatk_dictionary.log" "GATK sequence dictionary" \
    gatk --java-options "$GATK_JAVA_OPTIONS" CreateSequenceDictionary \
      -R "$REFERENCE_LINK" -O "$REFERENCE_DICT"
fi

[[ -s "$REFERENCE_LINK.fai" ]] || fail "samtools did not create a FASTA index" \
  "check $PROCESSED_DIR/logs/02_faidx.log and verify the FASTA is valid"
[[ -s "$BWA_PREFIX.bwt" ]] || fail "BWA index is incomplete" \
  "check $PROCESSED_DIR/logs/02_bwa_index.log and verify write permission"
if [[ "$RUN_GATK4" == 1 ]]; then
  [[ -s "$REFERENCE_DICT" ]] || fail "GATK sequence dictionary is incomplete" \
    "check $PROCESSED_DIR/logs/02_gatk_dictionary.log and verify write permission"
fi
echo "[2/9] OK: reference link, FASTA index, and BWA index are ready"
