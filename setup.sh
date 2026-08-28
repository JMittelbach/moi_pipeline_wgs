#!/usr/bin/env bash
set -Eeuo pipefail

PIPELINE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_NAME="moi_pipeline"

usage() {
  cat <<'EOF'
Usage: ./setup.sh [--help]

Create or update the Conda environment 'moi_pipeline', install the command-line
tools (including optional FreeBayes and GATK4 variant calling) and R MOI dependencies,
and check that the executables are available.
This script does not download FASTQs, reference genomes, indexes, targets, or
results. After setup, activate the environment and edit config/pipeline.env.

Options:
  -h, --help    show this help and exit without changing anything
EOF
}

case "${1:-}" in
  "" ) ;;
  -h|--help ) usage; exit 0 ;;
  * )
    echo "[ERROR] unknown setup option: $1" >&2
    usage >&2
    exit 2
    ;;
esac

echo "MOI pipeline setup"
echo "------------------"

if ! command -v conda >/dev/null 2>&1; then
  echo "[ERROR] Conda was not found in PATH." >&2
  echo "        Install Miniconda (recommended):" >&2
  echo "        https://docs.conda.io/projects/miniconda/en/latest/" >&2
  echo "        macOS/Linux quick install instructions:" >&2
  echo "        https://www.anaconda.com/docs/getting-started/miniconda/install" >&2
  echo "        After installation, open a new terminal and run ./setup.sh again." >&2
  exit 127
fi

if ! eval "$(conda shell.bash hook 2>/dev/null)"; then
  echo "[ERROR] Conda is on PATH but its shell integration could not be loaded." >&2
  echo "        Fix: run 'conda init bash' (or 'conda init zsh'), restart the terminal," >&2
  echo "        then run ./setup.sh again." >&2
  exit 1
fi

packages=(
  "python=3.11"
  "r-base=4.3"
  "openjdk=17"
  "pip"
  "setuptools"
  "wheel"
  "r-flexmix"
  "r-remotes"
  "r-biocmanager"
  "fastp"
  "bwa"
  "samtools"
  "bcftools"
  "htslib"
  "freebayes"
  "gatk4"
)

echo "[1/3] creating/updating Conda environment: $ENV_NAME"
if conda env list | awk -v name="$ENV_NAME" '$1 == name { found=1 } END { exit(found ? 0 : 1) }'; then
  if ! conda install -y -n "$ENV_NAME" -c conda-forge -c bioconda "${packages[@]}"; then
    echo "[ERROR] Conda could not update $ENV_NAME." >&2
    echo "        Fix: try 'conda clean --all', then rerun ./setup.sh." >&2
    echo "        If this says NoWritableEnvsDirError, choose a user-writable envs directory:" >&2
    echo "        conda config --add envs_dirs \"$HOME/.conda/envs\"" >&2
    echo "        If the solver remains stuck, install micromamba:" >&2
    echo "        https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html" >&2
    exit 1
  fi
else
  if ! conda create -y -n "$ENV_NAME" -c conda-forge -c bioconda "${packages[@]}"; then
    echo "[ERROR] Conda could not create $ENV_NAME." >&2
    echo "        Fix: update Conda ('conda update -n base -c defaults conda') and rerun ./setup.sh." >&2
    echo "        If this says NoWritableEnvsDirError, choose a user-writable envs directory:" >&2
    echo "        conda config --add envs_dirs \"$HOME/.conda/envs\"" >&2
    echo "        On Apple Silicon, use a native arm64 Miniconda installation." >&2
    exit 1
  fi
fi

# Some compiler activation hooks reference optional variables (for example
# GFORTRAN) without guarding them. Temporarily relax nounset for activation,
# then restore the script's strict mode immediately afterwards.
set +u
if conda activate "$ENV_NAME"; then
  activate_status=0
else
  activate_status=$?
fi
set -u
if [[ "$activate_status" -ne 0 ]]; then
  echo "[ERROR] Could not activate the Conda environment $ENV_NAME." >&2
  echo "        Fix: run 'conda init zsh', restart the terminal, then rerun ./setup.sh." >&2
  exit 1
fi
ENV_BIN="$CONDA_PREFIX/bin"
# Conda's macOS OpenJDK layout keeps Java below lib/jvm rather than linking it
# into bin. Export it explicitly so GATK works even when the host has no Java.
if [[ -x "$CONDA_PREFIX/lib/jvm/bin/java" ]]; then
  JAVA_HOME="$CONDA_PREFIX/lib/jvm"
  export JAVA_HOME
  export PATH="$JAVA_HOME/bin:$ENV_BIN:$PATH"
else
  export PATH="$ENV_BIN:$PATH"
fi
# Keep this setup invocation itself on the same toolchain that the pipeline
# will use. (The parent shell still needs `conda activate moi_pipeline`.)

echo "[2/3] installing MOI estimator (moimix; flexmix fallback is retained)"
if ! "$ENV_BIN/Rscript" "$PIPELINE_ROOT/scripts/setup_moimix.R"; then
  if "$ENV_BIN/Rscript" -e 'quit(status=if (requireNamespace("flexmix", quietly=TRUE)) 0 else 1)' >/dev/null 2>&1; then
    echo "[WARN] moimix could not be installed; flexmix fallback is available." >&2
    echo "       The pipeline will label its MOI result as a flexmix equivalent." >&2
  else
    echo "[ERROR] Neither moimix nor flexmix is available after setup." >&2
    echo "        Fix: inspect the R error above, then rerun ./setup.sh." >&2
    exit 1
  fi
fi

echo "[3/3] checking core executables"
for tool in python python3 R Rscript fastp bwa samtools bcftools bgzip tabix freebayes gatk; do
  if [[ ! -x "$ENV_BIN/$tool" ]]; then
    echo "[ERROR] $tool is missing from $ENV_NAME." >&2
    echo "        Fix: conda activate $ENV_NAME && conda install -c conda-forge -c bioconda $tool" >&2
    exit 1
  fi
done

JAVA_BIN="$ENV_BIN/java"
if [[ ! -x "$JAVA_BIN" && -x "$CONDA_PREFIX/lib/jvm/bin/java" ]]; then
  JAVA_BIN="$CONDA_PREFIX/lib/jvm/bin/java"
fi
if [[ ! -x "$JAVA_BIN" ]]; then
  echo "[ERROR] Java 17 is missing from $ENV_NAME." >&2
  echo "        Fix: conda activate $ENV_NAME && conda install -c conda-forge openjdk=17" >&2
  exit 1
fi

python_version="$($ENV_BIN/python --version 2>&1)"
r_version="$($ENV_BIN/Rscript -e 'cat(R.version.string)' 2>/dev/null)"
echo "Using Conda Python: $python_version"
echo "Using Conda R: $r_version"

chmod +x "$PIPELINE_ROOT"/run_pipeline.sh "$PIPELINE_ROOT"/scripts/*.sh "$PIPELINE_ROOT"/scripts/*.py
echo
echo "Setup complete. Activate the environment with:"
echo "  conda activate $ENV_NAME"
echo "Then edit config/pipeline.env and run:"
echo "  ./run_pipeline.sh build-metadata"
