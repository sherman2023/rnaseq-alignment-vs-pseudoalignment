#!/usr/bin/env bash
set -euo pipefail

# One-command, non-interactive end-to-end run:
#   - downloads micromamba (if needed)
#   - creates the conda env (in-repo) from envs/conda_env.yml
#   - downloads airway FASTQs (SRA)
#   - downloads references + builds STAR/Salmon indices
#   - runs Approach A (STAR+featureCounts)
#   - runs Approach B (Salmon)
#   - runs DESeq2 + A vs B comparison

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
THREADS="${THREADS:-$( (command -v nproc >/dev/null 2>&1 && nproc) || echo 8 )}"

# Keep all micromamba state inside the repo
export MAMBA_ROOT_PREFIX="${ROOT_DIR}/.micromamba"

have() { command -v "$1" >/dev/null 2>&1; }

# Determine micromamba platform string
uname_s="$(uname -s)"
uname_m="$(uname -m)"
PLATFORM=""
case "${uname_s}" in
  Linux)
    case "${uname_m}" in
      x86_64) PLATFORM="linux-64";;
      aarch64|arm64) PLATFORM="linux-aarch64";;
      ppc64le) PLATFORM="linux-ppc64le";;
      *) echo "ERROR: Unsupported Linux arch: ${uname_m}" >&2; exit 1;;
    esac
    ;;
  Darwin)
    case "${uname_m}" in
      x86_64) PLATFORM="osx-64";;
      arm64) PLATFORM="osx-arm64";;
      *) echo "ERROR: Unsupported macOS arch: ${uname_m}" >&2; exit 1;;
    esac
    ;;
  *)
    echo "ERROR: Unsupported OS: ${uname_s}" >&2
    exit 1
    ;;

esac

# Locate or install micromamba
MICROMAMBA=""
if have micromamba; then
  MICROMAMBA="$(command -v micromamba)"
else
  mkdir -p "${ROOT_DIR}/.tools"
  if [[ ! -x "${ROOT_DIR}/.tools/bin/micromamba" ]]; then
    echo "[run_all] downloading micromamba for ${PLATFORM}"
    if have curl; then
      curl -Ls "https://micro.mamba.pm/api/micromamba/${PLATFORM}/latest" | tar -xvj -C "${ROOT_DIR}/.tools" bin/micromamba
    elif have wget; then
      wget -qO- "https://micro.mamba.pm/api/micromamba/${PLATFORM}/latest" | tar -xvj -C "${ROOT_DIR}/.tools" bin/micromamba
    else
      echo "ERROR: Need curl or wget to download micromamba." >&2
      exit 1
    fi
    chmod +x "${ROOT_DIR}/.tools/bin/micromamba"
  fi
  MICROMAMBA="${ROOT_DIR}/.tools/bin/micromamba"
fi

echo "[run_all] micromamba: ${MICROMAMBA}"
echo "[run_all] THREADS: ${THREADS}"

ENV_NAME="rnaseq-avsb"

# Create env if missing
if ! "${MICROMAMBA}" env list | awk '{print $1}' | grep -qx "${ENV_NAME}"; then
  echo "[run_all] creating environment ${ENV_NAME}"
  "${MICROMAMBA}" create -y -f "${ROOT_DIR}/envs/conda_env.yml"
fi

MMRUN=("${MICROMAMBA}" run -n "${ENV_NAME}")

# Run pipeline
"${MMRUN[@]}" bash "${ROOT_DIR}/scripts/01_download_airway_fastq.sh"
"${MMRUN[@]}" bash "${ROOT_DIR}/scripts/00_setup_references.sh"
"${MMRUN[@]}" bash "${ROOT_DIR}/scripts/02A_star_featurecounts.sh"
"${MMRUN[@]}" bash "${ROOT_DIR}/scripts/02B_salmon_quant.sh"
"${MMRUN[@]}" Rscript "${ROOT_DIR}/analysis/03_deseq2_compare.R"

echo "[run_all] Complete. Key outputs in: ${ROOT_DIR}/results/deseq2/"
