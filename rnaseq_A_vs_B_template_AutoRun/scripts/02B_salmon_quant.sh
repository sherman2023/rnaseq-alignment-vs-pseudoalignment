#!/usr/bin/env bash
set -euo pipefail

# Approach B: Salmon quantification (transcriptome-level) -> tximport (in R) -> DESeq2
#
# Assumes paired-end FASTQs named:
#   data/fastq/<RUN>_1.fastq.gz
#   data/fastq/<RUN>_2.fastq.gz
#
# You must provide a Salmon index (transcriptome index).
# Default: ref/salmon_index/

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

THREADS="${THREADS:-8}"
SAMPLES_TSV="${SAMPLES_TSV:-${ROOT_DIR}/config/samples.tsv}"
SALMON_INDEX="${SALMON_INDEX:-${ROOT_DIR}/ref/salmon_index}"
GENE_MAP="${GENE_MAP:-${ROOT_DIR}/ref/tx2gene.tsv}"
FASTQ_DIR="${FASTQ_DIR:-${ROOT_DIR}/data/fastq}"

OUTDIR="${OUTDIR:-${ROOT_DIR}/results/approachB_salmon}"
mkdir -p "${OUTDIR}"

if [[ ! -f "${SAMPLES_TSV}" ]]; then
  echo "ERROR: samples table not found: ${SAMPLES_TSV}" >&2
  exit 1
fi
if [[ ! -d "${SALMON_INDEX}" ]]; then
  echo "ERROR: Salmon index directory not found: ${SALMON_INDEX}" >&2
  echo "Build or download a Salmon transcriptome index and set SALMON_INDEX accordingly." >&2
  exit 1
fi

echo "[B] THREADS=${THREADS}"
echo "[B] SAMPLES_TSV=${SAMPLES_TSV}"
echo "[B] SALMON_INDEX=${SALMON_INDEX}"
echo "[B] GENE_MAP=${GENE_MAP}"
echo "[B] FASTQ_DIR=${FASTQ_DIR}"
echo "[B] OUTDIR=${OUTDIR}"

GENE_MAP_ARGS=()
if [[ -f "${GENE_MAP}" ]]; then
  # If present, Salmon will also output quant.genes.sf (gene-level) for convenience.
  GENE_MAP_ARGS=(--geneMap "${GENE_MAP}")
fi

while IFS=$'\t' read -r run cell dex; do
  [[ "${run}" == "run" ]] && continue

  r1="${FASTQ_DIR}/${run}_1.fastq.gz"
  r2="${FASTQ_DIR}/${run}_2.fastq.gz"

  if [[ ! -f "${r1}" || ! -f "${r2}" ]]; then
    echo "ERROR: Missing FASTQs for ${run}:" >&2
    echo "  ${r1}" >&2
    echo "  ${r2}" >&2
    exit 1
  fi

  echo "[Salmon] ${run}"

  salmon quant \
    -i "${SALMON_INDEX}" \
    -l A \
    -1 "${r1}" \
    -2 "${r2}" \
    -p "${THREADS}" \
    --gcBias \
    "${GENE_MAP_ARGS[@]}" \
    -o "${OUTDIR}/${run}"
done < "${SAMPLES_TSV}"

echo "Done. Quant dirs: ${OUTDIR}/<RUN>/quant.sf"
