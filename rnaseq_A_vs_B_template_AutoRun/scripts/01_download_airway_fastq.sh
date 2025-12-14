#!/usr/bin/env bash
set -euo pipefail

# Download paired-end FASTQs for the airway dataset runs listed in config/samples.tsv
# Default behavior:
#   - stores .sra files in data/sra/
#   - stores gzipped FASTQs in data/fastq/
#
# Requirements: sra-tools (prefetch, fasterq-dump), pigz

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

THREADS="${THREADS:-8}"
SAMPLES_TSV="${SAMPLES_TSV:-${ROOT_DIR}/config/samples.tsv}"
SRA_DIR="${SRA_DIR:-${ROOT_DIR}/data/sra}"
FASTQ_DIR="${FASTQ_DIR:-${ROOT_DIR}/data/fastq}"

mkdir -p "${SRA_DIR}" "${FASTQ_DIR}"

if [[ ! -f "${SAMPLES_TSV}" ]]; then
  echo "ERROR: samples table not found: ${SAMPLES_TSV}" >&2
  exit 1
fi

mapfile -t RUNS < <(awk -F'\t' 'NR>1 {print $1}' "${SAMPLES_TSV}")

echo "[download] runs: ${#RUNS[@]}"
echo "[download] SRA_DIR:   ${SRA_DIR}"
echo "[download] FASTQ_DIR: ${FASTQ_DIR}"
echo "[download] THREADS:   ${THREADS}"

for run in "${RUNS[@]}"; do
  echo "== ${run} =="

  # 1) Download .sra (kept in project-local cache)
  prefetch -O "${SRA_DIR}" "${run}"

  # 2) Convert to FASTQ (paired-end)
  # prefetch typically creates: data/sra/SRRxxxx/SRRxxxx.sra
  sra_file="${SRA_DIR}/${run}/${run}.sra"
  if [[ -f "${sra_file}" ]]; then
    fasterq-dump --split-files --threads "${THREADS}" --outdir "${FASTQ_DIR}" "${sra_file}"
  else
    # fallback: fasterq-dump will search the default SRA cache
    fasterq-dump --split-files --threads "${THREADS}" --outdir "${FASTQ_DIR}" "${run}"
  fi

  # 3) Compress (saves disk significantly)
  for mate in 1 2; do
    fq="${FASTQ_DIR}/${run}_${mate}.fastq"
    if [[ -f "${fq}" ]]; then
      pigz -p "${THREADS}" "${fq}"
    fi
  done
done

echo "Done. FASTQs: ${FASTQ_DIR}/*.fastq.gz"
