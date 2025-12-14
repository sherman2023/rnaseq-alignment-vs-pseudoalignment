#!/usr/bin/env bash
set -euo pipefail

# Approach A: STAR (genome alignment) -> featureCounts (gene-level counts)
#
# Assumes paired-end FASTQs named:
#   data/fastq/<RUN>_1.fastq.gz
#   data/fastq/<RUN>_2.fastq.gz
#
# You must provide a STAR genome index and an annotation GTF.
# Defaults:
#   ref/star_index/      (STAR index)
#   ref/annotation.gtf   (GTF)

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

THREADS="${THREADS:-8}"
SAMPLES_TSV="${SAMPLES_TSV:-${ROOT_DIR}/config/samples.tsv}"
GENOME_DIR="${GENOME_DIR:-${ROOT_DIR}/ref/star_index}"
GTF="${GTF:-${ROOT_DIR}/ref/annotation.gtf}"
FASTQ_DIR="${FASTQ_DIR:-${ROOT_DIR}/data/fastq}"

OUTDIR="${OUTDIR:-${ROOT_DIR}/results/approachA_star}"
COUNTS_DIR="${COUNTS_DIR:-${ROOT_DIR}/results/approachA_featurecounts}"

mkdir -p "${OUTDIR}" "${COUNTS_DIR}"

if [[ ! -f "${SAMPLES_TSV}" ]]; then
  echo "ERROR: samples table not found: ${SAMPLES_TSV}" >&2
  exit 1
fi
if [[ ! -d "${GENOME_DIR}" ]]; then
  echo "ERROR: STAR index directory not found: ${GENOME_DIR}" >&2
  echo "Build or download a STAR genome index and set GENOME_DIR accordingly." >&2
  exit 1
fi
if [[ ! -f "${GTF}" ]]; then
  echo "ERROR: annotation GTF not found: ${GTF}" >&2
  echo "Place an annotation GTF at ref/annotation.gtf or set GTF accordingly." >&2
  exit 1
fi

echo "[A] THREADS=${THREADS}"
echo "[A] SAMPLES_TSV=${SAMPLES_TSV}"
echo "[A] GENOME_DIR=${GENOME_DIR}"
echo "[A] GTF=${GTF}"
echo "[A] FASTQ_DIR=${FASTQ_DIR}"
echo "[A] OUTDIR=${OUTDIR}"
echo "[A] COUNTS_DIR=${COUNTS_DIR}"

# Align each sample
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

  prefix="${OUTDIR}/${run}."
  echo "[STAR] ${run}"

  STAR \
    --genomeDir "${GENOME_DIR}" \
    --readFilesIn "${r1}" "${r2}" \
    --readFilesCommand zcat \
    --runThreadN "${THREADS}" \
    --outFileNamePrefix "${prefix}" \
    --outSAMtype BAM SortedByCoordinate

  bam="${OUTDIR}/${run}.Aligned.sortedByCoord.out.bam"
  if [[ ! -f "${bam}" ]]; then
    echo "ERROR: STAR did not produce expected BAM: ${bam}" >&2
    exit 1
  fi

  samtools index -@ "${THREADS}" "${bam}"
done < "${SAMPLES_TSV}"

# Count with featureCounts
mapfile -t BAMS < <(ls -1 "${OUTDIR}"/*.Aligned.sortedByCoord.out.bam 2>/dev/null || true)
if [[ "${#BAMS[@]}" -eq 0 ]]; then
  echo "ERROR: No BAM files found in ${OUTDIR}. Did STAR run successfully?" >&2
  exit 1
fi

echo "[featureCounts] Counting on ${#BAMS[@]} BAM files"

featureCounts \
  -T "${THREADS}" \
  -p -B -C \
  -a "${GTF}" \
  -o "${COUNTS_DIR}/featureCounts.counts.tsv" \
  "${BAMS[@]}"

echo "Done. Counts: ${COUNTS_DIR}/featureCounts.counts.tsv"
