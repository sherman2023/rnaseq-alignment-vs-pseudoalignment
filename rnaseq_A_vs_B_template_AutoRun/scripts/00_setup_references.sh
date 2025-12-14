#!/usr/bin/env bash
set -euo pipefail

# Download reference files (GENCODE + GRCh38) and build indices for STAR and Salmon.
#
# Non-interactive + idempotent: safe to re-run.
#
# Default pinned reference:
#   GENCODE human release 49 (GRCh38.p14)
#
# You may override by exporting GENCODE_RELEASE (e.g., 49).
#
# Outputs (default):
#   ref/genome.fa
#   ref/annotation.gtf
#   ref/transcripts.fa
#   ref/tx2gene.tsv
#   ref/star_index/
#   ref/salmon_index/

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

THREADS="${THREADS:-8}"
GENCODE_RELEASE="${GENCODE_RELEASE:-49}"

REF_DIR="${REF_DIR:-${ROOT_DIR}/ref}"
FASTQ_DIR="${FASTQ_DIR:-${ROOT_DIR}/data/fastq}"

mkdir -p "${REF_DIR}"

# ----------------------------
# Helpers
# ----------------------------
have() { command -v "$1" >/dev/null 2>&1; }

download() {
  local url="$1"
  local out="$2"
  if [[ -f "${out}" ]]; then
    echo "[ref] exists: ${out}"
    return 0
  fi
  echo "[ref] downloading: ${url}"
  if have curl; then
    curl -L --retry 5 --retry-delay 2 -o "${out}" "${url}"
  elif have wget; then
    wget -O "${out}" "${url}"
  else
    echo "ERROR: neither curl nor wget found (needed to download references)." >&2
    exit 1
  fi
}

decompress_gz() {
  local in_gz="$1"
  local out="$2"
  if [[ -f "${out}" ]]; then
    echo "[ref] exists: ${out}"
    return 0
  fi
  echo "[ref] decompressing: ${in_gz} -> ${out}"
  if have pigz; then
    pigz -dc "${in_gz}" > "${out}"
  else
    gzip -dc "${in_gz}" > "${out}"
  fi
}

# ----------------------------
# Determine read length (for STAR sjdbOverhang and Salmon k)
# ----------------------------
READLEN="${READLEN:-}"
if [[ -z "${READLEN}" ]]; then
  first_run="$(awk -F'\t' 'NR==2{print $1; exit}' "${ROOT_DIR}/config/samples.tsv" 2>/dev/null || true)"
  if [[ -n "${first_run}" && -f "${FASTQ_DIR}/${first_run}_1.fastq.gz" ]]; then
    READLEN="$(
      if have pigz; then
        pigz -dc "${FASTQ_DIR}/${first_run}_1.fastq.gz"
      else
        gzip -dc "${FASTQ_DIR}/${first_run}_1.fastq.gz"
      fi | awk 'NR==2 {print length($0); exit}'
    )"
  fi
fi

if [[ -z "${READLEN}" ]]; then
  READLEN=100
  echo "[ref] WARNING: Could not infer read length from FASTQ; defaulting READLEN=${READLEN}" >&2
else
  echo "[ref] inferred READLEN=${READLEN}"
fi

SJDB_OVERHANG=$((READLEN - 1))
if [[ "${SJDB_OVERHANG}" -lt 1 ]]; then
  SJDB_OVERHANG=100
fi

SALMON_K=31
if [[ "${READLEN}" -le "${SALMON_K}" ]]; then
  SALMON_K=$((READLEN - 1))
fi
if [[ "${SALMON_K}" -lt 15 ]]; then
  SALMON_K=15
fi
if (( SALMON_K % 2 == 0 )); then
  SALMON_K=$((SALMON_K - 1))
fi

echo "[ref] STAR sjdbOverhang=${SJDB_OVERHANG}"
echo "[ref] Salmon k=${SALMON_K}"

# ----------------------------
# Download reference files (GENCODE)
# ----------------------------
BASE_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${GENCODE_RELEASE}"

GENOME_GZ="${REF_DIR}/GRCh38.primary_assembly.genome.fa.gz"
GTF_GZ="${REF_DIR}/gencode.v${GENCODE_RELEASE}.primary_assembly.annotation.gtf.gz"
TX_GZ="${REF_DIR}/gencode.v${GENCODE_RELEASE}.transcripts.fa.gz"

download "${BASE_URL}/GRCh38.primary_assembly.genome.fa.gz" "${GENOME_GZ}"
download "${BASE_URL}/gencode.v${GENCODE_RELEASE}.primary_assembly.annotation.gtf.gz" "${GTF_GZ}"
download "${BASE_URL}/gencode.v${GENCODE_RELEASE}.transcripts.fa.gz" "${TX_GZ}"

# Decompress to standardized filenames used by downstream scripts

decompress_gz "${GENOME_GZ}" "${REF_DIR}/genome.fa"
decompress_gz "${GTF_GZ}" "${REF_DIR}/annotation.gtf"
decompress_gz "${TX_GZ}" "${REF_DIR}/transcripts.fa"

# Record provenance
cat > "${REF_DIR}/REFERENCE_SOURCES.txt" <<EOP
GENCODE_RELEASE=${GENCODE_RELEASE}
Genome (GRCh38 primary assembly): ${BASE_URL}/GRCh38.primary_assembly.genome.fa.gz
Annotation (GTF): ${BASE_URL}/gencode.v${GENCODE_RELEASE}.primary_assembly.annotation.gtf.gz
Transcripts (FASTA): ${BASE_URL}/gencode.v${GENCODE_RELEASE}.transcripts.fa.gz
EOP

# ----------------------------
# Build tx2gene.tsv from GENCODE transcript FASTA headers
# ----------------------------
TX2GENE="${REF_DIR}/tx2gene.tsv"
if [[ ! -f "${TX2GENE}" ]]; then
  echo "[ref] building tx2gene.tsv from transcripts.fa headers"
  # GENCODE transcript fasta headers look like:
  # >ENST...|ENSG...|...
  # We map transcript_id -> gene_id using the first two '|' delimited fields.
  awk -F'[| ]' '/^>/{tx=$1; sub(/^>/,"",tx); gene=$2; if (tx!="" && gene!="") print tx"\t"gene}' "${REF_DIR}/transcripts.fa" \
    | sort -u > "${TX2GENE}"
fi

# ----------------------------
# Build STAR index (genomeGenerate)
# ----------------------------
STAR_INDEX="${REF_DIR}/star_index"
if [[ ! -f "${STAR_INDEX}/SA" ]]; then
  echo "[ref] building STAR index at ${STAR_INDEX}"
  mkdir -p "${STAR_INDEX}"
  STAR --runThreadN "${THREADS}" \
       --runMode genomeGenerate \
       --genomeDir "${STAR_INDEX}" \
       --genomeFastaFiles "${REF_DIR}/genome.fa" \
       --sjdbGTFfile "${REF_DIR}/annotation.gtf" \
       --sjdbOverhang "${SJDB_OVERHANG}" \
       --genomeSAindexNbases 14
else
  echo "[ref] STAR index exists: ${STAR_INDEX}"
fi

# ----------------------------
# Build Salmon index
# ----------------------------
SALMON_INDEX="${REF_DIR}/salmon_index"
if [[ ! -f "${SALMON_INDEX}/hash.bin" ]]; then
  echo "[ref] building Salmon index at ${SALMON_INDEX}"
  mkdir -p "${SALMON_INDEX}"
  salmon index \
    --gencode \
    -t "${REF_DIR}/transcripts.fa" \
    -i "${SALMON_INDEX}" \
    -k "${SALMON_K}"
else
  echo "[ref] Salmon index exists: ${SALMON_INDEX}"
fi

echo "[ref] Done. References in ${REF_DIR}"
