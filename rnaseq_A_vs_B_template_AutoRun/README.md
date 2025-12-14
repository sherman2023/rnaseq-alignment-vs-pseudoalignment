# RNA-seq differential expression: alignment-based vs pseudoalignment/quasi-mapping

This repository is a **reproducible A vs B template** for bulk RNA-seq differential expression (DE) using:

- **Approach A (alignment-based)**: **STAR** (splice-aware genome alignment) → **featureCounts** → **DESeq2**
- **Approach B (pseudoalignment/quasi-mapping)**: **Salmon** → **tximport** → **DESeq2**

It is intended to support a final project paper comparing the two approaches, including compute/runtime tradeoffs and differences in DE results.

## Quickstart (one command)

From a fresh clone, run:

```bash
bash run_all.sh
```

What this does (non-interactive):
1. Downloads **micromamba** (if you don’t already have it in PATH).
2. Creates the conda environment **in the repo** (under `.micromamba/`) using `envs/conda_env.yml`.
3. Downloads the **airway** paired-end FASTQs listed in `config/samples.tsv` (via SRA tools).
4. Downloads **GENCODE (default: release 49) + GRCh38** reference files and builds:
   - `ref/star_index/` (STAR genome index)
   - `ref/salmon_index/` (Salmon transcriptome index)
   - `ref/tx2gene.tsv` (transcript → gene mapping)
5. Runs both pipelines (A and B) and then the DE + comparison script.

Key outputs end up in:
- `results/deseq2/`

### Notes
- This workflow downloads **large** reference and sequencing files and may take a long time on a laptop.
- You can control parallelism with `THREADS`, e.g.:
  ```bash
  THREADS=16 bash run_all.sh
  ```
- You can pin a different GENCODE release (must exist on the GENCODE FTP), e.g.:
  ```bash
  GENCODE_RELEASE=48 THREADS=16 bash run_all.sh
  ```

## Repository layout

- `config/samples.tsv` : sample sheet (run / cell line / treatment)
- `scripts/`           : bash + R helper scripts
- `analysis/`          : R script that runs DESeq2 for both approaches and compares outputs
- `ref/`               : downloaded references + indices (created by script)
- `data/`              : downloaded FASTQs / SRA cache
- `results/`           : outputs and plots

Large files are ignored by git via `.gitignore`.

## Step-by-step (manual)

If you prefer to run each step yourself (after you have activated the environment):

```bash
bash scripts/01_download_airway_fastq.sh
bash scripts/00_setup_references.sh
bash scripts/02A_star_featurecounts.sh
bash scripts/02B_salmon_quant.sh
Rscript analysis/03_deseq2_compare.R
```

## Environment

The environment file is at:
- `envs/conda_env.yml`

It installs STAR, Salmon, featureCounts (Subread), SRA tools, and the required R/Bioconductor packages.
