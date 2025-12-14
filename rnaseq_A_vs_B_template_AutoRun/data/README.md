# data/

This directory is intended for **downloaded input data** (FASTQ, SRA cache, etc.).

Recommended structure:
- `data/sra/`   : SRA downloads (from `prefetch`)
- `data/fastq/` : gzipped FASTQs (from `fasterq-dump`)

This repository's `.gitignore` is configured to avoid committing large sequencing files.
