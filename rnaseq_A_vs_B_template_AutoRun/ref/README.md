# ref/

This directory holds **reference files** and **indices** needed by STAR and Salmon.

Expected paths (defaults used by the scripts):
- `ref/annotation.gtf`      : gene annotation (GTF)
- `ref/star_index/`         : STAR genome index directory
- `ref/salmon_index/`       : Salmon transcriptome index directory
- `ref/tx2gene.tsv`         : two-column mapping: transcript_id<TAB>gene_id

The scripts do not download a reference automatically, because human references are large.
See the main `README.md` for example index-build commands.
