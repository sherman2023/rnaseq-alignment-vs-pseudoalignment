#!/usr/bin/env Rscript

# Build tx2gene.tsv from a GTF (Ensembl/Gencode style) using rtracklayer.
# Usage:
#   Rscript scripts/make_tx2gene_from_gtf.R ref/annotation.gtf ref/tx2gene.tsv

suppressPackageStartupMessages({
  library(rtracklayer)
  library(dplyr)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: make_tx2gene_from_gtf.R <annotation.gtf> <tx2gene.tsv>")
}
gtf_file <- args[[1]]
out_file <- args[[2]]

gtf <- import(gtf_file)

# Keep transcript records if available; otherwise use exons and unique pairs
if ("type" %in% names(gtf)) {
  gtf_tx <- gtf[gtf$type %in% c("transcript","exon")]
} else {
  gtf_tx <- gtf
}

m <- mcols(gtf_tx)
if (!all(c("transcript_id","gene_id") %in% colnames(m))) {
  stop("GTF does not contain transcript_id and gene_id attributes. Check your annotation.")
}

tx2gene <- tibble(
  tx = as.character(m$transcript_id),
  gene = as.character(m$gene_id)
) %>%
  distinct() %>%
  filter(!is.na(tx), !is.na(gene))

write_tsv(tx2gene, out_file, col_names = FALSE)
message("Wrote: ", out_file, " (", nrow(tx2gene), " transcript→gene mappings)")
