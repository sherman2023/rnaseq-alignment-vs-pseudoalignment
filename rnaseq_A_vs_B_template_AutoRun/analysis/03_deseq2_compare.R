#!/usr/bin/env Rscript

# Ensure paths are resolved relative to the repository root (one directory above this script).
# This allows you to run: `Rscript analysis/03_deseq2_compare.R` from any working directory.
get_repo_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 1) {
    script_path <- sub("^--file=", "", file_arg)
    script_dir <- dirname(normalizePath(script_path))
    return(normalizePath(file.path(script_dir, "..")))
  }
  return(getwd())
}
setwd(get_repo_root())

suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(tibble)
})

dir.create("results/deseq2", recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# Inputs / sample metadata
# ----------------------------
samples <- read_tsv("config/samples.tsv", show_col_types = FALSE) %>%
  mutate(cell = factor(cell),
         dex  = factor(dex, levels = c("untrt","trt")))

# ----------------------------
# Approach A: featureCounts -> DESeq2
# ----------------------------
fc_path <- "results/approachA_featurecounts/featureCounts.counts.tsv"
stopifnot(file.exists(fc_path))

fc <- read_tsv(fc_path, comment = "#", show_col_types = FALSE)
# featureCounts output: first columns are annotation; counts start at column 7
count_cols <- colnames(fc)[7:ncol(fc)]
count_mat <- as.matrix(fc[, count_cols])
rownames(count_mat) <- fc$Geneid

# Clean column names (featureCounts uses full BAM paths)
colnames(count_mat) <- gsub(".*/", "", colnames(count_mat))
colnames(count_mat) <- gsub("\\.Aligned\\.sortedByCoord\\.out\\.bam$", "", colnames(count_mat))

# Ensure sample order matches metadata
count_mat <- count_mat[, samples$run]

ddsA <- DESeqDataSetFromMatrix(countData = round(count_mat),
                              colData   = as.data.frame(samples),
                              design    = ~ cell + dex)
ddsA <- DESeq(ddsA)
resA <- results(ddsA, contrast = c("dex","trt","untrt")) %>%
  as.data.frame() %>%
  rownames_to_column("gene_id") %>%
  as_tibble()

write_tsv(resA, "results/deseq2/resA_featureCounts_DESeq2.tsv")

# ----------------------------
# Approach B: Salmon -> tximport -> DESeq2
# ----------------------------
# NOTE: To summarize Salmon to gene-level, tximport needs tx2gene (tx → gene mapping).
# This repository auto-generates ref/tx2gene.tsv from the same GENCODE transcript FASTA
# used to build the Salmon index (see scripts/00_setup_references.sh). In that default
# setup, transcript IDs match exactly and ignoreTxVersion should remain FALSE.
# (If you swap in a different annotation where transcript versions differ, you may need
# ignoreTxVersion=TRUE.)


tx2gene_path <- "ref/tx2gene.tsv"
stopifnot(file.exists(tx2gene_path))
tx2gene <- read_tsv(tx2gene_path, col_names = c("tx","gene"), show_col_types = FALSE)

files <- file.path("results/approachB_salmon", samples$run, "quant.sf")
names(files) <- samples$run
stopifnot(all(file.exists(files)))

txi <- tximport(files,
                type = "salmon",
                tx2gene = tx2gene,
                ignoreTxVersion = FALSE)

ddsB <- DESeqDataSetFromTximport(txi,
                                colData = as.data.frame(samples),
                                design  = ~ cell + dex)
ddsB <- DESeq(ddsB)
resB <- results(ddsB, contrast = c("dex","trt","untrt")) %>%
  as.data.frame() %>%
  rownames_to_column("gene_id") %>%
  as_tibble()

write_tsv(resB, "results/deseq2/resB_salmon_tximport_DESeq2.tsv")

# ----------------------------
# Compare results (A vs B)
# ----------------------------
comp <- resA %>%
  select(gene_id,
         log2FC_A = log2FoldChange,
         padj_A   = padj) %>%
  inner_join(resB %>% select(gene_id,
                             log2FC_B = log2FoldChange,
                             padj_B   = padj),
             by = "gene_id")

# Basic overlap stats
sigA <- comp %>% filter(!is.na(padj_A), padj_A < 0.05) %>% pull(gene_id) %>% unique()
sigB <- comp %>% filter(!is.na(padj_B), padj_B < 0.05) %>% pull(gene_id) %>% unique()

overlap <- length(intersect(sigA, sigB))
jaccard <- overlap / length(union(sigA, sigB))

stats <- tibble(
  n_genes_compared = nrow(comp),
  n_sig_A_padj_0.05 = length(sigA),
  n_sig_B_padj_0.05 = length(sigB),
  n_overlap = overlap,
  jaccard = jaccard,
  cor_log2FC = cor(comp$log2FC_A, comp$log2FC_B, use = "pairwise.complete.obs")
)

write_tsv(stats, "results/deseq2/AvsB_summary_stats.tsv")

# Scatter plot of log2FC
p_scatter <- ggplot(comp, aes(x = log2FC_A, y = log2FC_B)) +
  geom_point(alpha = 0.25, size = 0.8) +
  geom_abline(slope = 1, intercept = 0) +
  labs(title = "A vs B: log2 fold-changes (DESeq2)",
       x = "log2FC (STAR+featureCounts)",
       y = "log2FC (Salmon+tximport)") +
  theme_minimal(base_size = 12)

ggsave("results/deseq2/AvsB_log2FC_scatter.png", p_scatter, width = 6, height = 5, dpi = 300)

# Volcano plots (one per method) for quick visuals
volcano <- function(df, title, out_png) {
  df2 <- as.data.frame(df) %>%
    mutate(minusLog10Padj = -log10(padj)) %>%
    mutate(minusLog10Padj = ifelse(is.infinite(minusLog10Padj), NA, minusLog10Padj))
  p <- ggplot(df2, aes(x = log2FoldChange, y = minusLog10Padj)) +
    geom_point(alpha = 0.25, size = 0.8) +
    labs(title = title, x = "log2 fold-change", y = "-log10(adjusted p-value)") +
    theme_minimal(base_size = 12)
  ggsave(out_png, p, width = 6, height = 5, dpi = 300)
}

volcano(results(ddsA, contrast=c("dex","trt","untrt")),
        "Volcano: STAR+featureCounts → DESeq2",
        "results/deseq2/volcano_A.png")

volcano(results(ddsB, contrast=c("dex","trt","untrt")),
        "Volcano: Salmon+tximport → DESeq2",
        "results/deseq2/volcano_B.png")

# Save session info for reproducibility
sink("results/deseq2/sessionInfo.txt")
print(sessionInfo())
sink()

message("Done. See results/deseq2/ for tables and plots.")
