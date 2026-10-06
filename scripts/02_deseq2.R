# 02_deseq2.R ------------------------------------------------------------
# Differential expression: dexamethasone-treated vs untreated.
#
# Design: ~ cellLine + dexamethasone
#   The experiment is PAIRED (each of 4 cell lines was treated and untreated),
#   so cell line is included as a blocking factor. A naive model without it
#   is also fitted to show why this matters.

source("scripts/00_setup.R")
library(DESeq2)
library(apeglm)
library(AnnotationDbi)
library(org.Hs.eg.db)

# 1. Read data -----------------------------------------------------------------
counts_data <- read.csv(file.path(dir_data, "counts_data.csv"),
                        row.names = 1, check.names = FALSE)
sample_info <- read.csv(file.path(dir_data, "sample_info.csv"),
                        row.names = 1, stringsAsFactors = TRUE)

sample_info$dexamethasone <- relevel(sample_info$dexamethasone, ref = "untreated")

# Row names of colData must match column names of the counts, in the same order
stopifnot(all(colnames(counts_data) %in% rownames(sample_info)))
stopifnot(all(colnames(counts_data) == rownames(sample_info)))

# 2. Build DESeqDataSets -------------------------------------------------------
dds_naive <- DESeqDataSetFromMatrix(counts_data, sample_info, design = ~ dexamethasone)
dds       <- DESeqDataSetFromMatrix(counts_data, sample_info, design = ~ cellLine + dexamethasone)

# Pre-filter: keep genes with >= 10 reads in at least 4 samples
# (4 = size of the smallest group)
keep      <- rowSums(counts(dds) >= 10) >= 4
dds       <- dds[keep, ]
dds_naive <- dds_naive[keep, ]
message(sum(keep), " of ", length(keep), " genes retained after filtering")

# 3. Fit models ----------------------------------------------------------------
dds       <- DESeq(dds)
dds_naive <- DESeq(dds_naive)

resultsNames(dds)

# 4. Results -------------------------------------------------------------------
res       <- results(dds,       contrast = c("dexamethasone", "treated", "untreated"), alpha = padj_cut)
res_naive <- results(dds_naive, contrast = c("dexamethasone", "treated", "untreated"), alpha = padj_cut)

summary(res)

# Shrunken log2 fold changes (apeglm) for ranking / plotting.
# Keeping `res` means the p-values stay those of the Wald test.
res_shrunk <- lfcShrink(dds, coef = "dexamethasone_treated_vs_untreated",
                        res = res, type = "apeglm")

# 5. Annotate genes (Ensembl -> symbol / Entrez) --------------------------------
res_df <- as.data.frame(res) |>
  rownames_to_column("ensembl")

res_df$lfc_shrunk <- res_shrunk$log2FoldChange
res_df$symbol <- mapIds(org.Hs.eg.db, keys = res_df$ensembl, keytype = "ENSEMBL",
                        column = "SYMBOL", multiVals = "first")
res_df$entrez <- mapIds(org.Hs.eg.db, keys = res_df$ensembl, keytype = "ENSEMBL",
                        column = "ENTREZID", multiVals = "first")

res_df <- res_df |>
  mutate(
    direction = case_when(
      !is.na(padj) & padj < padj_cut & lfc_shrunk >=  lfc_cut ~ "Up",
      !is.na(padj) & padj < padj_cut & lfc_shrunk <= -lfc_cut ~ "Down",
      TRUE ~ "Not significant"
    )
  ) |>
  arrange(padj)

# 6. Design comparison: does blocking on cell line help? -------------------------
n_sig <- function(r) sum(r$padj < padj_cut, na.rm = TRUE)
design_comparison <- tibble(
  design = c("~ dexamethasone", "~ cellLine + dexamethasone"),
  genes_padj_lt_cutoff = c(n_sig(res_naive), n_sig(res))
)
print(design_comparison)

# 7. Save ----------------------------------------------------------------------
write.csv(res_df, file.path(dir_res, "deseq2_results_all.csv"), row.names = FALSE)
write.csv(filter(res_df, direction != "Not significant"),
          file.path(dir_res, "deseq2_results_significant.csv"), row.names = FALSE)
write.csv(design_comparison, file.path(dir_res, "design_comparison.csv"), row.names = FALSE)
writeLines(capture.output(summary(res)), file.path(dir_res, "deseq2_summary.txt"))

saveRDS(dds,    file.path(dir_res, "dds.rds"))
saveRDS(res_df, file.path(dir_res, "res_df.rds"))

writeLines(capture.output(sessionInfo()), file.path(dir_res, "sessionInfo.txt"))

message("Significant genes (padj < ", padj_cut, ", |shrunken LFC| >= ", lfc_cut, "): ",
        sum(res_df$direction != "Not significant"),
        "  (Up: ", sum(res_df$direction == "Up"),
        ", Down: ", sum(res_df$direction == "Down"), ")")
