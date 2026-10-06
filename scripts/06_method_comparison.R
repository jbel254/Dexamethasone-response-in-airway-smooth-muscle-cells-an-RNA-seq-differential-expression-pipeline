# 06_method_comparison.R -------------------------------------------------
# Cross-check DESeq2 against edgeR (quasi-likelihood) and limma-voom, all
# using the same paired design (~ cellLine + dexamethasone) and the same
# pre-filtered gene set. Agreement between methods increases confidence in
# the DE list.

source("scripts/00_setup.R")
library(DESeq2)
library(edgeR)
library(limma)

dds    <- readRDS(file.path(dir_res, "dds.rds"))
res_df <- readRDS(file.path(dir_res, "res_df.rds"))

counts_mat <- counts(dds)
meta       <- as.data.frame(colData(dds))
design     <- model.matrix(~ cellLine + dexamethasone, data = meta)
coef_name  <- "dexamethasonetreated"

# edgeR (QL F-test) ---------------------------------------------------------------------
y <- DGEList(counts = counts_mat)
y <- calcNormFactors(y)
y <- estimateDisp(y, design)
fit_er <- glmQLFit(y, design)
tt_er  <- topTags(glmQLFTest(fit_er, coef = coef_name), n = Inf)$table |>
  rownames_to_column("ensembl")

# limma-voom ----------------------------------------------------------------------------
v      <- voom(y, design)
fit_lv <- eBayes(lmFit(v, design))
tt_lv  <- topTable(fit_lv, coef = coef_name, number = Inf, sort.by = "none") |>
  rownames_to_column("ensembl")

# Significant sets (padj < cutoff and |LFC| >= cutoff, unshrunken LFC for all) -------------
sig_deseq <- res_df |>
  filter(!is.na(padj), padj < padj_cut, abs(log2FoldChange) >= lfc_cut) |> pull(ensembl)
sig_edger <- tt_er |>
  filter(FDR < padj_cut, abs(logFC) >= lfc_cut) |> pull(ensembl)
sig_limma <- tt_lv |>
  filter(adj.P.Val < padj_cut, abs(logFC) >= lfc_cut) |> pull(ensembl)

# Overlap table ---------------------------------------------------------------------------
all_ids <- union(union(sig_deseq, sig_edger), sig_limma)
membership <- tibble(
  ensembl = all_ids,
  DESeq2  = all_ids %in% sig_deseq,
  edgeR   = all_ids %in% sig_edger,
  limma   = all_ids %in% sig_limma
) |>
  mutate(pattern = paste0(ifelse(DESeq2, "DESeq2 ", ""),
                          ifelse(edgeR,  "edgeR ",  ""),
                          ifelse(limma,  "limma",   "")) |> str_trim() |> str_replace_all(" ", " + "))

overlap_counts <- membership |> count(pattern, sort = TRUE)
print(overlap_counts)
write.csv(overlap_counts, file.path(dir_res, "method_overlap.csv"), row.names = FALSE)

p_overlap <- ggplot(overlap_counts, aes(n, fct_reorder(pattern, n))) +
  geom_col(fill = "grey35") +
  geom_text(aes(label = n), hjust = -0.2, size = 3.5) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = "Number of genes", y = NULL,
       title = "Overlap of significant genes across methods",
       subtitle = paste0("padj/FDR < ", padj_cut, " and |log2FC| >= ", lfc_cut))
save_plot(p_overlap, "method_overlap.png", width = 8, height = 4)

# Fold-change agreement -------------------------------------------------------------------
lfc_cmp <- res_df |>
  select(ensembl, DESeq2 = log2FoldChange) |>
  inner_join(select(tt_er, ensembl, edgeR = logFC), by = "ensembl") |>
  inner_join(select(tt_lv, ensembl, limma = logFC), by = "ensembl")

cor_tbl <- tibble(
  comparison = c("DESeq2 vs edgeR", "DESeq2 vs limma-voom"),
  pearson_r  = c(cor(lfc_cmp$DESeq2, lfc_cmp$edgeR, use = "complete.obs"),
                 cor(lfc_cmp$DESeq2, lfc_cmp$limma, use = "complete.obs"))
)
print(cor_tbl)
write.csv(cor_tbl, file.path(dir_res, "method_lfc_correlation.csv"), row.names = FALSE)

p_lfc <- ggplot(lfc_cmp, aes(DESeq2, edgeR)) +
  geom_point(alpha = 0.3, size = 0.8) +
  geom_abline(colour = "firebrick") +
  labs(x = "DESeq2 log2FC", y = "edgeR log2FC", title = "Fold-change agreement")
save_plot(p_lfc, "method_lfc_scatter.png")

message("Method comparison complete.")
