# 04_visualise_results.R -------------------------------------------------
# Volcano plot, MA plot, top-gene heatmap and counts for known
# glucocorticoid-responsive genes (a biological sanity check).

source("scripts/00_setup.R")
library(DESeq2)
library(pheatmap)
library(ggrepel)

dds    <- readRDS(file.path(dir_res, "dds.rds"))
res_df <- readRDS(file.path(dir_res, "res_df.rds"))
vsd    <- readRDS(file.path(dir_res, "vsd.rds"))

plot_df <- res_df |>
  filter(!is.na(padj)) |>
  mutate(label = ifelse(is.na(symbol), ensembl, symbol))

# 1. Volcano plot -----------------------------------------------------------------
top_label <- plot_df |>
  filter(direction != "Not significant") |>
  slice_min(padj, n = 15, with_ties = FALSE)

p_volcano <- ggplot(plot_df, aes(lfc_shrunk, -log10(padj), colour = direction)) +
  geom_point(alpha = 0.6, size = 1.2) +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", colour = "grey40") +
  geom_hline(yintercept = -log10(padj_cut), linetype = "dashed", colour = "grey40") +
  geom_text_repel(data = top_label, aes(label = label),
                  size = 3, max.overlaps = Inf, show.legend = FALSE) +
  scale_colour_manual(values = col_dir) +
  labs(x = "log2 fold change (shrunken, treated vs untreated)",
       y = expression(-log[10]~adjusted~p),
       colour = NULL,
       title = "Dexamethasone response: volcano plot")
save_plot(p_volcano, "volcano.png", width = 7.5, height = 6)

# 2. MA plot -----------------------------------------------------------------------
p_ma <- ggplot(plot_df, aes(baseMean, lfc_shrunk, colour = direction)) +
  geom_point(alpha = 0.6, size = 1) +
  geom_hline(yintercept = 0, colour = "grey30") +
  scale_x_log10() +
  scale_colour_manual(values = col_dir) +
  labs(x = "Mean of normalised counts", y = "Shrunken log2 fold change",
       colour = NULL, title = "MA plot (apeglm-shrunken LFC)")
save_plot(p_ma, "ma_plot.png")

# 3. Heatmap of top 40 significant genes ----------------------------------------------
top_genes <- plot_df |>
  filter(direction != "Not significant") |>
  slice_min(padj, n = 40, with_ties = FALSE)

mat <- assay(vsd)[top_genes$ensembl, ]
rownames(mat) <- make.unique(top_genes$label)

annot <- as.data.frame(colData(vsd))[, c("dexamethasone", "cellLine")]

pheatmap(mat,
         scale = "row",
         annotation_col = annot,
         annotation_colors = list(dexamethasone = col_trt),
         show_colnames = FALSE,
         main = "Top 40 DE genes (row z-score of VST counts)",
         filename = file.path(dir_fig, "heatmap_top40.png"),
         width = 7, height = 9)

# 4. Known glucocorticoid-responsive genes ---------------------------------------------
# Expected to be up-regulated by dexamethasone.
known <- c("DUSP1", "KLF15", "PER1", "TSC22D3", "FKBP5", "CDKN1C")

known_ids <- plot_df |>
  filter(symbol %in% known) |>
  distinct(symbol, .keep_all = TRUE)

known_counts <- map_dfr(seq_len(nrow(known_ids)), function(i) {
  plotCounts(dds, gene = known_ids$ensembl[i],
             intgroup = c("dexamethasone", "cellLine"), returnData = TRUE) |>
    mutate(symbol = known_ids$symbol[i])
})

p_known <- ggplot(known_counts,
                  aes(dexamethasone, count, group = cellLine, colour = cellLine)) +
  geom_line(alpha = 0.6) +
  geom_point(size = 2.5) +
  scale_y_log10() +
  facet_wrap(~ symbol, scales = "free_y") +
  labs(x = NULL, y = "Normalised count (log10)", colour = "Cell line",
       title = "Known glucocorticoid-responsive genes",
       subtitle = "Lines connect paired samples from the same cell line")
save_plot(p_known, "known_gc_genes.png", width = 9, height = 6)

# 5. Top-gene tables ----------------------------------------------------------------------
res_df |>
  filter(direction == "Up") |> slice_min(padj, n = 20, with_ties = FALSE) |>
  select(symbol, ensembl, baseMean, lfc_shrunk, padj) |>
  write.csv(file.path(dir_res, "top20_up.csv"), row.names = FALSE)

res_df |>
  filter(direction == "Down") |> slice_min(padj, n = 20, with_ties = FALSE) |>
  select(symbol, ensembl, baseMean, lfc_shrunk, padj) |>
  write.csv(file.path(dir_res, "top20_down.csv"), row.names = FALSE)

message("Result figures written to ", dir_fig, "/")
