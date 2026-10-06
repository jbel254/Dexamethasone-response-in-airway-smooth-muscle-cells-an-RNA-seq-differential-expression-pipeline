# 03_qc.R ----------------------------------------------------------------
# Quality control: library sizes, PCA (before / after removing the cell-line
# effect), sample distances, dispersion estimates and p-value histogram.

source("scripts/00_setup.R")
library(DESeq2)
library(pheatmap)
library(limma)

dds    <- readRDS(file.path(dir_res, "dds.rds"))
res_df <- readRDS(file.path(dir_res, "res_df.rds"))

# Variance-stabilising transformation (blind = TRUE for unbiased QC)
vsd <- vst(dds, blind = TRUE)
saveRDS(vsd, file.path(dir_res, "vsd.rds"))

# 1. Library sizes ---------------------------------------------------------------
lib <- tibble(
  sample        = colnames(dds),
  million_reads = colSums(counts(dds)) / 1e6,
  dexamethasone = dds$dexamethasone,
  cellLine      = dds$cellLine
)

p_lib <- ggplot(lib, aes(sample, million_reads, fill = dexamethasone)) +
  geom_col() +
  scale_fill_manual(values = col_trt) +
  labs(x = NULL, y = "Million assigned reads", title = "Library sizes") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p_lib, "qc_library_sizes.png")

# 2. PCA -------------------------------------------------------------------------
pca_plot <- function(mat, meta, title) {
  pc  <- prcomp(t(mat))
  pve <- round(100 * pc$sdev^2 / sum(pc$sdev^2), 1)
  df  <- cbind(as.data.frame(pc$x[, 1:2]), meta)
  ggplot(df, aes(PC1, PC2, colour = dexamethasone, shape = cellLine)) +
    geom_point(size = 4) +
    scale_colour_manual(values = col_trt) +
    labs(x = paste0("PC1: ", pve[1], "% variance"),
         y = paste0("PC2: ", pve[2], "% variance"),
         title = title)
}

meta <- as.data.frame(colData(dds))[, c("cellLine", "dexamethasone")]

p_pca <- pca_plot(assay(vsd), meta, "PCA of VST counts")
save_plot(p_pca, "qc_pca.png")

# Remove cell-line effect (for visualisation only, not used for testing)
vsd_corrected <- removeBatchEffect(assay(vsd), batch = vsd$cellLine,
                                   design = model.matrix(~ dexamethasone, data = meta))
p_pca_corr <- pca_plot(vsd_corrected, meta, "PCA after removing cell-line effect")
save_plot(p_pca_corr, "qc_pca_cellline_removed.png")

# 3. Sample distance heatmap --------------------------------------------------------
sample_dist <- as.matrix(dist(t(assay(vsd))))
pheatmap(sample_dist,
         annotation_col = meta,
         clustering_distance_rows = as.dist(sample_dist),
         clustering_distance_cols = as.dist(sample_dist),
         main = "Sample-to-sample distances (VST)",
         filename = file.path(dir_fig, "qc_sample_distances.png"),
         width = 7, height = 6)

# 4. Dispersion estimates -----------------------------------------------------------
png(file.path(dir_fig, "qc_dispersion.png"), width = 6, height = 5, units = "in", res = 300)
plotDispEsts(dds, main = "Dispersion estimates")
dev.off()

# 5. P-value histogram ---------------------------------------------------------------
p_hist <- res_df |>
  filter(!is.na(pvalue), baseMean > 1) |>
  ggplot(aes(pvalue)) +
  geom_histogram(breaks = seq(0, 1, 0.05), fill = "grey40", colour = "white") +
  labs(x = "Raw p-value", y = "Genes",
       title = "P-value distribution",
       subtitle = "Expect a spike near 0 and a roughly flat remainder")
save_plot(p_hist, "qc_pvalue_histogram.png")

message("QC figures written to ", dir_fig, "/")
