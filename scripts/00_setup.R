# 00_setup.R -------------------------------------------------------------
# Shared settings sourced by every script. Run all scripts from the project
# root (the folder containing README.md), e.g.  Rscript scripts/02_deseq2.R

suppressPackageStartupMessages({
  library(tidyverse)
})

set.seed(42)

# Name-clash protection -----------------------------------------------------------
# Bioconductor packages (AnnotationDbi, clusterProfiler, IRanges, S4Vectors, ...)
# export functions with the same names as dplyr verbs, e.g. AnnotationDbi::select().
# Objects defined in the global environment always win over attached packages,
# so pinning the dplyr versions here protects every script, whatever is loaded later.
select <- dplyr::select
filter <- dplyr::filter
slice  <- dplyr::slice
rename <- dplyr::rename
count  <- dplyr::count
desc   <- dplyr::desc

# Directories ---------------------------------------------------------------
dir_data <- "data"
dir_res  <- "results"
dir_fig  <- "figures"
for (d in c(dir_data, dir_res, dir_fig)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

# Analysis thresholds (change here, re-run everything) ----------------------
padj_cut <- 0.05   # adjusted p-value cutoff
lfc_cut  <- 1      # absolute log2 fold-change cutoff

# Plot defaults ------------------------------------------------------------
theme_set(theme_bw(base_size = 12))
col_trt <- c(untreated = "#4C72B0", treated = "#DD8452")
col_dir <- c(Up = "#C44E52", Down = "#4C72B0", `Not significant` = "grey70")

save_plot <- function(p, name, width = 7, height = 5) {
  ggsave(file.path(dir_fig, name), plot = p,
         width = width, height = height, dpi = 300)
}
