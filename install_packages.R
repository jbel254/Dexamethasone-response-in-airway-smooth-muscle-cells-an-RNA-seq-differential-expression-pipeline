# install_packages.R -----------------------------------------------------
# One-off installation of everything the pipeline needs.
# For stricter reproducibility consider renv::init() afterwards.

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

cran <- c("tidyverse", "pheatmap", "ggrepel", "msigdbr", "knitr", "rmarkdown")
bioc <- c("airway", "DESeq2", "apeglm", "AnnotationDbi", "org.Hs.eg.db",
          "limma", "edgeR", "clusterProfiler", "enrichplot", "SummarizedExperiment")

install.packages(setdiff(cran, rownames(installed.packages())))
BiocManager::install(setdiff(bioc, rownames(installed.packages())), update = FALSE)
