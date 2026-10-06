# 01_prepare_data.R ------------------------------------------------------
# Extract the raw count matrix and sample metadata from the `airway`
# Bioconductor package (Himes et al. 2014: human airway smooth muscle cells
# treated with dexamethasone) and write them to data/ as plain CSV files.

source("scripts/00_setup.R")
library(airway)
library(SummarizedExperiment)

data(airway)
airway

# Sample metadata -------------------------------------------------------------
meta <- as.data.frame(colData(airway))

sample_info <- data.frame(
  cellLine      = factor(meta$cell),
  dexamethasone = factor(ifelse(meta$dex == "trt", "treated", "untreated"),
                         levels = c("untreated", "treated")),
  row.names     = rownames(meta)
)

print(sample_info)

# Raw counts -----------------------------------------------------------------
counts_data <- assay(airway, "counts")

# Sanity checks
stopifnot(identical(colnames(counts_data), rownames(sample_info)))
stopifnot(all(counts_data >= 0), all(counts_data == round(counts_data)))

write.csv(sample_info, file.path(dir_data, "sample_info.csv"), quote = FALSE)
write.csv(counts_data, file.path(dir_data, "counts_data.csv"), quote = FALSE)

message("Wrote ", nrow(counts_data), " genes x ", ncol(counts_data), " samples to ", dir_data, "/")
