# run_all.R --------------------------------------------------------------
# Run the full pipeline from the project root:  Rscript run_all.R

# Guard: all paths are relative to the project root
if (!file.exists("scripts/00_setup.R")) {
  stop("Working directory is '", getwd(), "'.\n",
       "Please setwd() to the project root (the folder containing run_all.R and scripts/) and try again.",
       call. = FALSE)
}

steps <- c(
  "scripts/01_prepare_data.R",
  "scripts/02_deseq2.R",
  "scripts/03_qc.R",
  "scripts/04_visualise_results.R",
  "scripts/05_enrichment.R",
  "scripts/06_method_comparison.R"
)

for (s in steps) {
  message("\n==== Running ", s, " ====")
  source(s, echo = FALSE)
}

message("\nDone. See results/ and figures/. Render the report with: quarto render report.qmd")
