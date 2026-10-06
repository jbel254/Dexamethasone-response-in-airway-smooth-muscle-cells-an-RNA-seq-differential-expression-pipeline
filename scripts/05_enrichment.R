# 05_enrichment.R --------------------------------------------------------
# Functional interpretation:
#   (a) GO Biological Process over-representation (up / down genes separately)
#   (b) KEGG over-representation (needs internet; skipped if unavailable)
#   (c) GSEA on all genes ranked by the Wald statistic, using MSigDB Hallmark

source("scripts/00_setup.R")
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(msigdbr)

res_df <- readRDS(file.path(dir_res, "res_df.rds"))

# Background = every gene that was actually tested
universe <- res_df |>
  filter(!is.na(padj), !is.na(entrez)) |>
  pull(entrez) |> unique()

gene_sets <- list(
  up   = res_df |> filter(direction == "Up",   !is.na(entrez)) |> pull(entrez) |> unique(),
  down = res_df |> filter(direction == "Down", !is.na(entrez)) |> pull(entrez) |> unique()
)
message("Up: ", length(gene_sets$up), " | Down: ", length(gene_sets$down),
        " | Background: ", length(universe))

has_rows <- function(x) !is.null(x) && nrow(as.data.frame(x)) > 0

# (a) GO Biological Process ---------------------------------------------------------------
for (nm in names(gene_sets)) {
  ego <- enrichGO(gene          = gene_sets[[nm]],
                  universe      = universe,
                  OrgDb         = org.Hs.eg.db,
                  keyType       = "ENTREZID",
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
  if (has_rows(ego)) {
    ego <- simplify(ego, cutoff = 0.7)   # collapse redundant GO terms
    write.csv(as.data.frame(ego), file.path(dir_res, paste0("go_bp_", nm, ".csv")),
              row.names = FALSE)
    save_plot(dotplot(ego, showCategory = 15) + ggtitle(paste("GO BP -", nm, "regulated genes")),
              paste0("go_bp_", nm, ".png"), width = 8, height = 7)
  } else {
    message("No significant GO BP terms for: ", nm)
  }
}

# (b) KEGG (optional, queries the KEGG web service) -----------------------------------------
kegg <- tryCatch(
  enrichKEGG(gene = gene_sets$up, universe = universe, organism = "hsa", pvalueCutoff = 0.05),
  error = function(e) { message("KEGG skipped: ", conditionMessage(e)); NULL }
)
if (has_rows(kegg)) {
  write.csv(as.data.frame(kegg), file.path(dir_res, "kegg_up.csv"), row.names = FALSE)
  save_plot(dotplot(kegg, showCategory = 15) + ggtitle("KEGG - up-regulated genes"),
            "kegg_up.png", width = 8, height = 6)
}

# (c) GSEA: Hallmark gene sets --------------------------------------------------------------
# Rank by Wald statistic; one entry per gene symbol (largest |stat| wins)
ranked <- res_df |>
  filter(!is.na(stat), !is.na(symbol)) |>
  group_by(symbol) |>
  slice_max(abs(stat), n = 1, with_ties = FALSE) |>
  ungroup()

gene_list <- sort(setNames(ranked$stat, ranked$symbol), decreasing = TRUE)

# msigdbr changed its argument name (category -> collection) in newer versions
hallmark <- tryCatch(
  msigdbr(species = "Homo sapiens", collection = "H"),
  error = function(e) msigdbr(species = "Homo sapiens", category = "H")
)
term2gene <- hallmark |> distinct(gs_name, gene_symbol)

set.seed(42)
gsea <- GSEA(geneList      = gene_list,
             TERM2GENE     = term2gene,
             pvalueCutoff  = 0.05,
             pAdjustMethod = "BH",
             eps           = 0,
             seed          = TRUE,
             verbose       = FALSE)

if (has_rows(gsea)) {
  gsea_df <- as.data.frame(gsea)
  write.csv(gsea_df, file.path(dir_res, "gsea_hallmark.csv"), row.names = FALSE)

  p_nes <- gsea_df |>
    mutate(ID = str_remove(ID, "^HALLMARK_"),
           ID = fct_reorder(ID, NES),
           sign = ifelse(NES > 0, "Up in treated", "Down in treated")) |>
    ggplot(aes(NES, ID, fill = sign)) +
    geom_col() +
    scale_fill_manual(values = c(`Up in treated` = "#C44E52", `Down in treated` = "#4C72B0")) +
    labs(x = "Normalised enrichment score", y = NULL, fill = NULL,
         title = "GSEA: MSigDB Hallmark pathways (FDR < 0.05)")
  save_plot(p_nes, "gsea_hallmark_nes.png", width = 8, height = 6)

  # Classic running-score plots for the top 3 gene sets by |NES|
  top_ids <- gsea_df |> slice_max(abs(NES), n = 3, with_ties = FALSE) |> pull(ID)
  p_run <- gseaplot2(gsea, geneSetID = top_ids, pvalue_table = FALSE)
  save_plot(p_run, "gsea_running_scores.png", width = 8, height = 6)
} else {
  message("No significant Hallmark gene sets.")
}

message("Enrichment results written to ", dir_res, "/ and ", dir_fig, "/")
