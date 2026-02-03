# 05_auc_left_vs_right.R
# Purpose: AUCell analyses of Nodal signaling in left vs right cardiomyocytes at 16 ss, 18 ss, and 20 ss.

# Setup
source(here("scripts", "00_setup.R"))
set.seed(123)

# Output folders
results_dir <- here("results", "auc_left_vs_right")
figures_dir <- here("figures", "auc_left_vs_right")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)




# 05_auc_left_vs_right.R
# Purpose:
#   - AUCell scoring of Nodal gene sets in Left vs Right cardiomyocytes
#   - Statistical comparison (Wilcoxon rank-sum) of AUCell AUC distributions
#
# Expects:
#   - 01_seurat_16ss.R has been run and produced an LR-labeled cardiomyocyte object
#     (recommended: save as results/16ss/heart16s_cluster3_LR.rds)
#
# Outputs:
#   - results/auc_left_vs_right/tables/
#   - figures/auc_left_vs_right/

# -------------------------------#
# Setup
# -------------------------------#
source(here("scripts", "00_setup.R"))
set.seed(123)

# Output folders
results_dir <- here("results", "auc_left_vs_right")
figures_dir <- here("figures", "auc_left_vs_right")
tables_dir  <- file.path(results_dir, "tables")

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

message("Running AUCell Left vs Right comparisons...")

# -------------------------------#
# Load LR-labeled cardiomyocytes
# -------------------------------#
# Update this path if you saved it somewhere else in 01_seurat_16ss.R
lr_obj_path <- here("results", "16ss", "heart16s_cluster3_LR.rds")

if (!file.exists(lr_obj_path)) {
  stop(
    paste0(
      "\nMissing input RDS: ", lr_obj_path,
      "\nPlease save the LR-labeled cardiomyocyte object in 01_seurat_16ss.R, e.g.:",
      "\n  saveRDS(heart16s_final_cluster3only, file = file.path(results_dir, 'heart16s_cluster3_LR.rds'))"
    ),
    call. = FALSE
  )
}

heart16s_final_cluster3only <- readRDS(lr_obj_path)

# Ensure identities are set + ordered
Idents(heart16s_final_cluster3only) <- factor(
  Idents(heart16s_final_cluster3only),
  levels = c("Lpositive", "Rnegative", "3")
)
print(table(Idents(heart16s_final_cluster3only)))

# -------------------------------#
# Helper function
# -------------------------------#
run_auc_lr <- function(
  seu,
  gene_set_name,
  gene_vector,
  label_for_plot,
  stage_label = "16 ss",
  colors_lr = c("Lpositive" = "#6998CA", "Rnegative" = "#F2A165")
) {
  # AUCell requires a ranking built from an expression matrix
  exprMatrix <- GetAssayData(seu, slot = "counts")

  geneSets <- list()
  geneSets[[gene_set_name]] <- gene_vector

  cells_rankings <- AUCell_buildRankings(exprMatrix, nCores = 1, plotStats = TRUE)
  cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)

  # Histogram (save)
  p_hist <- AUCell_plotHist(cells_AUC)

  ggsave(
    filename = file.path(figures_dir, paste0("AUCell_hist_", gene_set_name, "_", gsub(" ", "", stage_label), ".png")),
    plot = p_hist,
    width = 7, height = 5, dpi = 300
  )

  # Add AUC scores to object metadata
  auc_df <- as.data.frame(t(getAUC(cells_AUC)))
  seu <- AddMetaData(seu, auc_df)

  # Violin plot (Left vs Right only)
  p_vln <- VlnPlot(
    seu,
    features = gene_set_name,
    idents = c("Lpositive", "Rnegative"),
    cols = unname(colors_lr)
  ) +
    theme(
      legend.position = "none",
      plot.title = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x  = element_text(angle = 0, hjust = 0.5)
    ) +
    labs(y = label_for_plot) +
    scale_x_discrete(labels = c("Lpositive" = "Left", "Rnegative" = "Right")) +
    xlab(stage_label) +
    theme(
      axis.title.x = element_text(size = 12, face = "bold"),
      axis.title.y = element_text(size = 12)
    )

  ggsave(
    filename = file.path(figures_dir, paste0("AUCell_vln_", gene_set_name, "_", gsub(" ", "", stage_label), ".png")),
    plot = p_vln,
    width = 6, height = 4, dpi = 300
  )

  # Stats: Wilcoxon test (Right vs Left)
  stats_df <- auc_df
  stats_df$cluster <- Idents(seu)

  p_value <- wilcox.test(
    stats_df[stats_df$cluster == "Rnegative", gene_set_name],
    stats_df[stats_df$cluster == "Lpositive", gene_set_name]
  )$p.value

  list(
    updated_seu = seu,
    p_value = p_value,
    stats_df = stats_df
  )
}

# -------------------------------#
# 1) AUCell: Left-exclusive Nodal targets
# -------------------------------#
res_leftonly <- run_auc_lr(
  seu = heart16s_final_cluster3only,
  gene_set_name = "Active_Nodal_LeftOnly",
  gene_vector   = c("lft2", "pitx2", "ndr2", "lft1"),
  label_for_plot = "Nodal signaling\n(Left-exclusive target genes)",
  stage_label = "16 ss"
)

message("Wilcoxon p-value (Left-exclusive): ", signif(res_leftonly$p_value, 3))

# -------------------------------#
# 2) AUCell: All canonical Nodal targets
# -------------------------------#
res_all <- run_auc_lr(
  seu = res_leftonly$updated_seu, # carry forward metadata
  gene_set_name = "Active_Nodal",
  gene_vector   = c("lft2", "pitx2", "ndr2", "lft1", "acta1b", "has2", "bmp4"),
  label_for_plot = "Nodal signaling\n(All known target genes)",
  stage_label = "16 ss"
)

message("Wilcoxon p-value (All canonical): ", signif(res_all$p_value, 3))

# -------------------------------#
# Save stats table + updated object
# -------------------------------#
pvals_out <- data.frame(
  stage = "16 ss",
  gene_set = c("Active_Nodal_LeftOnly", "Active_Nodal"),
  p_value = c(res_leftonly$p_value, res_all$p_value)
)

write.csv(
  pvals_out,
  file = file.path(tables_dir, "AUCell_LR_Wilcoxon_pvalues.csv"),
  row.names = FALSE
)

saveRDS(
  res_all$updated_seu,
  file = file.path(results_dir, "heart16s_cluster3_LR_with_AUCell.rds")
)

message("05_auc_left_vs_right.R complete. Outputs saved to results/auc_left_vs_right and figures/auc_left_vs_right.")
