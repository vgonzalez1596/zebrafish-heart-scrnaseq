# 05_auc_left_vs_right.R
# Purpose: AUCell analyses of Nodal signaling in left vs right cardiomyocytes at 16 ss, 18 ss, and 20 ss.

# Setup
source(here("scripts", "00_setup.R"))
set.seed(123)

results_dir <- here("results", "auc_left_vs_right")
figures_dir <- here("figures", "auc_left_vs_right")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# Gene sets
gene_sets <- list(
  Active_Nodal_LeftOnly = c("lft2", "pitx2", "ndr2", "lft1"),
  Active_Nodal_AllCanonical = c(
    "lft2", "pitx2", "ndr2", "lft1",
    "acta1b", "has2", "bmp4"
  )
)

# AUCell + plot + stats
run_auc_lr <- function(seu_lr, gene_set_name, gene_set, stage_label, prefix) {

  Idents(seu_lr) <- factor(Idents(seu_lr),
                            levels = c("Lpositive", "Rnegative"))

  exprMatrix <- GetAssayData(seu_lr, slot = "counts")

  cells_rankings <- AUCell_buildRankings(exprMatrix, nCores = 1,
                                          plotStats = FALSE)

  cells_AUC <- AUCell_calcAUC(
    list(setNames(list(gene_set), gene_set_name)),
    cells_rankings
  )

  auc_df <- as.data.frame(t(getAUC(cells_AUC)))
  seu_lr <- AddMetaData(seu_lr, auc_df)

  # Plot label
  ylab_txt <- if (gene_set_name == "Active_Nodal_LeftOnly") {
    "Nodal signaling\n(Left-exclusive targets)"
  } else {
    "Nodal signaling\n(All canonical targets)"
  }

  p_vln <- VlnPlot(
    seu_lr,
    features = gene_set_name,
    idents = c("Lpositive", "Rnegative"),
    cols = c("#6998CA", "#F2A165")
  ) +
    theme(
      legend.position = "none",
      plot.title = element_blank(),
      axis.text.x = element_text(angle = 0, hjust = 0.5)
    ) +
    labs(y = ylab_txt) +
    scale_x_discrete(labels = c("Lpositive" = "Left",
                                 "Rnegative" = "Right")) +
    xlab(stage_label)

  ggsave(
    file.path(figures_dir,
              paste0(prefix, "_", gene_set_name, ".png")),
    p_vln, width = 6, height = 4, dpi = 300
  )

  # Wilcoxon test
  auc_df$cluster <- Idents(seu_lr)

  p_value <- wilcox.test(
    auc_df[auc_df$cluster == "Rnegative", gene_set_name],
    auc_df[auc_df$cluster == "Lpositive", gene_set_name]
  )$p.value

  data.frame(
    stage = stage_label,
    gene_set = gene_set_name,
    n_left = sum(auc_df$cluster == "Lpositive"),
    n_right = sum(auc_df$cluster == "Rnegative"),
    wilcox_p_value = p_value,
    stringsAsFactors = FALSE
  )
}

# Load LR objects
datasets <- list(
  `16 ss` = list(
    obj = readRDS(here("results","16ss",
                       "heart16s_cluster3_LR.rds")),
    prefix = "16ss"
  ),
  `18 ss` = list(
    obj = readRDS(here("results","18ss",
                       "heart18s_cluster2_LR.rds")),
    prefix = "18ss"
  ),
  `20 ss` = list(
    obj = readRDS(here("results","20ss",
                       "heart20s_cluster3_LR.rds")),
    prefix = "20ss"
  )
)

# Run AUCell for all datasets
stats_all <- list()

for (stage_label in names(datasets)) {

  seu_lr <- datasets[[stage_label]]$obj
  prefix <- datasets[[stage_label]]$prefix

  for (gs_name in names(gene_sets)) {

    stats_row <- run_auc_lr(
      seu_lr,
      gs_name,
      gene_sets[[gs_name]],
      stage_label,
      prefix
    )

    stats_all[[paste(stage_label, gs_name)]] <- stats_row
  }
}

stats_df <- do.call(rbind, stats_all)

write.csv(
  stats_df,
  file.path(results_dir,
            "auc_left_vs_right_wilcox_stats.csv"),
  row.names = FALSE
)

print(stats_df)

message("05_auc_left_vs_right.R complete. Outputs saved to results/auc_left_vs_right and figures/auc_left_vs_right.")
