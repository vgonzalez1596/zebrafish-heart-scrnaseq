# 03_seurat_20ss.R
# Purpose: standard Seurat workflow for 20 somite stage (20 ss) + LR assignment + DE + AUCell signaling
#
# Expects: 00_setup.R has been run (or sourced) and created:
#   - heart20s.data  (Read10X matrix)
#
# Outputs (created by this script):
#   results/20ss/            (tables + R objects)
#   figures/20ss/            (plots)

# -------------------------------#
# 0) Setup
# -------------------------------#

source(here("scripts", "00_setup.R"))  # adjust if your repo uses a different path
set.seed(123)

# Output folders (match 01)
results_dir <- here("results", "20ss")
figures_dir <- here("figures", "20ss")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

message("Running 20 ss Seurat workflow...")

# Small helper for saving ggplots consistently
save_plot <- function(p, filename, width = 7, height = 5) {
  ggsave(
    filename = file.path(figures_dir, filename),
    plot = p,
    width = width,
    height = height,
    units = "in",
    dpi = 300
  )
}

# -------------------------------#
# 1) Create Seurat object + QC
# -------------------------------#

heart20s <- CreateSeuratObject(
  counts = heart20s.data,
  project = "20somdata",
  min.cells = 3,
  min.features = 200
)

heart20s[["percent.mt"]] <- PercentageFeatureSet(heart20s, pattern = "^mt-")

# (Optional) QC summaries (kept for reproducibility/logging)
qc_summary <- list(
  nFeature_RNA_median = median(heart20s@meta.data$nFeature_RNA),
  nFeature_RNA_mad    = mad(heart20s@meta.data$nFeature_RNA,
                            center = median(heart20s@meta.data$nFeature_RNA),
                            low = FALSE, high = FALSE),
  percent_mt_median   = median(heart20s@meta.data$percent.mt),
  percent_mt_mad      = mad(heart20s@meta.data$percent.mt,
                            center = median(heart20s@meta.data$percent.mt),
                            low = FALSE, high = FALSE)
)
saveRDS(qc_summary, file.path(results_dir, "qc_summary_20ss.rds"))

# Filter (thresholds from your original notebook)
heart20s_filtered <- subset(
  heart20s,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 5820 &
    percent.mt < 5.22 &
    nCount_RNA < 100000
)

saveRDS(heart20s_filtered, file.path(results_dir, "heart20s_filtered.rds"))

# -------------------------------#
# 2) Normalize + HVGs + Scale + PCA + Neighbors/Clusters + UMAP
# -------------------------------#

heart20s_norm <- NormalizeData(
  heart20s_filtered,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

heart20s_norm <- FindVariableFeatures(
  heart20s_norm,
  selection.method = "vst",
  nfeatures = 2000
)

top10_hvgs <- head(VariableFeatures(heart20s_norm), 10)
write.csv(
  data.frame(top10_hvgs = top10_hvgs),
  file.path(results_dir, "top10_variable_features_20ss.csv"),
  row.names = FALSE
)

all.genes <- rownames(heart20s_norm)
heart20s_norm <- ScaleData(heart20s_norm, features = all.genes)

heart20s_norm <- RunPCA(heart20s_norm, features = VariableFeatures(object = heart20s_norm))

# dims + resolution from your original notebook
heart20s_norm <- FindNeighbors(heart20s_norm, dims = 1:13)
heart20s_norm <- FindClusters(heart20s_norm, resolution = 0.2)

heart20s_final <- RunUMAP(heart20s_norm, dims = 1:13)

saveRDS(heart20s_final, file.path(results_dir, "heart20s_final.rds"))

# -------------------------------#
# 3) UMAPs + cluster naming
# -------------------------------#

# UMAP with all clusters, named
heart20s_final_clustersnamed <- heart20s_final
new_cluster_names <- c(
  "Mixed Cell Types",
  "Mesoderm",
  "Macrophage",
  "Cardiomyocyte",
  "Central Nervous System",
  "Neural Crest",
  "Somite",
  "Ectoderm"
)
names(new_cluster_names) <- levels(heart20s_final_clustersnamed)
heart20s_final_clustersnamed <- RenameIdents(heart20s_final_clustersnamed, new_cluster_names)

p_umap_all <- DimPlot(heart20s_final_clustersnamed, reduction = "umap", pt.size = 0.5)
save_plot(p_umap_all, "umap_all_clusters_20ss.png", width = 7.5, height = 6)

# UMAP highlighting cardiac cluster (cluster "3" in your original object) in green
p_umap_heart <- DimPlot(
  heart20s_final,
  reduction = "umap",
  cols = c("darkgrey", "darkgrey", "darkgrey", "#0a9f0a", "darkgrey", "darkgrey", "darkgrey", "darkgrey")
) + NoLegend()
save_plot(p_umap_heart, "umap_heart_highlight_20ss.png", width = 7.5, height = 6)

# -------------------------------#
# 4) Cluster markers + heatmap
# -------------------------------#

heart20s_final_markers <- FindAllMarkers(heart20s_final, only.pos = TRUE)
write.csv(
  heart20s_final_markers,
  file.path(results_dir, "cluster_markers_heart20ss.csv"),
  row.names = FALSE
)

top10 <- heart20s_final_markers %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = 10, wt = avg_log2FC)

p_heat <- DoHeatmap(heart20s_final, features = top10$gene) +
  NoLegend() +
  theme(strip.text = element_text(angle = 0))

save_plot(p_heat, "heatmap_top10_markers_20ss.png", width = 10, height = 8)

# -------------------------------#
# 5) Cardiac marker violins + “uncharacterized” gene panels
# -------------------------------#

# Reorder + rename clusters for violin plots (match your notebook)
heart20s_final_reordered <- SetIdent(
  heart20s_final,
  value = factor(Idents(heart20s_final), levels = c("3", "0", "1", "2", "4", "5", "6", "7"))
)

new_cluster_names_short <- c(
  "3" = "H",
  "0" = "MCT",
  "1" = "M",
  "2" = "MA",
  "4" = "CNS",
  "5" = "NC",
  "6" = "S",
  "7" = "EC"
)
heart20s_final_reordered <- RenameIdents(heart20s_final_reordered, new_cluster_names_short)

# Canonical cardiac markers
cardiac_markers <- c("mef2ca", "myl7", "nkx2.5", "ttn.1")
for (g in cardiac_markers) {
  p <- VlnPlot(
    heart20s_final_reordered,
    features = g,
    cols = c("#0a9f0a", rep("darkgrey", 8))
  ) +
    NoLegend() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
      plot.title = element_text(face = "bold.italic")
    ) +
    xlab("Cluster")

  save_plot(p, paste0("vln_", g, "_20ss.png"), width = 7.5, height = 4.5)
}

# Uncharacterized/novel expression panels (as in your notebook)
unknown_sets <- list(
  unknown_function = c("si:ch211-131k2.3", "si:ch211-221f10.2", "si:dkey-184p18.2"),
  unknown_tx_regs  = c("dcaf6", "mxd1", "nfat5b"),
  unknown_cytoskel = c("coro6", "fhod3a", "lmod2b")
)

for (set_name in names(unknown_sets)) {
  genes <- unknown_sets[[set_name]]

  plots <- VlnPlot(
    heart20s_final_reordered,
    features = genes,
    cols = c("#0a9f0a", rep("darkgrey", 7)),
    combine = FALSE
  )

  plots <- lapply(plots, function(p) {
    p + NoLegend() +
      theme(
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
        plot.title = element_text(face = "bold.italic")
      ) +
      xlab("Cluster")
  })

  p_panel <- wrap_plots(plots)
  save_plot(p_panel, paste0("vln_panel_", set_name, "_20ss.png"), width = 10, height = 4)
}

# -------------------------------#
# 6) Left/Right assignment + DE in cardiac cluster
# -------------------------------#

# Cardiac cluster is "3" in your original 20ss notebook
heart20s_final_cluster3only <- subset(heart20s_final, idents = c("3"))

# (Optional) inspect expression distributions
expr_summaries <- list(
  pitx2 = summary(FetchData(heart20s_final_cluster3only, vars = "pitx2")[,1]),
  lft2  = summary(FetchData(heart20s_final_cluster3only, vars = "lft2")[,1]),
  lft1  = summary(FetchData(heart20s_final_cluster3only, vars = "lft1")[,1]),
  ndr2  = summary(FetchData(heart20s_final_cluster3only, vars = "ndr2")[,1])
)
saveRDS(expr_summaries, file.path(results_dir, "lr_gene_expression_summaries_20ss.rds"))

# Thresholds from your notebook:
# Mean expression: lft1 = 0.07007, lft2 = 1.760, ndr2 = 0.09333, pitx2 = 0.1417
left_sided_genes_positive <- WhichCells(
  heart20s_final_cluster3only,
  expression = pitx2 > 0.1417 | lft2 > 1.760 | lft1 > 0.07007 | ndr2 > 0.09333
)
heart20s_final_cluster3only <- SetIdent(
  heart20s_final_cluster3only,
  cells = left_sided_genes_positive,
  value = "Lpositive"
)

left_sided_genes_negative <- WhichCells(
  heart20s_final_cluster3only,
  expression = lft2 == 0 & pitx2 == 0 & lft1 == 0 & ndr2 == 0
)
heart20s_final_cluster3only <- SetIdent(
  heart20s_final_cluster3only,
  cells = left_sided_genes_negative,
  value = "Rnegative"
)

# Sanity check table
lr_table <- table(Idents(heart20s_final_cluster3only))
write.csv(
  as.data.frame(lr_table),
  file.path(results_dir, "lr_identity_counts_20ss.csv"),
  row.names = FALSE
)

# DE: Left vs Right
heart20s_LR_DEanalysis_markers <- FindMarkers(
  heart20s_final_cluster3only,
  ident.1 = "Lpositive",
  ident.2 = "Rnegative"
)

heart20s_LR_DEanalysis_markersFILTERED <- subset(
  heart20s_LR_DEanalysis_markers,
  p_val_adj < 0.05 & abs(avg_log2FC) > 0.25
)

write.csv(
  heart20s_LR_DEanalysis_markersFILTERED,
  file.path(results_dir, "LR_DEanalysis_heart20ss_filtered.csv"),
  row.names = TRUE
)

saveRDS(heart20s_final_cluster3only, file.path(results_dir, "heart20s_final_cluster3only_LR.rds"))

# -------------------------------#
# 7) AUCell signaling across clusters (20 ss)
# -------------------------------#

exprMatrix <- GetAssayData(heart20s_final, slot = "counts")

# Default: Nodal (as in notebook)
geneSets <- list(Active_Nodal = c("lft2", "pitx2", "ndr2", "lft1"))

# To measure BMP signaling instead:
# geneSets <- list(Active_BMP = c("nkx2.5", "hand2", "tbx20", "tbx2b"))
#
# To measure FGF signaling instead:
# geneSets <- list(Active_FGF = c("pea3", "erm", "gata4", "spry4"))

cells_rankings <- AUCell_buildRankings(exprMatrix, nCores = 1, plotStats = TRUE)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)

auc_df <- as.data.frame(t(getAUC(cells_AUC)))
heart20s_final_auc <- AddMetaData(heart20s_final, auc_df)

# Reorder + rename for plotting (same mapping as above)
heart20s_final_auc_reordered <- SetIdent(
  heart20s_final_auc,
  value = factor(Idents(heart20s_final_auc), levels = c("3", "0", "1", "2", "4", "5", "6", "7"))
)
heart20s_final_auc_reordered <- RenameIdents(heart20s_final_auc_reordered, new_cluster_names_short)

p_auc <- VlnPlot(
  heart20s_final_auc_reordered,
  features = names(geneSets),
  cols = c("#0a9f0a", rep("darkgrey", 7))
) +
  theme(
    legend.position = "none",
    plot.title = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.title.y = element_text(size = 16)
  ) +
  labs(y = "Nodal signaling at 20 ss") +
  xlab("Cluster")

save_plot(p_auc, "aucell_nodal_by_cluster_20ss.png", width = 7.5, height = 5)

saveRDS(heart20s_final_auc, file.path(results_dir, "heart20s_final_with_AUCell_metadata.rds"))

message("Done! Outputs written to: ", results_dir, " and ", figures_dir)

