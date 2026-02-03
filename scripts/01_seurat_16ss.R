# 01_seurat_16ss.R
# Purpose:
#   - Standard Seurat workflow for 16 somite stage (16 ss)
#   - Cluster annotation + marker discovery
#   - Define Left/Right cardiomyocytes based on Nodal target expression
#   - AUCell scoring for pathway activity across clusters

# Setup
source(here("scripts", "00_setup.R"))
set.seed(123)

# Output folders
results_dir <- here("results", "16ss")
figures_dir <- here("figures", "16ss")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# Create Seurat object + QC
heart16s <- CreateSeuratObject(
  counts = heart16s.data,
  project = "16somdata",
  min.cells = 3,
  min.features = 200
)

heart16s[["percent.mt"]] <- PercentageFeatureSet(heart16s, pattern = "^mt-")

# QC summary
qc_summary <- list(
  nFeature_RNA_median = median(heart16s@meta.data$nFeature_RNA),
  nFeature_RNA_mad    = mad(heart16s@meta.data$nFeature_RNA, center = median(heart16s@meta.data$nFeature_RNA),
                            low = FALSE, high = FALSE),
  percent_mt_median   = median(heart16s@meta.data$percent.mt),
  percent_mt_mad      = mad(heart16s@meta.data$percent.mt, center = median(heart16s@meta.data$percent.mt),
                            low = FALSE, high = FALSE)
)
print(qc_summary)

# QC thresholds
heart16s_filtered <- subset(
  heart16s,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 9096 &
    percent.mt < 10.46 &
    nCount_RNA < 100000
)

# Save filtered object
saveRDS(heart16s_filtered, file = file.path(results_dir, "heart16s_filtered.rds"))

# Normalize + HVGs + Scale + PCA/UMAP/Clustering
heart16s_filtered <- NormalizeData(heart16s_filtered, normalization.method = "LogNormalize", scale.factor = 10000)
heart16s_filtered <- FindVariableFeatures(heart16s_filtered, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(heart16s_filtered)
heart16s_filtered <- ScaleData(heart16s_filtered, features = all.genes)

heart16s_filtered <- RunPCA(heart16s_filtered, features = VariableFeatures(object = heart16s_filtered))

# Dims/resolution
heart16s_filtered <- FindNeighbors(heart16s_filtered, dims = 1:13)
heart16s_filtered <- FindClusters(heart16s_filtered, resolution = 0.12)
heart16s_final <- RunUMAP(heart16s_filtered, dims = 1:13)

saveRDS(heart16s_final, file = file.path(results_dir, "heart16s_final.rds"))

# Cluster naming + UMAP plots
heart16s_final_clustersnamed <- heart16s_final

new_cluster_names <- c(
  "Mixed Cell Types", "Mesoderm", "Central Nervous System", "Cardiomyocyte",
  "Epidermis", "Somite", "Vasculature"
)
names(new_cluster_names) <- levels(heart16s_final_clustersnamed)
heart16s_final_clustersnamed <- RenameIdents(heart16s_final_clustersnamed, new_cluster_names)

p_umap_all <- DimPlot(heart16s_final_clustersnamed, reduction = "umap", pt.size = 0.5)

ggsave(
  filename = file.path(figures_dir, "16ss_umap_all_clusters.png"),
  plot = p_umap_all,
  width = 7, height = 5, dpi = 300
)

# Cardiac cluster highlight
p_umap_cardiac <- DimPlot(
  heart16s_final,
  reduction = "umap",
  cols = c("darkgrey", "darkgrey", "darkgrey", "#0a9f0a", "darkgrey", "darkgrey", "darkgrey", "darkgrey")
) + NoLegend()

ggsave(
  filename = file.path(figures_dir, "16ss_umap_cardiac_highlight.png"),
  plot = p_umap_cardiac,
  width = 7, height = 5, dpi = 300
)

# Cluster markers + heatmap
heart16s_final.markers <- FindAllMarkers(heart16s_final, only.pos = TRUE)
write.csv(
  heart16s_final.markers,
  file = file.path(results_dir, "cluster_markers_heart16ss.csv"),
  row.names = TRUE
)

top10 <- heart16s_final.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

p_heatmap <- DoHeatmap(heart16s_final, features = top10$gene) +
  NoLegend() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

ggsave(
  filename = file.path(figures_dir, "16ss_heatmap_top10_markers_per_cluster.png"),
  plot = p_heatmap,
  width = 10, height = 7, dpi = 300
)

# Violin plots: canonical cardiac markers
heart16s_final_reordered <- SetIdent(
  heart16s_final,
  value = factor(Idents(heart16s_final), levels = c("3", "0", "1", "2", "4", "5", "6"))
)

new_cluster_names_short <- c(
  "3" = "H",
  "0" = "MCT",
  "1" = "M",
  "2" = "CNS",
  "4" = "EP",
  "5" = "S",
  "6" = "V"
)
heart16s_final_reordered <- RenameIdents(heart16s_final_reordered, new_cluster_names_short)

marker_genes <- c("mef2ca", "myl7", "nkx2.5", "ttn.1")

for (g in marker_genes) {
  p_vln <- VlnPlot(
    heart16s_final_reordered,
    features = g,
    cols = c(
      "#0a9f0a", "darkgrey", "darkgrey", "darkgrey", "darkgrey",
      "darkgrey", "darkgrey", "darkgrey"
    )
  ) +
    NoLegend() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 8),
      plot.title  = element_text(face = "bold.italic")
    ) +
    xlab("Cluster")

  ggsave(
    filename = file.path(figures_dir, paste0("16ss_vln_", g, ".png")),
    plot = p_vln,
    width = 8, height = 4, dpi = 300
  )
}

# Left/Right assignment within cardiomyocytes
heart16s_final_cluster3only <- subset(heart16s_final, idents = c("3"))

# Left genes mean expression
pitx2_mean <- 0.1444
lft2_mean  <- 0.2881
lft1_mean  <- 0.07302
ndr2_mean  <- 0.06474

left_sided_genes_positive <- WhichCells(
  heart16s_final_cluster3only,
  expression = pitx2 > pitx2_mean | lft2 > lft2_mean | lft1 > lft1_mean | ndr2 > ndr2_mean
)
heart16s_final_cluster3only <- SetIdent(
  heart16s_final_cluster3only,
  cells = left_sided_genes_positive,
  value = "Lpositive"
)

left_sided_genes_negative <- WhichCells(
  heart16s_final_cluster3only,
  expression = lft2 == 0 & pitx2 == 0 & lft1 == 0 & ndr2 == 0
)
heart16s_final_cluster3only <- SetIdent(
  heart16s_final_cluster3only,
  cells = left_sided_genes_negative,
  value = "Rnegative"
)

# Save LR object
saveRDS(
  heart16s_final_cluster3only,
  file = file.path(results_dir, "heart16s_cardiomyocytes_LR_labeled.rds")
)

# LR differential expression
heart16s_LR_DEanalysis.markers <- FindMarkers(
  heart16s_final_cluster3only,
  ident.1 = "Lpositive",
  ident.2 = "Rnegative"
)

heart16s_LR_DEanalysis.markersFILTERED <- subset(
  heart16s_LR_DEanalysis.markers,
  p_val_adj < 0.05 & abs(avg_log2FC) > 0.25
)

write.csv(
  heart16s_LR_DEanalysis.markersFILTERED,
  file = file.path(results_dir, "LR_DEanalysis_heart16ss_filtered.csv")
)

# AUCell: pathway activity across clusters
exprMatrix <- GetAssayData(heart16s_final, slot = "counts")

geneSets <- list(
  Active_Nodal = c("lft2", "pitx2", "ndr2", "lft1")
)
# To measure BMP signaling instead:
# geneSets <- list(Active_BMP = c("nkx2.5", "hand2", "tbx20", "tbx2b"))
# To measure FGF signaling instead:
# geneSets <- list(Active_FGF = c("pea3", "erm", "gata4", "spry4"))

cells_rankings <- AUCell_buildRankings(exprMatrix, nCores = 1, plotStats = TRUE)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)

auc_df <- as.data.frame(t(getAUC(cells_AUC)))
heart16s_final <- AddMetaData(heart16s_final, auc_df)

# reorder again for plotting
heart16s_final_reordered <- SetIdent(
  heart16s_final,
  value = factor(Idents(heart16s_final), levels = c("3", "0", "1", "2", "4", "5", "6"))
)
heart16s_final_reordered <- RenameIdents(heart16s_final_reordered, new_cluster_names_short)

p_auc <- VlnPlot(
  heart16s_final_reordered,
  features = names(geneSets),
  cols = c("#0a9f0a", rep("darkgrey", 7))
) +
  theme(
    legend.position = "none",
    plot.title = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.title.y = element_text(size = 16)
  ) +
  labs(y = "Nodal signaling at 16 ss") +
  xlab("Cluster")

ggsave(
  filename = file.path(figures_dir, "AUCell_nodal_by_cluster.png"),
  plot = p_auc,
  width = 7, height = 5, dpi = 300
)

# Save updated object w/ AUCell metadata
saveRDS(heart16s_final, file = file.path(results_dir, "heart16s_final_with_AUCell.rds"))

message("01_seurat_16ss.R complete. Outputs saved to results/16ss and figures/16ss.")
