# 02_seurat_18ss.R
# Purpose: 
#   - Standard Seurat workflow for 18 somite stage (18 ss)
#   - Cluster annotation + marker discovery
#   - Define Left/Right cardiomyocytes based on Nodal target expression
#   - AUCell scoring for pathway activity across clusters

# Setup
source(here("scripts", "00_setup.R"))
set.seed(123)

# Output folders
results_dir <- here("results", "18ss")
figures_dir <- here("figures", "18ss")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# Create Seurat object + QC
heart18s <- CreateSeuratObject(
  counts = heart18s.data,
  project = "18somdata",
  min.cells = 3,
  min.features = 200
)

heart18s[["percent.mt"]] <- PercentageFeatureSet(heart18s, pattern = "^mt-")

# QC summary
qc_summary <- list(
  nFeature_RNA_median = median(heart18s@meta.data$nFeature_RNA),
  nFeature_RNA_mad    = mad(heart18s@meta.data$nFeature_RNA, center = median(heart18s@meta.data$nFeature_RNA),
                            low = FALSE, high = FALSE),
  percent_mt_median   = median(heart18s@meta.data$percent.mt),
  percent_mt_mad      = mad(heart18s@meta.data$percent.mt, center = median(heart18s@meta.data$percent.mt),
                            low = FALSE, high = FALSE)
)
print(qc_summary)

# QC thresholds
heart18s_filtered <- subset(
  heart18s,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 7974 &
    percent.mt < 7.79 &
    nCount_RNA < 100000
)

# Save filtered object
saveRDS(heart18s_filtered, file = file.path(results_dir, "heart18s_filtered.rds"))

# Normalize + HVGs + Scale + PCA/UMAP/Clustering
heart18s_filtered <- NormalizeData(heart18s_filtered, normalization.method = "LogNormalize", scale.factor = 10000)
heart18s_filtered <- FindVariableFeatures(heart18s_filtered, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(heart18s_filtered)
heart18s_filtered <- ScaleData(heart18s_filtered, features = all.genes)

heart18s_filtered <- RunPCA(heart18s_filtered, features = VariableFeatures(object = heart18s_filtered))

# Dims/resolution
heart18s_filtered <- FindNeighbors(heart18s_filtered, dims = 1:16)
heart18s_filtered <- FindClusters(heart18s_filtered, resolution = 0.15)
heart18s_final <- RunUMAP(heart18s_filtered, dims = 1:16)

saveRDS(heart18s_final, file = file.path(results_dir, "heart18s_final.rds"))

# Cluster naming + UMAP plots
heart18s_final_clustersnamed <- heart18s_final

new_cluster_names <- c(
  "Mixed Cell Types",
  "Pharyngeal Arch",
  "Cardiomyocyte",
  "Brain",
  "Neural Crest",
  "Epidermis",
  "Somite",
  "Hatching Gland",
  "Vasculature",
  "Blood Cells"
)
names(new_cluster_names) <- levels(heart18s_final_clustersnamed)
heart18s_final_clustersnamed <- RenameIdents(heart18s_final_clustersnamed, new_cluster_names)

p_umap_all <- DimPlot(heart18s_final_clustersnamed, reduction = "umap", pt.size = 0.5)

ggsave(
  filename = file.path(figures_dir, "18ss_umap_all_clusters.png"),
  plot = p_umap_all,
  width = 7, height = 5, dpi = 300
)

# Cardiac cluster highlight
p_umap_cardiac <- DimPlot(
  heart18s_final,
  reduction = "umap",
  cols = c(
    "darkgrey", "darkgrey", "#0a9f0a", "darkgrey", "darkgrey",
    "darkgrey", "darkgrey", "darkgrey", "darkgrey", "darkgrey")
) + NoLegend()

ggsave(
  filename = file.path(figures_dir, "18ss_umap_cardiac_highlight.png"),
  plot = p_umap_cardiac,
  width = 7, height = 5, dpi = 300
)

# Cluster markers + heatmap
heart18s_final.markers <- FindAllMarkers(heart18s_final, only.pos = TRUE)
write.csv(
  heart18s_final.markers,
  file = file.path(results_dir, "cluster_markers_heart18ss.csv"),
  row.names = TRUE
)

top10 <- heart18s_final.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

p_heatmap <- DoHeatmap(heart18s_final, features = top10$gene) +
  NoLegend() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

ggsave(
  filename = file.path(figures_dir, "18ss_heatmap_top10_markers_per_cluster.png"),
  plot = p_heatmap,
  width = 10, height = 7, dpi = 300
)

# Violin plots: canonical cardiac markers
heart18s_final_reordered <- SetIdent(
  heart18s_final,
  value = factor(Idents(heart18s_final), levels = c("2", "0", "1", "3", "4", "5", "6", "7", "8", "9"))
)

new_cluster_names_short <- c(
  "2" = "H",
  "0" = "MCT",
  "1" = "PA",
  "3" = "CNS",
  "4" = "NC",
  "5" = "EP",
  "6" = "S",
  "7" = "HG",
  "8" = "V",
  "9" = "BC"
)
heart18s_final_reordered <- RenameIdents(heart18s_final_reordered, new_cluster_names_short)

marker_genes <- c("mef2ca", "myl7", "nkx2.5", "ttn.1")

for (g in marker_genes) {
  p_vln <- VlnPlot(
    heart18s_final_reordered,
    features = g,
    cols = c(
      "#0a9f0a", "darkgrey", "darkgrey", "darkgrey", "darkgrey",
      "darkgrey", "darkgrey", "darkgrey", "darkgrey", "darkgrey"
    )
  ) +
    NoLegend() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 8),
      plot.title  = element_text(face = "bold.italic")
    ) +
    xlab("Cluster")

  ggsave(
    filename = file.path(figures_dir, paste0("18ss_vln_", g, ".png")),
    plot = p_vln,
    width = 8, height = 4, dpi = 300
  )
}

# Left/Right assignment within cardiomyocytes
heart18s_final_cluster2only <- subset(heart18s_final, idents = c("2"))

summary(FetchData(object = heart18s_final_cluster2only, vars = "pitx2"))
summary(FetchData(object = heart18s_final_cluster2only, vars = "lft2"))
summary(FetchData(object = heart18s_final_cluster2only, vars = "lft1"))
summary(FetchData(object = heart18s_final_cluster2only, vars = "ndr2"))
# Mean expression: lft1 = 0.3455, lft2 = 0.9976, ndr2 = 0.2029, pitx2 = 0.3872

left_sided_genes_positive <- WhichCells(
  heart18s_final_cluster2only,
  expression = pitx2 > 0.3872 | lft2 > 0.9976 | lft1 > 0.3455 | ndr2 > 0.2029
)
heart18s_final_cluster2only <- SetIdent(
  heart18s_final_cluster2only,
  cells = left_sided_genes_positive,
  value = "Lpositive"
)

left_sided_genes_negative <- WhichCells(
  heart18s_final_cluster2only,
  expression = lft2 == 0 & pitx2 == 0 & lft1 == 0 & ndr2 == 0
)
heart18s_final_cluster2only <- SetIdent(
  heart18s_final_cluster2only,
  cells = left_sided_genes_negative,
  value = "Rnegative"
)

# LR differential expression
heart18s_LR_DEanalysis.markers <- FindMarkers(
  heart18s_final_cluster2only,
  ident.1 = "Lpositive",
  ident.2 = "Rnegative"
)

heart18s_LR_DEanalysis.markersFILTERED <- subset(
  heart18s_LR_DEanalysis.markers,
  p_val_adj < 0.05 & abs(avg_log2FC) > 0.25
)

write.csv(
  heart18s_LR_DEanalysis.markersFILTERED,
  file = file.path(results_dir, "LR_DEanalysis_heart18ss_filtered.csv"),
  row.names = TRUE
)

# AUCell: pathway activity across clusters
exprMatrix <- GetAssayData(heart18s_final, slot = "counts")

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
heart18s_final <- AddMetaData(heart18s_final, auc_df)

# Reorder + rename clusters for plotting
heart18s_final_reordered <- SetIdent(
  heart18s_final,
  value = factor(Idents(heart18s_final), levels = c("2", "0", "1", "3", "4", "5", "6", "7", "8", "9"))
)
heart18s_final_reordered <- RenameIdents(heart18s_final_reordered, new_cluster_names_short)

p_auc <- VlnPlot(
  heart18s_auc_reordered,
  features = names(geneSets),
  cols = c(
    "#0a9f0a", "darkgrey", "darkgrey", "darkgrey", "darkgrey",
    "darkgrey", "darkgrey", "darkgrey", "darkgrey", "darkgrey"
  )
) +
  theme(
    legend.position = "none",
    plot.title = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
    axis.title.y = element_text(size = 16)
  ) +
  labs(y = "Nodal signaling at 18 ss") +
  xlab("Cluster")

ggsave(
  filename = file.path(figures_dir, "18ss_AUCell_Nodal_by_cluster.png"),
  plot = p_auc,
  width = 7, height = 5, dpi = 300
)

saveRDS(heart18s_final, file = file.path(results_dir, "heart18s_final_with_AUCell.rds"))

message("02_seurat_18ss.R complete. Outputs saved to results/18ss and figures/18ss.")
