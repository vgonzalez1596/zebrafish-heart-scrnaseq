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


# -------------------------------#
# 1) Create Seurat object + QC
# -------------------------------#

heart18s <- CreateSeuratObject(
  counts = heart18s.data,
  project = "18somdata",
  min.cells = 3,
  min.features = 200
)

heart18s[["percent.mt"]] <- PercentageFeatureSet(heart18s, pattern = "^mt-")

# Optional QC summaries (kept for parity with your original code)
median(heart18s@meta.data$nFeature_RNA)
mad(heart18s@meta.data$nFeature_RNA, center = median(heart18s@meta.data$nFeature_RNA), low = FALSE, high = FALSE)
median(heart18s@meta.data$percent.mt)
mad(heart18s@meta.data$percent.mt, center = median(heart18s@meta.data$percent.mt), low = FALSE, high = FALSE)

# Filter thresholds (from your original analysis)
heart18s_filtered <- subset(
  heart18s,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 7974 &
    percent.mt < 7.79 &
    nCount_RNA < 100000
)

saveRDS(heart18s_filtered, file.path(rds_dir, "heart18s_filtered.rds"))

# -------------------------------#
# 2) Normalize + HVGs + Scale + PCA
# -------------------------------#

heart18s_norm <- NormalizeData(
  heart18s_filtered,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

heart18s_norm <- FindVariableFeatures(
  heart18s_norm,
  selection.method = "vst",
  nfeatures = 2000
)

top10 <- head(VariableFeatures(heart18s_norm), 10)

all.genes <- rownames(heart18s_norm)
heart18s_norm <- ScaleData(heart18s_norm, features = all.genes)

heart18s_norm <- RunPCA(heart18s_norm, features = VariableFeatures(object = heart18s_norm))

saveRDS(heart18s_norm, file.path(rds_dir, "heart18s_norm_pca.rds"))

# -------------------------------#
# 3) Neighbors + Clusters + UMAP
# -------------------------------#

heart18s_nn <- FindNeighbors(heart18s_norm, dims = 1:16)
heart18s_nn <- FindClusters(heart18s_nn, resolution = 0.15)
heart18s_final <- RunUMAP(heart18s_nn, dims = 1:16)

saveRDS(heart18s_final, file.path(rds_dir, "heart18s_final_umap.rds"))

# UMAP: all clusters (named)
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
  filename = file.path(fig_dir, "18ss_umap_all_clusters_named.png"),
  plot = p_umap_all,
  width = 7, height = 5, dpi = 300
)

# UMAP: highlight cardiac cluster in green, others grey
# (matches your original colors vector ordering)
p_umap_heart <- DimPlot(
  heart18s_final,
  reduction = "umap",
  cols = c(
    "darkgrey", "darkgrey", "#0a9f0a", "darkgrey", "darkgrey",
    "darkgrey", "darkgrey", "darkgrey", "darkgrey", "darkgrey"
  )
) + NoLegend()

ggsave(
  filename = file.path(fig_dir, "18ss_umap_cardiac_highlight.png"),
  plot = p_umap_heart,
  width = 7, height = 5, dpi = 300
)

# -------------------------------#
# 4) Cluster markers + heatmap
# -------------------------------#

heart18s_final.markers <- FindAllMarkers(heart18s_final, only.pos = TRUE)

write.csv(
  heart18s_final.markers,
  file = file.path(tbl_dir, "Clustermarkers_heart18ss.csv"),
  row.names = TRUE
)

top10_markers <- heart18s_final.markers %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = 10, wt = avg_log2FC)

p_heat <- DoHeatmap(heart18s_final, features = top10_markers$gene) +
  NoLegend() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

ggsave(
  filename = file.path(fig_dir, "18ss_heatmap_top10_markers_per_cluster.png"),
  plot = p_heat,
  width = 10, height = 7, dpi = 300
)

# -------------------------------#
# 5) Canonical cardiac marker violins
# -------------------------------#

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
    filename = file.path(fig_dir, paste0("18ss_vln_", g, ".png")),
    plot = p_vln,
    width = 8, height = 4, dpi = 300
  )
}

# -------------------------------#
# 6) Assign left/right within cardiac cluster + LR DE
# -------------------------------#

# Cardiac cluster is "2" in 18 ss (per your original)
heart18s_final_cluster2only <- subset(heart18s_final, idents = c("2"))

# Optional expression summaries (parity with original)
summary(FetchData(object = heart18s_final_cluster2only, vars = "pitx2"))
summary(FetchData(object = heart18s_final_cluster2only, vars = "lft2"))
summary(FetchData(object = heart18s_final_cluster2only, vars = "lft1"))
summary(FetchData(object = heart18s_final_cluster2only, vars = "ndr2"))

# Thresholds from your original mean-expression notes
# Mean expression: lft1 = 0.3455, lft2 = 0.9976, ndr2 = 0.2029, pitx2 = 0.3872
left_pos <- WhichCells(
  heart18s_final_cluster2only,
  expression = pitx2 > 0.3872 | lft2 > 0.9976 | lft1 > 0.3455 | ndr2 > 0.2029
)
heart18s_final_cluster2only <- SetIdent(
  heart18s_final_cluster2only,
  cells = left_pos,
  value = "Lpositive"
)

right_neg <- WhichCells(
  heart18s_final_cluster2only,
  expression = lft2 == 0 & pitx2 == 0 & lft1 == 0 & ndr2 == 0
)
heart18s_final_cluster2only <- SetIdent(
  heart18s_final_cluster2only,
  cells = right_neg,
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
  file = file.path(tbl_dir, "LR_DEanalysis_heart18ss.csv"),
  row.names = TRUE
)

saveRDS(
  heart18s_final_cluster2only,
  file.path(rds_dir, "heart18s_cardiac_cluster_LR_labeled.rds")
)

message("LR label counts:")
print(table(Idents(heart18s_final_cluster2only)))

# -------------------------------#
# 7) AUCell signaling across all clusters (18 ss)
# -------------------------------#

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

# Reorder + rename clusters (same as earlier) for plotting
heart18s_auc_reordered <- SetIdent(
  heart18s_final,
  value = factor(Idents(heart18s_final), levels = c("2", "0", "1", "3", "4", "5", "6", "7", "8", "9"))
)
heart18s_auc_reordered <- RenameIdents(heart18s_auc_reordered, new_cluster_names_short)

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
  filename = file.path(fig_dir, "18ss_aucell_nodal_by_cluster.png"),
  plot = p_auc,
  width = 8, height = 4.5, dpi = 300
)

saveRDS(heart18s_final, file.path(rds_dir, "heart18s_final_with_AUCell.rds"))

message("02_seurat_18ss.R complete. Outputs saved to: ", out_dir)

