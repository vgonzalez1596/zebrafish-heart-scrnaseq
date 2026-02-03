# 03_seurat_20ss.R
# Purpose:
#   - Standard Seurat workflow for 20 somite stage (20 ss)
#   - Cluster annotation + marker discovery
#   - Define Left/Right cardiomyocytes based on Nodal target expression
#   - AUCell scoring for pathway activity across clusters

# Setup
source(here("scripts", "00_setup.R"))
set.seed(123)

# Output folders
results_dir <- here("results", "20ss")
figures_dir <- here("figures", "20ss")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# Create Seurat object + QC
heart20s <- CreateSeuratObject(
  counts = heart20s.data,
  project = "20somdata",
  min.cells = 3,
  min.features = 200
)

heart20s[["percent.mt"]] <- PercentageFeatureSet(heart20s, pattern = "^mt-")

# QC summary
qc_summary <- list(
  nFeature_RNA_median = median(heart20s@meta.data$nFeature_RNA),
  nFeature_RNA_mad    = mad(heart20s@meta.data$nFeature_RNA, center = median(heart20s@meta.data$nFeature_RNA),
                            low = FALSE, high = FALSE),
  percent_mt_median   = median(heart20s@meta.data$percent.mt),
  percent_mt_mad      = mad(heart20s@meta.data$percent.mt, center = median(heart20s@meta.data$percent.mt),
                            low = FALSE, high = FALSE)
)
print(qc_summary)

# QC thresholds
heart20s_filtered <- subset(
  heart20s,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 5820 &
    percent.mt < 5.22 &
    nCount_RNA < 100000
)

# Save filtered object
saveRDS(heart20s_filtered, file = file.path(results_dir, "heart20s_filtered.rds"))

# Normalize + HVGs + Scale + PCA/UMAP/Clustering
heart20s_filtered <- NormalizeData(heart20s_filtered, normalization.method = "LogNormalize", scale.factor = 10000)
heart20s_filtered <- FindVariableFeatures(heart20s_filtered, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(heart20s_filtered)
heart20s_filtered <- ScaleData(heart20s_filtered, features = all.genes)

heart20s_filtered <- RunPCA(heart20s_filtered, features = VariableFeatures(object = heart20s_filtered))

# Dims/resolution
heart20s_filtered <- FindNeighbors(heart20s_filtered, dims = 1:13)
heart20s_filtered <- FindClusters(heart20s_filtered, resolution = 0.2)
heart20s_final <- RunUMAP(heart20s_filtered, dims = 1:13)

saveRDS(heart20s_final, file = file.path(results_dir, "heart20s_final.rds"))

# Cluster naming + UMAP plots
heart20s_final_clustersnamed <- heart20s_final

new_cluster_names <- c(
  "Mixed Cell Types", "Mesoderm", "Macrophage", "Cardiomyocyte",
  "Central Nervous System", "Neural Crest", "Somite", "Ectoderm"
)
names(new_cluster_names) <- levels(heart20s_final_clustersnamed)
heart20s_final_clustersnamed <- RenameIdents(heart20s_final_clustersnamed, new_cluster_names)

p_umap_all <- DimPlot(heart20s_final_clustersnamed, reduction = "umap", pt.size = 0.5)

ggsave(
  filename = file.path(figures_dir, "20ss_umap_all_clusters.png"),
  plot = p_umap_all,
  width = 7, height = 5, dpi = 300
)

# Cardiac cluster highlight
p_umap_cardiac <- DimPlot(
  heart20s_final,
  reduction = "umap",
  cols = c("darkgrey", "darkgrey", "darkgrey", "#0a9f0a", "darkgrey", "darkgrey", "darkgrey", "darkgrey")
) + NoLegend()

ggsave(
  filename = file.path(figures_dir, "20ss_umap_cardiac_highlight.png"),
  plot = p_umap_cardiac,
  width = 7, height = 5, dpi = 300
)

# Cluster markers + heatmap
heart20s_final.markers <- FindAllMarkers(heart20s_final, only.pos = TRUE)
write.csv(
  heart20s_final.markers,
  file = file.path(results_dir, "cluster_markers_heart20ss.csv"),
  row.names = TRUE
)

top10 <- heart20s_final.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

p_heatmap <- DoHeatmap(heart20s_final, features = top10$gene) +
  NoLegend() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

ggsave(
  filename = file.path(figures_dir, "20ss_heatmap_top10_markers_per_cluster.png"),
  plot = p_heatmap,
  width = 10, height = 7, dpi = 300
)

# Violin plots: canonical cardiac markers
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

marker_genes <- c("mef2ca", "myl7", "nkx2.5", "ttn.1")

for (g in marker_genes) {
  p_vln <- VlnPlot(
    heart20s_final_reordered,
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
    filename = file.path(figures_dir, paste0("20ss_vln_", g, ".png")),
    plot = p_vln,
    width = 8, height = 4, dpi = 300
  )
}

# Novel/uncharacterized genes
gene_sets <- list(
  unknown_function = c("si:ch211-131k2.3", "si:ch211-221f10.2", "si:dkey-184p18.2"),
  uncharacterized_regulators = c("dcaf6", "mxd1", "nfat5b"),
  uncharacterized_cytoskeleton_ecm = c("coro6", "fhod3a", "lmod2b")
)

for (set_name in names(gene_sets)) {
  genes <- gene_sets[[set_name]]

  plots <- VlnPlot(
    heart20s_final_reordered,
    features = genes,
    cols = c(
      "#0a9f0a", "darkgrey", "darkgrey", "darkgrey", "darkgrey",
      "darkgrey", "darkgrey", "darkgrey"
    )
    combine = FALSE
  )

  plots <- lapply(plots, function(p) {
    p + NoLegend() +
      theme(
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
        plot.title  = element_text(face = "bold.italic")
      ) +
      xlab("Cluster")
  })

  p_panel <- wrap_plots(plots)

  ggsave(
    filename = file.path(figures_dir, paste0("20ss_vln_panel_", set_name, ".png")),
    plot = p_panel,
    width = 12, height = 4, dpi = 300
  )
}

# Left/Right assignment within cardiomyocytes
heart20s_final_cluster3only <- subset(heart20s_final, idents = c("3"))

summary(FetchData(object = heart20s_final_cluster3only, vars = "pitx2"))
summary(FetchData(object = heart20s_final_cluster3only, vars = "lft2"))
summary(FetchData(object = heart20s_final_cluster3only, vars = "lft1"))
summary(FetchData(object = heart20s_final_cluster3only, vars = "ndr2"))
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

saveRDS(
  heart20s_final_cluster3only,
  file = file.path(results_dir, "heart20s_cluster3_LR.rds")
)

# LR differential expression
heart20s_LR_DEanalysis.markers <- FindMarkers(
  heart20s_final_cluster3only,
  ident.1 = "Lpositive",
  ident.2 = "Rnegative"
)

heart20s_LR_DEanalysis.markersFILTERED <- subset(
  heart20s_LR_DEanalysis.markers,
  p_val_adj < 0.05 & abs(avg_log2FC) > 0.25
)

write.csv(
  heart20s_LR_DEanalysis.markersFILTERED,
  file = file.path(results_dir, "LR_DEanalysis_heart20ss_filtered.csv")
)

# AUCell: pathway activity across clusters
exprMatrix <- GetAssayData(heart20s_final, slot = "counts")

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
heart20s_final <- AddMetaData(heart20s_final, auc_df)

# Reorder + rename clusters for plotting
heart20s_final_reordered <- SetIdent(
  heart20s_final,
  value = factor(Idents(heart20s_final), levels = c("3", "0", "1", "2", "4", "5", "6", "7"))
)
heart20s_final_reordered <- RenameIdents(heart20s_final_reordered, new_cluster_names_short)

p_auc <- VlnPlot(
  heart20s_final_reordered,
  features = names(geneSets),
  cols = c(
      "#0a9f0a", "darkgrey", "darkgrey", "darkgrey", "darkgrey",
      "darkgrey", "darkgrey", "darkgrey"
    )
) +
  theme(
    legend.position = "none",
    plot.title = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.title.y = element_text(size = 16)
  ) +
  labs(y = "Nodal signaling at 20 ss") +
  xlab("Cluster")

ggsave(
  filename = file.path(figures_dir, "20ss_AUCell_Nodal_by_cluster.png"),
  plot = p_auc,
  width = 7, height = 5, dpi = 300
)

# Save updated object w/ AUCell metadata
saveRDS(heart20s_final, file = file.path(results_dir, "heart20s_final_with_AUCell.rds"))

message("03_seurat_20ss.R complete. Outputs saved to results/20ss and figures/20ss.")
