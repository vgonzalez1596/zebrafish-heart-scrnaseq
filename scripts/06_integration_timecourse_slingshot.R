# 06_integration_timecourse_slingshot.R
# Purpose:
#   - Integrate the 16 ss, 18 ss, and 20 ss datasets
#   - Visualize timecourse trends in gene expression
#   - Run Slingshot pseudotime on integrated cardiomyocytes
#   - AUCell pathway scoring across time in cardiomyocytes

# Setup
source(here("scripts", "00_setup.R"))
set.seed(123)

# Output folders
results_dir <- here("results", "integration")
figures_dir <- here("figures", "integration")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# Seurat integration and analysis, standard workflow
heart16s <- CreateSeuratObject(counts = heart16s.data, project = "16somdata", min.cells = 3, min.features = 200)
heart16s[["percent.mt"]] <- PercentageFeatureSet(heart16s, pattern = "^mt-")
heart16s_filtered <- subset(heart16s, subset = nFeature_RNA > 200 & nFeature_RNA < 9096 & percent.mt < 10.46 & nCount_RNA < 100000)
heart18s <- CreateSeuratObject(counts = heart18s.data, project = "18somdata", min.cells = 3, min.features = 200)
heart18s[["percent.mt"]] <- PercentageFeatureSet(heart18s, pattern = "^mt-")
heart18s_filtered <- subset(heart18s, subset = nFeature_RNA > 200 & nFeature_RNA < 7974 & percent.mt < 7.79 & nCount_RNA < 100000)
heart20s <- CreateSeuratObject(counts = heart20s.data, project = "20somdata", min.cells = 3, min.features = 200)
heart20s[["percent.mt"]] <- PercentageFeatureSet(heart20s, pattern = "^mt-")
heart20s_filtered <- subset(heart20s, subset = nFeature_RNA > 200 & nFeature_RNA < 5820 & percent.mt < 5.22 & nCount_RNA < 100000)
#Save the filtered objects
saveRDS(heart16s_filtered, file = file.path(results_dir, "heart16s_filtered_for_integration.rds"))
saveRDS(heart18s_filtered, file = file.path(results_dir, "heart18s_filtered_for_integration.rds"))
saveRDS(heart20s_filtered, file = file.path(results_dir, "heart20s_filtered_for_integration.rds"))

# Merge and integrate the datasets
heartmerged <- merge(heart16s_filtered, y = c(heart18s_filtered, heart20s_filtered), add.cell.ids = c("16som", "18som", "20som"), project = "heartmerged_16_18_20som")
unique(sapply(X = strsplit(colnames(heartmerged), split = "_"), FUN = "[", 1))
heartmerged.list <- SplitObject(heartmerged, split.by = "orig.ident")
heartmerged.list <- lapply(X = heartmerged.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
integration_features <- SelectIntegrationFeatures(object.list = heartmerged.list)
heart.anchors <- FindIntegrationAnchors(object.list = heartmerged.list, anchor.features = integration_features)
heart_integrated <- IntegrateData(anchorset = heart.anchors)
#Save raw integrated object
saveRDS(heart_integrated, file = file.path(results_dir, "heart_integrated_raw.rds"))

# Processing integrated dataset
DefaultAssay(heart_integrated) <- "integrated"
heart_integrated <- ScaleData(heart_integrated)
heart_integrated <- RunPCA(heart_integrated, npcs = 30)
heart_integrated <- RunUMAP(heart_integrated, reduction = "pca", dims = 1:25)
heart_integrated <- FindNeighbors(heart_integrated, reduction = "pca", dims = 1:25)
heart_integrated <- FindClusters(heart_integrated, resolution = 0.25)
# Save processed integrated object
saveRDS(heart_integrated, file = file.path(results_dir, "heart_integrated_clustered.rds"))
p_umap_by_time <- DimPlot(heart_integrated, reduction = "umap", group.by = "orig.ident")
ggsave(
  filename = file.path(figures_dir, "integration_umap_by_timepoint.png"),
  plot = p_umap_by_time,
  width = 7, height = 5, dpi = 300
)
DefaultAssay(heart_integrated) <- "RNA"
heart_integrated_heartonly <- subset(heart_integrated, idents = c("2"))
# Save integrated dataset containing only cardiomyocytes
saveRDS(heart_integrated_heartonly, file = file.path(results_dir, "heart_integrated_heartonly_subset.rds"))

# Defining left and right cells
summary(FetchData(object = heart_integrated_heartonly, vars = 'pitx2'))
summary(FetchData(object = heart_integrated_heartonly, vars = 'lft2'))
summary(FetchData(object = heart_integrated_heartonly, vars = 'lft1'))
summary(FetchData(object = heart_integrated_heartonly, vars = 'ndr2'))
# Mean expression: lft1 = 0.1176, lft2 = 0.9449, ndr2 = 0.09055, pitx2 = 0.184
left_sided_genes_positive <- WhichCells(heart_integrated_heartonly, expression = pitx2 > .184 | lft2 > .9449 | lft1 > .1176 | ndr2 > .09055)
heart_integrated_heartonly <- SetIdent(heart_integrated_heartonly, cells = left_sided_genes_positive, value = 'Lpositive')
left_sided_genes_negative <- WhichCells(heart_integrated_heartonly, expression = lft2 == 0 & pitx2 == 0 & lft1 == 0 & ndr2 == 0)
heart_integrated_heartonly <- SetIdent(heart_integrated_heartonly, cells = left_sided_genes_negative, value = 'Rnegative')
# Save integrated dataset containing only cardiomyocytes with left and right identities
saveRDS(heart_integrated_heartonly, file = file.path(results_dir, "heart_integrated_heartonly_with_LR_idents.rds"))

# Visualizing gene expression in the cardiac cells of the integrated dataset via scatter plot with a LOESS trendline
gene_of_interest <- "lft2" #Replace with your gene of interest
color_ident1 <- "#3F78C1"
color_ident2 <- "#D37538"
expression_data <- FetchData(heart_integrated_heartonly, vars = c(gene_of_interest, "orig.ident"))
expression_data$Ident <- Idents(heart_integrated_heartonly)
idents_to_keep <- c("Lpositive", "Rnegative")
expression_data <- expression_data %>% filter(Ident %in% idents_to_keep)
expression_data$orig.ident <- factor(expression_data$orig.ident, levels = c("16somdata", "18somdata", "20somdata"))
p_timecourse <- ggplot(expression_data, aes(x = orig.ident, y = !!rlang::sym(gene_of_interest), color = Ident, group = Ident)) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 1.2) +
  scale_color_manual(values = c("Lpositive" = color_ident1, "Rnegative" = color_ident2)) +
  labs(
    title = paste(gene_of_interest),
    x = "Time point",
    y = "Expression Level"
  ) +
  scale_x_discrete(labels = c("16somdata" = "16 ss", "18somdata" = "18 ss", "20somdata" = "20 ss")) +
  ylim(c(0, 5)) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.line = element_line(color = "black"),
    plot.title = element_text(face = "bold.italic", size = 16, hjust = 0.5),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14)
  )
ggsave(
  filename = file.path(figures_dir, paste0("integration_timecourse_loess_", gene_of_interest, ".png")),
  plot = p_timecourse,
  width = 8, height = 5, dpi = 300
)

# Visualizing epithelial and mesenchymal genes in the cardiac cells of the integrated dataset with feature plots and blended feature plots
epigenes <- c("cdh1", "epcam", "krt8", "krt18a.1")
epiplot_list <- lapply(epigenes, function(gene) {
  FeaturePlot(heart_integrated_heartonly, features = gene, cols = c("lightgray", "#005AB5")) +
    theme(plot.title = element_text(face = "bold.italic"))
})
p_epi <- wrap_plots(epiplot_list, ncol = 2, nrow = 2)
ggsave(
  filename = file.path(figures_dir, "integration_heartonly_epithelial_featureplots.png"),
  plot = p_epi,
  width = 10, height = 8, dpi = 300
)
mesgenes <- c("cdh2", "snai1b", "twist1b", "zeb1b")
mesplot_list <- lapply(mesgenes, function(gene) {
  FeaturePlot(heart_integrated_heartonly, features = gene, cols = c("lightgray", "#DC3220")) +
    theme(plot.title = element_text(face = "bold.italic"))
})
p_mes <- wrap_plots(mesplot_list, ncol = 2, nrow = 2)
ggsave(
  filename = file.path(figures_dir, "integration_heartonly_mesenchymal_featureplots.png"),
  plot = p_mes,
  width = 10, height = 8, dpi = 300
)
blend_plot_EMT <- FeaturePlot(
  heart_integrated_heartonly,
  features = c("krt8", "cdh2"),
  blend = TRUE,
  cols = c("#005AB5", "#DC3220")
)
blend_plot_EMT[[1]] <- blend_plot_EMT[[1]] + ggtitle("krt8") + theme(plot.title = element_text(face = "bold.italic"))
blend_plot_EMT[[2]] <- blend_plot_EMT[[2]] + ggtitle("cdh2") + theme(plot.title = element_text(face = "bold.italic"))
blend_plot_EMT[[3]] <- blend_plot_EMT[[3]] + ggtitle("krt8 and cdh2") + theme(plot.title = element_text(face = "bold.italic"))
blend_plot_EMT[[4]] <- blend_plot_EMT[[4]] + labs(x = "krt8", y = "cdh2") + theme(axis.title.x = element_text(face = "italic"), axis.title.y = element_text(face = "italic"))
p_blend <- blend_plot_EMT[[1]] + blend_plot_EMT[[2]] + blend_plot_EMT[[3]] + blend_plot_EMT[[4]] + plot_layout(ncol = 4)
ggsave(
  filename = file.path(figures_dir, "integration_heartonly_blended_krt8_cdh2.png"),
  plot = p_blend,
  width = 16, height = 4.5, dpi = 300
)
epi_mes_pairs <- list(c("krt8", "cdh2"), c("krt18a.1", "zeb1b"), c("cdh1", "snai1b"))
threshold <- 0.1
featurescatter_epi_mes <- lapply(epi_mes_pairs, function(pair) {
  df <- FetchData(heart_integrated_heartonly, vars = pair)
  heart_integrated_heartonly$DoublePositive_tmp <- ifelse(df[, 1] > threshold & df[, 2] > threshold, "Double+", "Other")
  percent_dp <- mean(heart_integrated_heartonly$DoublePositive_tmp == "Double+") * 100
  percent_text <- paste0("Double+ = ", round(percent_dp, 1), "%")
  FeatureScatter(
    heart_integrated_heartonly,
    feature1 = pair[1],
    feature2 = pair[2],
    group.by = "DoublePositive_tmp",
    pt.size = 0.4
  ) +
    scale_color_manual(values = c("Other" = "lightgrey", "Double+" = "#dc8cd4")) +
    labs(title = paste(pair[1], "vs", pair[2]), x = pair[1], y = pair[2]) +
    theme(
      plot.title = element_text(face = "bold.italic"),
      legend.position = "none",
      axis.title.x = element_text(face = "italic"),
      axis.title.y = element_text(face = "italic")
    ) +
    annotate(
      "text",
      x = max(df[, 1], na.rm = TRUE) * 0.7,
      y = max(df[, 2], na.rm = TRUE) * 0.9,
      label = percent_text,
      size = 4,
      fontface = "bold",
      color = "#dc8cd4"
    )
})
p_scatter <- wrap_plots(featurescatter_epi_mes, ncol = 3)
ggsave(
  filename = file.path(figures_dir, "integration_heartonly_featurescatter_epi_mes.png"),
  plot = p_scatter,
  width = 14, height = 5, dpi = 300
)

# Carrying out slingshot pseudotime analysis on the cardiac cells of the integrated dataset
heart_integrated_heartonly <- NormalizeData(heart_integrated_heartonly)
heart_integrated_heartonly <- FindVariableFeatures(heart_integrated_heartonly)
heart_integrated_heartonly <- ScaleData(heart_integrated_heartonly)
heart_integrated_heartonly <- RunPCA(heart_integrated_heartonly)
heart_integrated_heartonly <- FindNeighbors(heart_integrated_heartonly, dims = 1:30, reduction = "pca")
heart_integrated_heartonly <- FindClusters(heart_integrated_heartonly, cluster.name = "unintegrated_clusters")
heart_integrated_heartonly <- RunUMAP(heart_integrated_heartonly, dims = 1:30, reduction = "pca")
p_heart_umap_by_time <- DimPlot(heart_integrated_heartonly, group.by = "orig.ident")
ggsave(
  filename = file.path(figures_dir, "slingshot_heartonly_umap_by_timepoint.png"),
  plot = p_heart_umap_by_time,
  width = 7, height = 5, dpi = 300
)
heart_integrated_heartonly[["RNA"]] <- as(object = heart_integrated_heartonly[["RNA"]], Class = "Assay")
heart_integrated_heartonly$seurat_clusters <- heart_integrated_heartonly$orig.ident
sce <- as.SingleCellExperiment(heart_integrated_heartonly, assay = "RNA")
sce_slingshot <- slingshot(
  sce,
  reducedDim = "PCA",
  clusterLabels = colData(sce)$orig.ident,
  start.clus = "16somdata"
)

# Pseudotime plot with original dataset identities
pca_mat <- reducedDim(sce_slingshot, "PCA")
df_pca <- data.frame(
  PC1 = pca_mat[, 1],
  PC2 = pca_mat[, 2],
  orig.ident = colData(sce_slingshot)$orig.ident
)
p_pca_timepoint <- ggplot(df_pca, aes(x = PC1, y = PC2, color = orig.ident)) +
  geom_point(size = 0.6, alpha = 0.8) +
  scale_color_manual(
    values = c("16somdata" = "#f8766d",
               "18somdata" = "#00ba38",
               "20somdata" = "#619cff"),
    labels = c("16somdata" = "16 ss",
               "18somdata" = "18 ss",
               "20somdata" = "20 ss")
  ) +
  labs(x = "PC_1", y = "PC_2", color = "Timepoint") +
  theme_classic() +
  theme(
    legend.position = c(0.98, 0.98),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key = element_rect(fill = "white")
  )
ggsave(
  filename = file.path(figures_dir, "slingshot_pca_timepoint_colored.png"),
  plot = sce_slingshot_plot,
  width = 7, height = 5, dpi = 300
)

# Pseudotime plot with pseudotime coloring
pseudotime_values <- slingPseudotime(sce_slingshot, na = FALSE)
pt <- pseudotime_values[, 1]
pt_col <- viridis(100)[as.numeric(cut(pt, breaks = 100))]
png(
  filename = file.path(out_dir, "slingshot_pseudotime_colored.png"),
  width = 1800, height = 1400, res = 250
)
par(mar = c(5, 4, 4, 6), bty = "n")
plot(
  pca_mat,
  col = pt_col,
  pch = 16,
  cex = 0.5,
  xlab = "PC_1",
  ylab = "PC_2",
  main = NULL
)
lines(SlingshotDataSet(sce_slingshot), lwd = 2, col = "black")
breaks <- seq(min(pt, na.rm = TRUE), max(pt, na.rm = TRUE), length.out = 6)
legend(
  "topright",
  legend = round(breaks, 2),
  fill = viridis(6),
  title = "Pseudotime",
  border = "black",
  xpd = TRUE,
  inset = c(-0.2, 0),
  cex = 0.8
)
dev.off()

# Pseudotime FeaturePlots
features <- c("nkx2.5", "myl7", "acta1b", "bmp4", "has2", "lft1", "lft2", "ndr2", "pitx2")
plots <- list()
for (feature in features) {
  p <- FeaturePlot(heart_integrated_heartonly, features = feature, reduction = "pca") +
    ggtitle(feature) +
    theme(plot.title = element_text(face = "bold.italic", size = 14))
  plots[[feature]] <- p
}
p_pt_features <- plot_grid(plotlist = plots, align = "hv", ncol = 5)
ggsave(
  filename = file.path(figures_dir, "slingshot_pseudotime_featureplots_pca.png"),
  plot = p_pt_features,
  width = 16, height = 9, dpi = 300
)

# Visualizing Nodal signaling in the heart clusters of the integrated dataset via AUCell
DefaultAssay(heart_integrated_heartonly) <- "RNA"
exprMatrix <- GetAssayData(heart_integrated_heartonly, slot = "data")

# 1) Nodal
nodal_genes <- c("lft1", "lft2", "ndr2", "pitx2")
geneSets <- list("Nodal_Signaling" = nodal_genes)
cells_rankings <- AUCell_buildRankings(exprMatrix, nCores = 1, plotStats = TRUE)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
auc_scores <- as.data.frame(t(getAUC(cells_AUC)))
heart_integrated_heartonly$nodal_activity <- auc_scores$Nodal_Signaling
p_nodal <- VlnPlot(
  heart_integrated_heartonly,
  features = "nodal_activity",
  group.by = "orig.ident",
  pt.size = 0.1,
  cols = c("darkgrey", "darkgrey", "darkgrey")
) + ylab("Nodal signaling response") +
  xlab("Time point") +
  scale_x_discrete(labels = c("16somdata" = "16 ss", "18somdata" = "18 ss", "20somdata" = "20 ss")) +
  theme(
    legend.position = "none",
    plot.title = element_blank(),
    axis.text.x = element_text(angle = 0, vjust = 1, hjust = 0.5)
  ) +
  coord_cartesian(ylim = c(0, 0.8))
ggsave(
  filename = file.path(figures_dir, "integration_AUCell_Nodal_by_timepoint.png"),
  plot = p_nodal,
  width = 7, height = 5, dpi = 300
)

# 2) BMP
bmp_genes <- c("nkx2.5", "hand2", "tbx20", "tbx2b")
geneSets <- list("BMP_Signaling" = bmp_genes)
cells_rankings <- AUCell_buildRankings(exprMatrix, nCores = 1, plotStats = TRUE)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
auc_scores <- as.data.frame(t(getAUC(cells_AUC)))
heart_integrated_heartonly$bmp_activity <- auc_scores$BMP_Signaling
p_bmp <- VlnPlot(
  heart_integrated_heartonly,
  features = "bmp_activity",
  group.by = "orig.ident",
  pt.size = 0.1,
  cols = c("darkgrey", "darkgrey", "darkgrey")
) + ylab("BMP signaling response") +
  xlab("Time point") +
  scale_x_discrete(labels = c("16somdata" = "16 ss", "18somdata" = "18 ss", "20somdata" = "20 ss")) +
  theme(
    legend.position = "none",
    plot.title = element_blank(),
    axis.text.x = element_text(angle = 0, vjust = 1, hjust = 0.5)
  ) +
  coord_cartesian(ylim = c(0, 0.8))
ggsave(
  filename = file.path(figures_dir, "integration_AUCell_BMP_by_timepoint.png"),
  plot = p_bmp,
  width = 7, height = 5, dpi = 300
)

# 3) FGF
fgf_genes <- c("pea3", "erm", "gata4", "spry4")
geneSets <- list("FGF_Signaling" = fgf_genes)
cells_rankings <- AUCell_buildRankings(exprMatrix, nCores = 1, plotStats = TRUE)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
auc_scores <- as.data.frame(t(getAUC(cells_AUC)))
heart_integrated_heartonly$fgf_activity <- auc_scores$FGF_Signaling
p_fgf <- VlnPlot(
  heart_integrated_heartonly,
  features = "fgf_activity",
  group.by = "orig.ident",
  pt.size = 0.1,
  cols = c("darkgrey", "darkgrey", "darkgrey")
) + ylab("FGF signaling response") +
  xlab("Time point") +
  scale_x_discrete(labels = c("16somdata" = "16 ss", "18somdata" = "18 ss", "20somdata" = "20 ss")) +
  theme(
    legend.position = "none",
    plot.title = element_blank(),
    axis.text.x = element_text(angle = 0, vjust = 1, hjust = 0.5)
  ) +
  coord_cartesian(ylim = c(0, 0.8))
ggsave(
  filename = file.path(figures_dir, "integration_AUCell_FGF_by_timepoint.png"),
  plot = p_fgf,
  width = 7, height = 5, dpi = 300
)

# 4) Hedgehog
hedgehog_genes <- c("ptch1", "ptch2", "gli1", "gli2a", "gli2b", "gli3")
geneSets <- list("Hedgehog_Signaling" = hedgehog_genes)
cells_rankings <- AUCell_buildRankings(exprMatrix, nCores = 1, plotStats = TRUE)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
auc_scores <- as.data.frame(t(getAUC(cells_AUC)))
heart_integrated_heartonly$hedgehog_activity <- auc_scores$Hedgehog_Signaling
p_hh <- VlnPlot(
  heart_integrated_heartonly,
  features = "hedgehog_activity",
  group.by = "orig.ident",
  pt.size = 0.1,
  cols = c("darkgrey", "darkgrey", "darkgrey")
) + ylab("Hedgehog signaling response") +
  xlab("Time point") +
  scale_x_discrete(labels = c("16somdata" = "16 ss", "18somdata" = "18 ss", "20somdata" = "20 ss")) +
  theme(
    legend.position = "none",
    plot.title = element_blank(),
    axis.text.x = element_text(angle = 0, vjust = 1, hjust = 0.5)
  ) +
  coord_cartesian(ylim = c(0, 0.8))
ggsave(
  filename = file.path(figures_dir, "integration_AUCell_Hedgehog_by_timepoint.png"),
  plot = p_hh,
  width = 7, height = 5, dpi = 300
)

# 5) Notch
notch_genes <- c("wnt9a", "wif1", "notum1b", "msx1b", "msx3")
geneSets <- list("Notch_Signaling" = notch_genes)
cells_rankings <- AUCell_buildRankings(exprMatrix, nCores = 1, plotStats = TRUE)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
auc_scores <- as.data.frame(t(getAUC(cells_AUC)))
heart_integrated_heartonly$notch_activity <- auc_scores$Notch_Signaling
p_notch <- VlnPlot(
  heart_integrated_heartonly,
  features = "notch_activity",
  group.by = "orig.ident",
  pt.size = 0.1,
  cols = c("darkgrey", "darkgrey", "darkgrey")
) + ylab("Notch signaling response") +
  xlab("Time point") +
  scale_x_discrete(labels = c("16somdata" = "16 ss", "18somdata" = "18 ss", "20somdata" = "20 ss")) +
  theme(
    legend.position = "none",
    plot.title = element_blank(),
    axis.text.x = element_text(angle = 0, vjust = 1, hjust = 0.5)
  ) +
  coord_cartesian(ylim = c(0, 0.8))
ggsave(
  filename = file.path(figures_dir, "integration_AUCell_Notch_by_timepoint.png"),
  plot = p_notch,
  width = 7, height = 5, dpi = 300
)

# 6) Retinoic acid
ra_genes <- c("hoxb5b", "hoxb5a", "cyp26a1")
geneSets <- list("Retinoic_Acid_Signaling" = ra_genes)
cells_rankings <- AUCell_buildRankings(exprMatrix, nCores = 1, plotStats = TRUE)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
auc_scores <- as.data.frame(t(getAUC(cells_AUC)))
heart_integrated_heartonly$ra_activity <- auc_scores$Retinoic_Acid_Signaling
p_ra <- VlnPlot(
  heart_integrated_heartonly,
  features = "ra_activity",
  group.by = "orig.ident",
  pt.size = 0.1,
  cols = c("darkgrey", "darkgrey", "darkgrey")
) + ylab("Retinoic acid signaling response") +
  xlab("Time point") +
  scale_x_discrete(labels = c("16somdata" = "16 ss", "18somdata" = "18 ss", "20somdata" = "20 ss")) +
  theme(
    legend.position = "none",
    plot.title = element_blank(),
    axis.text.x = element_text(angle = 0, vjust = 1, hjust = 0.5)
  ) +
  coord_cartesian(ylim = c(0, 0.8))
ggsave(
  filename = file.path(figures_dir, "integration_AUCell_RA_by_timepoint.png"),
  plot = p_ra,
  width = 7, height = 5, dpi = 300
)

# 7) Wnt
wnt_genes <- c("ccnd1", "axin2", "myca")
geneSets <- list("Wnt_Signaling" = wnt_genes)
cells_rankings <- AUCell_buildRankings(exprMatrix, nCores = 1, plotStats = TRUE)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
auc_scores <- as.data.frame(t(getAUC(cells_AUC)))
heart_integrated_heartonly$wnt_activity <- auc_scores$Wnt_Signaling
p_wnt <- VlnPlot(
  heart_integrated_heartonly,
  features = "wnt_activity",
  group.by = "orig.ident",
  pt.size = 0.1,
  cols = c("darkgrey", "darkgrey", "darkgrey")
) + ylab("Wnt signaling response") +
  xlab("Time point") +
  scale_x_discrete(labels = c("16somdata" = "16 ss", "18somdata" = "18 ss", "20somdata" = "20 ss")) +
  theme(
    legend.position = "none",
    plot.title = element_blank(),
    axis.text.x = element_text(angle = 0, vjust = 1, hjust = 0.5)
  ) +
  coord_cartesian(ylim = c(0, 0.8))
ggsave(
  filename = file.path(figures_dir, "integration_AUCell_Wnt_by_timepoint.png"),
  plot = p_wnt,
  width = 7, height = 5, dpi = 300
)

# Save final heart-only object with AUCell metadata
saveRDS(heart_integrated_heartonly, file = file.path(results_dir, "heart_integrated_heartonly_with_LR_slingshot_AUCell.rds"))

message("06_integration_timecourse_slingshot.R complete. Outputs saved to results/integration and figures/integration.")
