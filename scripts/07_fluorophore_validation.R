# 07_fluorophore_validation.R
# Purpose: Validate that cardiac cells express transgenic fluorophore transcripts (ZsYellow / EGFP) via blended FeaturePlots overlaying fluorophore + canonical cardiac marker

# Data should be downloaded from GEO GSE311898 and placed in:
#   data/raw/16ss_fluorophore_files/
#   data/raw/18ss_fluorophore_files/
#   data/raw/20ss_fluorophore_files/

# Setup
source(here("scripts", "00_setup.R"))
set.seed(123)

# Output folders
results_dir <- here("results", "fluorophore_validation")
figures_dir <- here("figures", "fluorophore_validation")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------#
# 0) Load fluorophore Cell Ranger matrices
# -------------------------------#

data_dir <- here("data", "raw")

data_dir_16_fluoro <- file.path(data_dir, "16ss_fluorophore_files")
data_dir_18_fluoro <- file.path(data_dir, "18ss_fluorophore_files")
data_dir_20_fluoro <- file.path(data_dir, "20ss_fluorophore_files")

check_dir <- function(path) {
  if (!dir.exists(path)) {
    stop(
      paste(
        "\nMissing data folder:",
        path,
        "\nPlace Cell Ranger outputs in data/raw/*_fluorophore_files/ and re-run."
      ),
      call. = FALSE
    )
  }
}

check_dir(data_dir_16_fluoro)
check_dir(data_dir_18_fluoro)
check_dir(data_dir_20_fluoro)

heart16s_fluoro.data <- Read10X(data.dir = data_dir_16_fluoro)
heart18s_fluoro.data <- Read10X(data.dir = data_dir_18_fluoro)
heart20s_fluoro.data <- Read10X(data.dir = data_dir_20_fluoro)

message("Fluorophore datasets successfully loaded.")

# Seurat analysis, standard workflow, 16 ss fluorophore dataset 
heart16s <- CreateSeuratObject(counts = heart16s.data, project = "16somdata", min.cells = 3, min.features = 200)
heart16s[["percent.mt"]] <- PercentageFeatureSet(heart16s, pattern = "^mt-")
heart16s_filtered <- subset(heart16s, subset = nFeature_RNA > 200 & nFeature_RNA < 9096 & percent.mt < 10.46 & nCount_RNA < 100000)
heart16s_filtered_normalized <- NormalizeData(heart16s_filtered, normalization.method = "LogNormalize", scale.factor = 10000)
heart16s_filtered_normalized <- FindVariableFeatures(heart16s_filtered_normalized, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(heart16s_filtered_normalized), 10)
all.genes <- rownames(heart16s_filtered_normalized)
heart16s_filtered_normalized <- ScaleData(heart16s_filtered_normalized, features = all.genes)
heart16s_filtered_normalized <- RunPCA(heart16s_filtered_normalized, features = VariableFeatures(object = heart16s_filtered_normalized))
heart16s_filtered_normalized_PCA13 <- FindNeighbors(heart16s_filtered_normalized, dims = 1:13)
heart16s_filtered_normalized_PCA13_res0.12 <- FindClusters(heart16s_filtered_normalized_PCA13, resolution = .12)
heart16s_final <- RunUMAP(heart16s_filtered_normalized_PCA13_res0.12, dims = 1:13)
#Save final fluorophore dataset
saveRDS(heart16s_final, file = file.path(results_dir, "heart16s_fluorophore_final.rds"))

# Feature plots for ZsYellow and nkx2.5
blend_plot_16s <- FeaturePlot(
  heart16s_final,
  features = c("ZsYellow", "nkx2.5"),
  blend = TRUE,
  cols = c("#D41159", "#1A85FF")
)
blend_plot_16s[[1]] <- blend_plot_16s[[1]] + ggtitle("zsyellow") + theme(plot.title = element_text(face = "bold.italic"))
blend_plot_16s[[2]] <- blend_plot_16s[[2]] + ggtitle("nkx2.5") + theme(plot.title = element_text(face = "bold.italic"))
blend_plot_16s[[3]] <- blend_plot_16s[[3]] + ggtitle("zsyellow and nkx2.5") + theme(plot.title = element_text(face = "bold.italic"))
blend_plot_16s[[4]] <- blend_plot_16s[[4]] + labs(x = "zsyellow", y = "nkx2.5") + theme(axis.title.x = element_text(face = "italic"), axis.title.y = element_text(face = "italic"))
p_16 <- blend_plot_16s[[1]] + blend_plot_16s[[2]] + blend_plot_16s[[3]] + blend_plot_16s[[4]] + plot_layout(ncol = 4)
ggsave(
  filename = file.path(figures_dir, "16ss_fluorophore_blend_ZsYellow_nkx2p5.png"),
  plot = p_16,
  width = 14, height = 4, dpi = 300
)

# Seurat analysis, standard workflow, 18 ss fluorophore dataset 
heart18s <- CreateSeuratObject(counts = heart18s.data, project = "18somdata", min.cells = 3, min.features = 200)
heart18s[["percent.mt"]] <- PercentageFeatureSet(heart18s, pattern = "^mt-")
heart18s_filtered <- subset(heart18s, subset = nFeature_RNA > 200 & nFeature_RNA < 7974 & percent.mt < 7.79 & nCount_RNA < 100000)
heart18s_filtered_normalized <- NormalizeData(heart18s_filtered, normalization.method = "LogNormalize", scale.factor = 10000)
heart18s_filtered_normalized <- FindVariableFeatures(heart18s_filtered_normalized, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(heart18s_filtered_normalized), 10)
all.genes <- rownames(heart18s_filtered_normalized)
heart18s_filtered_normalized <- ScaleData(heart18s_filtered_normalized, features = all.genes)
heart18s_filtered_normalized <- RunPCA(heart18s_filtered_normalized, features = VariableFeatures(object = heart18s_filtered_normalized))
heart18s_filtered_normalized_PCA16 <- FindNeighbors(heart18s_filtered_normalized, dims = 1:16)
heart18s_filtered_normalized_PCA16_res0.15 <- FindClusters(heart18s_filtered_normalized_PCA16, resolution = .15)
heart18s_final <- RunUMAP(heart18s_filtered_normalized_PCA16_res0.15, dims = 1:16)
#Save final fluorophore dataset
saveRDS(heart18s_final, file = file.path(results_dir, "heart18s_fluorophore_final.rds"))

# Feature plots for ZsYellow and nkx2.5
blend_plot_18s <- FeaturePlot(
  heart18s_final,
  features = c("ZsYellow", "nkx2.5"),
  blend = TRUE,
  cols = c("#D41159", "#1A85FF")
)
blend_plot_18s[[1]] <- blend_plot_18s[[1]] + ggtitle("zsyellow") + theme(plot.title = element_text(face = "bold.italic"))
blend_plot_18s[[2]] <- blend_plot_18s[[2]] + ggtitle("nkx2.5") + theme(plot.title = element_text(face = "bold.italic"))
blend_plot_18s[[3]] <- blend_plot_18s[[3]] + ggtitle("zsyellow and nkx2.5") + theme(plot.title = element_text(face = "bold.italic"))
blend_plot_18s[[4]] <- blend_plot_18s[[4]] + labs(x = "zsyellow", y = "nkx2.5") + theme(axis.title.x = element_text(face = "italic"), axis.title.y = element_text(face = "italic"))
p_18 <- blend_plot_18s[[1]] + blend_plot_18s[[2]] + blend_plot_18s[[3]] + blend_plot_18s[[4]] + plot_layout(ncol = 4)
ggsave(
  filename = file.path(figures_dir, "18ss_fluorophore_blend_ZsYellow_nkx2p5.png"),
  plot = p_18,
  width = 14, height = 4, dpi = 300
)

# Seurat analysis, standard workflow, 20 ss fluorophore dataset 
heart20s <- CreateSeuratObject(counts = heart20s.data, project = "20somdata", min.cells = 3, min.features = 200)
heart20s[["percent.mt"]] <- PercentageFeatureSet(heart20s, pattern = "^mt-")
heart20s_filtered <- subset(heart20s, subset = nFeature_RNA > 200 & nFeature_RNA < 5820 & percent.mt < 5.22 & nCount_RNA < 100000)
heart20s_filtered_normalized <- NormalizeData(heart20s_filtered, normalization.method = "LogNormalize", scale.factor = 10000)
heart20s_filtered_normalized <- FindVariableFeatures(heart20s_filtered_normalized, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(heart20s_filtered_normalized), 10)
all.genes <- rownames(heart20s_filtered_normalized)
heart20s_filtered_normalized <- ScaleData(heart20s_filtered_normalized, features = all.genes)
heart20s_filtered_normalized <- RunPCA(heart20s_filtered_normalized, features = VariableFeatures(object = heart20s_filtered_normalized))
heart20s_filtered_normalized_PCA13 <- FindNeighbors(heart20s_filtered_normalized, dims = 1:13)
heart20s_filtered_normalized_PCA13_res0.2 <- FindClusters(heart20s_filtered_normalized_PCA13, resolution = .2)
heart20s_final <- RunUMAP(heart20s_filtered_normalized_PCA13_res0.2, dims = 1:13)
#Save final fluorophore dataset
saveRDS(heart20s_final, file = file.path(results_dir, "heart20s_fluorophore_final.rds"))

# Feature plots for EGFP and myl7
blend_plot_20s <- FeaturePlot(
  heart20s_final,
  features = c("EGFP", "myl7"),
  blend = TRUE,
  cols = c("#D41159", "#1A85FF")
)
blend_plot_20s[[1]] <- blend_plot_20s[[1]] + ggtitle("egfp") + theme(plot.title = element_text(face = "bold.italic"))
blend_plot_20s[[2]] <- blend_plot_20s[[2]] + ggtitle("myl7") + theme(plot.title = element_text(face = "bold.italic"))
blend_plot_20s[[3]] <- blend_plot_20s[[3]] + ggtitle("egfp and myl7") + theme(plot.title = element_text(face = "bold.italic"))
blend_plot_20s[[4]] <- blend_plot_20s[[4]] + labs(x = "egfp", y = "myl7") + theme(axis.title.x = element_text(face = "italic"), axis.title.y = element_text(face = "italic"))
p_20 <- blend_plot_20s[[1]] + blend_plot_20s[[2]] + blend_plot_20s[[3]] + blend_plot_20s[[4]] + plot_layout(ncol = 4)
ggsave(
  filename = file.path(figures_dir, "20ss_fluorophore_blend_EGFP_myl7.png"),
  plot = p_20,
  width = 14, height = 4, dpi = 300
)

message("07_fluorophore_validation.R complete. Outputs saved to results/fluorophore_validation and figures/fluorophore_validation.")
