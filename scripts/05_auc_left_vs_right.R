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

# Load in the Seurat objects with left and right defined
heart16s_final_cluster3only <- readRDS(here("results", "16ss", "heart16s_cluster3_LR.rds"))
heart18s_final_cluster2only <- readRDS(here("results","18ss","heart18s_cluster2_LR.rds"))
heart20s_final_cluster3only <- readRDS(here("results","20ss","heart20s_cluster3_LR.rds"))

# AUCell analysis of left-exclusive canonical Nodal targets at 16 ss
# To instead analyze at 18 ss, replace all instances of heart16s_final_cluster3only with heart18s_final_cluster2only
# To instead analyze at 20 ss, replace all instances of heart16s_final_cluster3only with heart20s_final_cluster3only
Idents(heart16s_final_cluster3only) <- factor(Idents(heart16s_final_cluster3only), levels = c("Lpositive", "Rnegative", "3"))
table(Idents(heart16s_final_cluster3only))
exprMatrix <- GetAssayData(heart16s_final_cluster3only, slot = "counts")
geneSets <- list(Active_Nodal_LeftOnly = c("lft2", "pitx2", "ndr2", "lft1")) #this includes genes that are expressed on the left only
cells_rankings <- AUCell_buildRankings(exprMatrix, nCores = 1, plotStats = TRUE)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
AUCell_plotHist(cells_AUC)
auc_df <- as.data.frame(t(getAUC(cells_AUC)))
heart16s_final_cluster3only <- AddMetaData(heart16s_final_cluster3only, auc_df)
p_auc <- VlnPlot(heart16s_final_cluster3only, features = names(geneSets), idents = c("Lpositive", "Rnegative"), cols = c("#6998CA", "#F2A165")) + 
  theme(legend.position = "none", plot.title = element_blank(), axis.title.x = element_blank(), axis.text.x = element_text(angle = 0, hjust = 0.5)) +  
  labs(y = "Nodal signaling\n(Left-exclusive target genes)") +  
  scale_x_discrete(labels = c("Lpositive" = "Left", "Rnegative" = "Right")) + 
  xlab("16 ss")
ggsave(
  filename = file.path(figures_dir, "16ss_AUCell_Nodal_left_vs_right_LeftExclusive.png"),
  plot = p_auc,
  width = 7, height = 5, dpi = 300
)
auc_df <- as.data.frame(t(getAUC(cells_AUC)))
auc_df$cluster <- Idents(heart16s_final_cluster3only)
gene_set_name <- "Active_Nodal_LeftOnly"
p_value <- wilcox.test(
  auc_df[auc_df$cluster == "Rnegative", gene_set_name],
  auc_df[auc_df$cluster == "Lpositive", gene_set_name]
)$p.value
print(p_value)
write.csv(data.frame(stage = "16ss", gene_set = gene_set_name, p_value = p_value), file = file.path(results_dir, "16ss_AUCell_LeftOnly_wilcox_pvalue.csv"), row.names = FALSE)
saveRDS(heart16s_final_cluster3only, file = file.path(results_dir, "heart16s_cluster3_LR_with_AUCell_LeftOnly.rds"))

# AUCell analysis of all 7 canonical Nodal targets at 16 ss
# To instead analyze at 18 ss, replace all instances of heart16s_final_cluster3only with heart18s_final_cluster2only
# To instead analyze at 20 ss, replace all instances of heart16s_final_cluster3only with heart20s_final_cluster3only
Idents(heart16s_final_cluster3only) <- factor(Idents(heart16s_final_cluster3only), levels = c("Lpositive", "Rnegative", "3"))
table(Idents(heart16s_final_cluster3only))
exprMatrix <- GetAssayData(heart16s_final_cluster3only, slot = "counts")
geneSets <- list(Active_Nodal = c("lft2", "pitx2", "ndr2", "lft1", "acta1b", "has2", "bmp4")) #this includes all canonical genes
cells_rankings <- AUCell_buildRankings(exprMatrix, nCores = 1, plotStats = TRUE)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
AUCell_plotHist(cells_AUC)
auc_df <- as.data.frame(t(getAUC(cells_AUC)))
heart16s_final_cluster3only <- AddMetaData(heart16s_final_cluster3only, auc_df)
p_auc <- VlnPlot(heart16s_final_cluster3only, features = names(geneSets), idents = c("Lpositive", "Rnegative"), cols = c("#6998CA", "#F2A165")) + 
  theme(legend.position = "none", plot.title = element_blank(), axis.title.x = element_blank(), axis.text.x = element_text(angle = 0, hjust = 0.5)) +  
  labs(y = "Nodal signaling\n(All known target genes)") +  
  scale_x_discrete(labels = c("Lpositive" = "Left", "Rnegative" = "Right")) + 
  xlab("16 ss")
ggsave(
  filename = file.path(figures_dir, "16ss_AUCell_Nodal_left_vs_right_AllCanonical.png"),
  plot = p_auc,
  width = 7, height = 5, dpi = 300
)
auc_df <- as.data.frame(t(getAUC(cells_AUC)))
auc_df$cluster <- Idents(heart16s_final_cluster3only)
gene_set_name <- "Active_Nodal"
p_value <- wilcox.test(
  auc_df[auc_df$cluster == "Rnegative", gene_set_name],
  auc_df[auc_df$cluster == "Lpositive", gene_set_name]
)$p.value
print(p_value)
write.csv(data.frame(stage = "16ss", gene_set = gene_set_name, p_value = p_value), file = file.path(results_dir, "16ss_AUCell_AllCanonical_wilcox_pvalue.csv"), row.names = FALSE)
saveRDS(heart16s_final_cluster3only, file = file.path(results_dir, "heart16s_cluster3_LR_with_AUCell_AllCanonical.rds"))

# As a control, we carry out AUCell analysis at 20 ss with a set of confirmed cardiac cone markers (ZFIN) that are not involved in LR patterning or Nodal signaling
Idents(heart20s_final_cluster3only) <- factor(Idents(heart20s_final_cluster3only), levels = c("Lpositive", "Rnegative", "3"))
table(Idents(heart20s_final_cluster3only))
exprMatrix <- GetAssayData(heart20s_final_cluster3only, slot = "counts")
geneSets <- list(Cardiac_Cone_Genes = c("myh6", "myh7", "ttn.1", "tbx5a", "mef2ca", "gata5", "tbx20"))
cells_rankings <- AUCell_buildRankings(exprMatrix, nCores = 1, plotStats = TRUE)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
AUCell_plotHist(cells_AUC)
auc_df <- as.data.frame(t(getAUC(cells_AUC)))
heart20s_final_cluster3only <- AddMetaData(heart20s_final_cluster3only, auc_df)
p_auc <- VlnPlot(heart20s_final_cluster3only, features = names(geneSets), idents = c("Lpositive", "Rnegative"), cols = c("#6998CA", "#F2A165")) + 
  theme(legend.position = "none", plot.title = element_blank(), axis.title.x = element_blank(), axis.text.x = element_text(angle = 0, hjust = 0.5)) +  
  labs(y = "Expression of cardiac cone markers\nunrelated to Nodal signaling") +  
  scale_x_discrete(labels = c("Lpositive" = "Left", "Rnegative" = "Right")) + 
  xlab("20 ss")
ggsave(
  filename = file.path(figures_dir, "20ss_AUCell_Nodal_left_vs_right_ConeMarkers.png"),
  plot = p_auc,
  width = 7, height = 5, dpi = 300
)
auc_df <- as.data.frame(t(getAUC(cells_AUC)))
auc_df$cluster <- Idents(heart20s_final_cluster3only)
gene_set_name <- "Cardiac_Cone_Genes"
p_value <- wilcox.test(
  auc_df[auc_df$cluster == "Rnegative", gene_set_name],
  auc_df[auc_df$cluster == "Lpositive", gene_set_name]
)$p.value
print(p_value)
write.csv(data.frame(stage = "20ss", gene_set = gene_set_name, p_value = p_value), file = file.path(results_dir, "20ss_AUCell_ConeMarkers_wilcox_pvalue.csv"), row.names = FALSE)
saveRDS(heart20s_final_cluster3only, file = file.path(results_dir, "heart20s_cluster3_LR_with_AUCell_AllCanonical.rds"))

message("05_auc_left_vs_right.R complete. Outputs saved to results/auc_left_vs_right and figures/auc_left_vs_right.")
