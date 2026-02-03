# 04_lr_deg_overlap_venn.R
# Purpose: Visualize overlap between Left/Right DEGs across 16 ss, 18 ss, and 20 ss via Euler/Venn-style diagram

# Setup
source(here("scripts", "00_setup.R"))
set.seed(123)

# Output folders
results_dir <- here("results", "lr_deg_overlap")
figures_dir <- here("figures", "lr_deg_overlap")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# LR DEG gene lists
genes_16s <- c(
  "lft2","elovl6","pitx2","lft1","ndr2","tmem88a","arhgap42b",
  "abca1b","pald1a","bmpr2b","jak2a","skia","isl2a","swap70b","ldlrad4a","spaw"
)

genes_18s <- c(
  "lft2","elovl6","lft1","pitx2","ndr2","ephb2b","jak2a",
  "tgfb1a","abca1b","pald1a","arhgap42b","lmna","trib3",
  "ldlrad4a","slc12a2","itgb1b.2","stk17al","swap70b"
)

genes_20s <- c(
  "lft2","mibp2","elovl6","acta1a","actc2","crip2l","apln",
  "csrp1a","nkx2.5","tnni2a.4","apobec2a","dbn1","pitx2",
  "ptmab","hmgb3a","slc12a2","actc1a","tmem88a","ephb2b",
  "tgfb1a","hnrnpa0a","ubc","bckdk","sult2st1","abca1b",
  "ndr2","gng13b","popdc1","prkar1aa","lft1","pold4",
  "p2rx3b","ryr3","jak2a","pald1a","tanc1b","raraa",
  "mef2cb","gse1","diaph2","cntfr","hdac4","cacna2d2a","aff2"
)

# Fit Venn diagram and plot
sets <- list(`16 ss` = genes_16s, `18 ss` = genes_18s, `20 ss` = genes_20s)
fit <- euler(sets)
plot_euler <- function() {
  plot(
    fit,
    fills = c("#E56997", "#FBC740", "#66D2D6"),
    edges = list(col = "black", lwd = 1.5),
    labels = list(font = 2, cex = 1.2),
    quantities = list(cex = 1.1),
    legend = FALSE
  )
}

# Save Venn diagram as PNG
png(file.path(figures_dir, "eulerr_venn_diagram.png"), width = 1800, height = 1800, res = 300)
plot_euler()
dev.off()

message("04_lr_deg_overlap_venn.R complete. Outputs saved to results/lr_deg_overlap and figures/lr_deg_overlap.")
