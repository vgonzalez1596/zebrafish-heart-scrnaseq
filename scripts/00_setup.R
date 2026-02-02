# 00_setup.R
# Purpose: load packages + define paths + load data

# Data should be downloaded from GEO GSE311898 and placed in:
#   data/raw/16ss_files/
#   data/raw/18ss_files/
#   data/raw/20ss_files/

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(patchwork)
  library(ggplot2)
  library(ggtext)
  library(ggsci)
  library(eulerr)
  library(viridis)
  library(cowplot)
  library(tidyr)
  library(tidyverse)
  library(slingshot)
  library(tradeSeq)
  library(AUCell)
  library(GSEABase)
  library(SummarizedExperiment)
  library(here)
})

# Base data directory
data_dir <- here("data", "raw")

data_dir_16 <- file.path(data_dir, "16ss_files")
data_dir_18 <- file.path(data_dir, "18ss_files")
data_dir_20 <- file.path(data_dir, "20ss_files")

# Safety checks
check_dir <- function(path) {
  if (!dir.exists(path)) {
    stop(
      paste(
        "\nMissing data folder:",
        path,
        "\nDownload GEO data (GSE311898) and place it here."
      ),
      call. = FALSE
    )
  }
}

check_dir(data_dir_16)
check_dir(data_dir_18)
check_dir(data_dir_20)

# Load datasets
heart16s.data <- Read10X(data.dir = data_dir_16)
heart18s.data <- Read10X(data.dir = data_dir_18)
heart20s.data <- Read10X(data.dir = data_dir_20)

message("All datasets successfully loaded.")
