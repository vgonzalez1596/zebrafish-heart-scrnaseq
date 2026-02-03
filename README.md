# Featured Publication: Single-cell transcriptomics reveals Nodal-responsive cytoskeletal programs driving asymmetric heart morphogenesis.

This repository contains the R code used to analyse the single-cell RNA-seq datasets from the above publication. 
*Gonzalez et al. Nature Communications (under revision)*

## Project Goal
To characterize cardiomyocyte transcriptional programs and signaling-dependent gene expression changes across developmental timepoints in the zebrafish heart using single-cell RNA sequencing.

## Analysis Workflow
This repository contains scripts used to:

• Perform scRNA-seq quality control and filtering  
• Normalize and integrate datasets across developmental timepoints  
• Cluster and annotate cardiomyocyte populations  
• Perform differential gene expression analysis  
• Score signaling pathway and gene set activity  
• Perform pseudotime and trajectory inference  
• Generate publication-quality visualizations  

## Example Output

![UMAP](figures/umap.png)

---

## Tools & Packages
Seurat, Slingshot, AUCell, tidyverse, ggplot2, and related R-based single-cell analysis tools.

---

## Repository Structure



## Data availability
This repository does not include raw data due to file size limits.

Raw sequencing data and processed count matrices are available from: **NCBI GEO accession GSE311898**. 

To reproduce the analysis:

1. Download the processed barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz files for each dataset from GEO.

2. Extract and organize files in their corresponding location in `data/raw/`.

3. Run analysis scripts in `scripts/`.
