# zebrafish-heart-scrnaseq
Code for Featured Publication: Single-cell transcriptomics reveals Nodal-responsive cytoskeletal programs driving asymmetric heart morphogenesis.

## Data availability
This repository does not include raw data due to file size limits.
Raw sequencing data and processed count matrices are available from: **NCBI GEO accession GSE311898**. 

To reproduce the analysis:

1. Download the processed barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz files for each dataset from GEO.

2. Extract and organize files as:

  data/raw/16ss_files/
  barcodes.tsv.gz
  features.tsv.gz
  matrix.mtx.gz

  data/raw/18ss_files/
  ...

  data/raw/20ss_files/
  ...

  data/raw/16ss__fluorophore_files/
  ...

  data/raw/18ss__fluorophore_files/
  ...

  data/raw/20ss__fluorophore_files/
  ...

3. Run analysis scripts in `scripts/`.
