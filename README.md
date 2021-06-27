## RISC


### Overview
Integrated analysis of single cell RNA-sequencing (scRNA-seq) data from multiple batches or studies is often necessary in order to learn functional changes in cellular states upon experimental perturbations or cell type relationships in a developmental lineage. Here we introduce a new algorithm (RPCI) that uses the gene-eigenvectors from a reference dataset to establish a global frame for projecting all datasets, with a clear advantage in preserving genuine gene expression differences in matching cell types between samples, such as those present in cells at distinct developmental stages or in perturbated vs control studies. This R package “RISC” (Robust Integration of Sinlgle Cell RNA-seq data implements the RPCI algorithm, with additional functions for scRNA-seq analysis, such as clustering cells, identifying cluster marker genes, detecting differentially expressed genes between experimental conditions, and importantly outputting integrated gene expression values for downstream data analysis.


#### Install dependent packages:
```
install.packages(c("Matrix", "matrixStats", "irlba", "doParallel", "foreach", "data.table", "Rtsne", "umap", "MASS", "pbmcapply", "Rcpp", "RcppEigen", "densityClust", "FNN", "igraph", "RColorBrewer", "ggplot2", "gridExtra", "pheatmap"))
```

#### Install RISC:
```
install_github("https://github.com/bioinfoDZ/RISC.git")
```
The RISC package can also be downloaded and installed mannually
<a href="https://github.com/bioinfoDZ/RISC/blob/master/RISC_1.0.tar.gz" download="RISC_1.0.tar.gz">Link</a>
```
install.packages("/Path/to/RISC_1.0.tar.gz", repos = NULL, type = "source")
```


### vignettes
Here we provide a vignettes which shows the key steps in analyzing example scRNA-seq datasets from the basal or squamous carcinoma patients before and after anti-PD-1 therapy (GSE123813). 
<a href="https://github.com/bioinfoDZ/RISC/blob/master/GSE123813_Vignette.pdf" download="GSE123813_Vignette.pdf">Link</a>

We also provide an example how to convert Seurat object to RISC object (to use the new features, please reinstall RISC package)
<a href="https://github.com/bioinfoDZ/RISC/blob/master/Seurat_to_RISC.pdf" download="Seurat_to_RISC.pdf">Link</a>

#### Notice, RISC package is developed in R (v3.6.3), we test this vignette in the same R version.


#### Contents:
(1) RISC package: "RISC_1.0.tar.gz" <br />
(2) Vignette for GSE123813: "GSE123813_Vignette.pdf" <br />
(3) GSE123813 directory contains the informaiton of cell-type, patients and treatment.
file position, "/GSE123813/Raw_Data/bcc_annotation.tsv"


