## RISC


### Overview
Integrated analysis of single cell RNA-sequencing (scRNA-seq) data from multiple batches or studies is often necessary in order to learn functional changes in cellular states upon experimental perturbations or cell type relationships in a developmental lineage. Here we introduce a new algorithm (RPCI) that uses the gene-eigenvectors from a reference dataset to establish a global frame for projecting all datasets, with a clear advantage in preserving genuine gene expression differences in matching cell types between samples, such as those present in cells at distinct developmental stages or in perturbated vs control studies. This R package “RISC” (Robust Integration of Sinlgle Cell RNA-seq data implements the RPCI algorithm, with additional functions for scRNA-seq analysis, such as clustering cells, identifying cluster marker genes, detecting differentially expressed genes between experimental conditions, and importantly outputting integrated gene expression values for downstream data analysis.


#### RISC v1.6 update
This version mainly solves the problems that are caused by dependent package updates

#### RISC v1.5 update
Changes from the last release (v1.0) <br />
(1) Replace dependent "RcppEigen" with "RcppArmadillo", fully support sparse matrix in core functions. <br />
(2) Replace dependent "pbmcapply" with "pbapply" <br />
(3) Optimize "scMultiIntegrate" function and reduce memory-consuming; the new RISC release can support integration of datasets with >1.5 million cells and 10,000 genes. <br />
(4) When data integration, all the genes expressed in indivudal datasets will be reserved in the integrated data. The genes shared expressed across samples will be labeled in "rowdata" of RISC object. <br />
(5) Convert "logcounts" in the integrated RISC object "object@assay$logcount" from a large matrix to a list including multiple logcounts matrices, each corrected matrix for the corresponding individual data sets. To output full integrated matrix, mat0 = do.call(cbind, object@assay$logcount) <br />
(6) Change function name "readscdata" -> "readsc" <br />
(7) Change function name "read10Xgenomics" -> "read10X_mtx" <br />
(8) Parameter names in some functions are changed. <br />

Added new functions <br />
(1) In "scMarker" and "AllMarker" functions, add Wilcoxon Rank Sum and Signed Rank model. <br />
(2) In "scMarker", "AllMarker" and "scDEG" functions, add pseudo-cell (bin cells to generate meta-cells) option to detect marker genes. <br />
(3) Add "slot" parameter in "DimPlot" function, external dimension reduction results can be added in RISC object, e.g. add phate results (phate0) to RISC object obj0@DimReduction$cell.phate = phate0; DimPlot(obj0, slot = "cell.phate", colFactor = 'Group', size = 2, label = TRUE) <br />
(4) Add "read10X_h5" function for 10X Genomics h5 file. <br />

Removed old functions <br />
(1) delete "readHTSeqdata" function. <br />


#### Install dependent packages:
```
install.packages(c("Matrix", "sparseMatrixStats", "Matrix.utils", "irlba", "doParallel", "foreach", "data.table", "Rtsne", "umap", "MASS", "pbapply", "Rcpp", "RcppArmadillo", "densityClust", "FNN", "igraph", "RColorBrewer", "ggplot2", "gridExtra", "pheatmap", "hdf5r"))
```

#### Install RISC:
```
install_github("https://github.com/bioinfoDZ/RISC.git")
```
The RISC package can also be downloaded and installed mannually
<a href="https://github.com/bioinfoDZ/RISC/blob/master/RISC_1.6.0.tar.gz" download="RISC_1.6.0.tar.gz">Link</a>
```
install.packages("/Path/to/RISC_1.6.0.tar.gz", repos = NULL, type = "source")
```


### vignettes
Here we provide a vignettes which shows the key steps in analyzing example scRNA-seq datasets from the basal or squamous carcinoma patients before and after anti-PD-1 therapy (GSE123813). 

#### RISC   v1.0 <a href="https://github.com/bioinfoDZ/RISC/blob/master/GSE123813_Vignette_RISC_v1.0.pdf" download="GSE123813_Vignette_RISC_v1.0.pdf">Link</a>
#### RISC   v1.6 <a href="https://github.com/bioinfoDZ/RISC/blob/master/GSE123813_Vignette_RISC_v1.6.pdf" download="GSE123813_Vignette_RISC_v1.6.pdf">Link</a>

We also provide an example how to convert Seurat object to RISC object (to use the new features, please reinstall RISC package)

#### RISC v1.0   <a href="https://github.com/bioinfoDZ/RISC/blob/master/Seurat_to_RISC_RISC_v1.0.pdf" download="Seurat_to_RISC_RISC_v1.0.pdf">Link</a>

#### Notice, RISC v1.6 package is developed in R (v4.2.2), we test this vignette in the same R version.


#### Contents:
(1) RISC package: "RISC_1.6.0.tar.gz" <br />
(2) Vignette for GSE123813: "GSE123813_Vignette_RISC_v1.6.pdf" <br />
(3) GSE123813 directory contains the informaiton of cell-type, patients and treatment.
file position, "/GSE123813/Raw_Data/bcc_annotation.tsv" <br />

Old RISC version: "RISC_1.0.tar.gz"


### Citation:
Liu Y, Tao W, Zhou B, Zheng D (2021) Robust integration of multiple single-cell RNA sequencing datasets using a single reference space. 
 <a href="https://doi.org/10.1038/s41587-021-00859-x">Nat Biotechnol 39(7):877-884.</a>
