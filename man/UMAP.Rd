% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Reduce_Dimension.R
\name{scUMAP}
\alias{scUMAP}
\title{Dimension Reduction.}
\usage{
scUMAP(object, npc = 20, embedding = 2, use = "PCA", seed = 123)
}
\arguments{
\item{object}{RISC object: a framework dataset.}

\item{npc}{The number of the PCs (or the PLS) using for UMAP, the default is 20, 
but need to be modified by the users. The PCA for individual dataset, while PLS 
for the integrated data.}

\item{embedding}{The number of components UMAP output.}

\item{use}{What components used for UMAP: PCA or PLS.}

\item{seed}{The random seed to keep tSNE result consistent.}
}
\value{
RISC single cell dataset, the DimReduction slot.
}
\description{
The UMAP is calculated based on the eigenvectors of single cell dataset, and the 
user can select the eigenvectors manually. Of note, the selected eigenvectors 
directly affect UMAP values. 
For the integrated data (the result of "scMultiIntegrate" funciton), RISC utilizes
the PCR output "PLS" to calculate the UMAP, therefore, the user has to input "PLS"
in "use = ", instead of the default parameter "PCA".
}
\examples{
# RISC object
obj0 = raw.mat[[3]]
obj0 = scPCA(obj0, npc = 10)
obj0 = scUMAP(obj0, npc = 3)
DimPlot(obj0, slot = "cell.umap", colFactor = 'Group', size = 2)
}
\references{
Becht et al., Nature Biotech. (2018)
}
