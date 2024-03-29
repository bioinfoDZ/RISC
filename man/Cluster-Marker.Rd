% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ClusterMarker.R
\name{scMarker}
\alias{scMarker}
\title{Find Cluster Markers}
\usage{
scMarker(
  object,
  cluster = 1,
  positive = TRUE,
  frac = 0.25,
  log2FC = 0.5,
  Padj = 0.05,
  latent.factor = NULL,
  method = "QP",
  min.cells = 10,
  ncore = 1
)
}
\arguments{
\item{object}{RISC object: a framework dataset.}

\item{cluster}{Select the cluster that we want to detect cluster marker genes.}

\item{positive}{Whether only output the cluster markers with positive log2FC.}

\item{frac}{A fraction cutoff, the marker genes expressed at least a 
cutoff fraction of all the cells.}

\item{log2FC}{The cutoff of log2 Fold-change for differentially expressed marker 
genes.}

\item{Padj}{The cutoff of the adjusted P-value.}

\item{latent.factor}{The latent factor from coldata, which represents number 
values or factors, and only one latent factor can be inputed.}

\item{method}{Which method is used to identify cluster markers, three options: 'NB' 
for Negative Binomial model, 'QP' for QuasiPoisson model, and 'Wilcox' for Wilcoxon 
Rank Sum and Signed Rank model.}

\item{min.cells}{The minimum cells for each cluster to calculate marker genes.}

\item{ncore}{The multiple cores for parallel calculating.}
}
\description{
This is the basic function in RISC, it can identify the cluster markers by
comparing samples in the selected cluster to the samples of the rest clusters.
Therefore, it is possible one gene labeled as a marker in more than one clusters.
Two methods are employed in RISC, one is based on Negative Binomial model 
while the other using QuasiPoisson model.
}
\details{
Because log2 cannot handle counts with value 0, we use log1p to calculate average 
values of counts and log2 to format fold-change.
}
\examples{
# RISC object
obj0 = raw.mat[[3]]
obj0 = scPCA(obj0, npc = 10)
obj0 = scUMAP(obj0, npc = 3)
obj0 = scCluster(obj0, slot = "cell.umap", k = 3, method = 'density')
DimPlot(obj0, slot = "cell.umap", colFactor = 'Cluster', size = 2)
marker1 = scMarker(obj0, cluster = 1, method = 'QP', min.cells = 3)
}
\references{
Paternoster et al., Criminology (1997)

Berk et al., Journal of Quantitative Criminology (2008)

Liu et al., Nature Biotech. (2021)
}
