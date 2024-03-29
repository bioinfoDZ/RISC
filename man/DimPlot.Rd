% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Graph.R
\name{tSNEPlot}
\alias{tSNEPlot}
\alias{DimPlot}
\title{Dimension Reduction Plots}
\usage{
DimPlot(
  object,
  slot = "cell.umap",
  colFactor = NULL,
  genes = NULL,
  legend = TRUE,
  Colors = NULL,
  size = 0.5,
  Alpha = 0.8,
  plot.ncol = NULL,
  exp.range = NULL,
  exp.col = "firebrick2",
  label = FALSE,
  adjust.label = 0.25,
  label.font = 5
)
}
\arguments{
\item{object}{RISC object: a framework dataset.}

\item{slot}{The dimension_reduction slot for drawing the plots. The default is 
"cell.umap" under RISC object "DimReduction" item for UMAP plot, but the customer 
can add new dimension_reduction method under DimReduction and use it.}

\item{colFactor}{Use the factor (column name) in the coldata to make a 
dimension_reduction plot, but each time only one column name can be inputted.}

\item{genes}{Use the gene expression values (gene symbol) to make dimension 
reduction plot, each time more than one genes can be inputted.}

\item{legend}{Whether a legend shown at dimension_reduction plot.}

\item{Colors}{The users can use their own colors (color vector). The default of 
the "tSNEPlot" funciton will assign colors automatically.}

\item{size}{Choose the size of dots at dimension_reduction plot, the default 
size is 0.5.}

\item{Alpha}{Whether show transparency of individual points, the default is 0.8.}

\item{plot.ncol}{If the users input more than one genes, the arrangement of 
multiple dimension_reduction plot depends on this parameter.}

\item{exp.range}{The gene expression cutoff for plot, e.g. "c(0, 1.5)" for 
expression level between 0 and 1.5.}

\item{exp.col}{The gradient color for gene expression.}

\item{label}{Whether label the clusters or cell populations in the plot.}

\item{adjust.label}{The adjustment of the label position.}

\item{label.font}{The font size for the label.}
}
\description{
The Dimension Reduction plots are widespread in scRNA-seq data analysis. Here, 
the "DimPlot" function not only can make plots for factor labels of individual 
cells but also can show gene expression values of each cell.
}
\examples{
# RISC object
obj0 = raw.mat[[3]]
obj0 = scPCA(obj0, npc = 10)
obj0 = scUMAP(obj0, npc = 3)
DimPlot(obj0, slot = "cell.umap", colFactor = 'Group', size = 2, label = TRUE)
DimPlot(obj0, genes = c('Gene718', 'Gene325', 'Gene604'), size = 2)
}
\references{
Wickham, H. (2016)

Auguie, B. (2015)
}
