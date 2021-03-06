% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Integrating.R
\name{SIMPLS}
\alias{SIMPLS}
\title{Integration Algorithm SIMPLS}
\usage{
SIMPLS(X, Y, npcs = 10, seed = 123)
}
\arguments{
\item{X}{The reference matrix, row for genes and column for cells.}

\item{Y}{The target matrix, row for genes and column for cells.}

\item{npcs}{The number of the PCs used for data integration.}

\item{seed}{The random seed to keep consistent result.}
}
\description{
The partial least square (PLS) with SIMPLS algorithm is an extension of the 
multiple linear regression model and considered as bilinear factor models. 
Instead of embedding the reference and target matrices into a hyperplane 
of maximum variance, the PLS utilizes a linear regression to project the 
reference and target matrices into a new place. The SIMPLS algorithm provides 
the regularization procedure for PLS. The matrices need to be centered before 
SIMPLS integraton.
}
\examples{
## SIMPLS with two matrices
mat0 = as.matrix(raw.mat[[1]])
coldata0 = as.data.frame(raw.mat[[2]])
coldata1 = coldata0[coldata0$Batch0 == 'Batch1',]
coldata2 = coldata0[coldata0$Batch0 == 'Batch4',]
mat1 = mat0[,rownames(coldata1)]
mat2 = mat0[,rownames(coldata2)]
SIM0 = SIMPLS(mat1, mat2, npcs = 4)
}
\references{
De-Jong et al. (1993)
}
