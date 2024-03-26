## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(RISC)
library(RColorBrewer)

data("raw.mat")
mat0 = raw.mat[[1]]
coldata0 = raw.mat[[2]]

coldata1 = coldata0[coldata0$Batch0 == 'Batch1',]
coldata2 = coldata0[coldata0$Batch0 == 'Batch2',]
coldata3 = coldata0[coldata0$Batch0 == 'Batch3',]
coldata4 = coldata0[coldata0$Batch0 == 'Batch4',]
coldata5 = coldata0[coldata0$Batch0 == 'Batch5',]
coldata6 = coldata0[coldata0$Batch0 == 'Batch6',]
mat1 = mat0[,rownames(coldata1)]
mat2 = mat0[,rownames(coldata2)]
mat3 = mat0[,rownames(coldata3)]
mat4 = mat0[,rownames(coldata4)]
mat5 = mat0[,rownames(coldata5)]
mat6 = mat0[,rownames(coldata6)]

## -----------------------------------------------------------------------------
sce1 = readsc(count = mat1, cell = coldata1, gene = data.frame(Symbol = rownames(mat1), row.names = rownames(mat1)), is.filter = FALSE)
sce2 = readsc(count = mat2, cell = coldata2, gene = data.frame(Symbol = rownames(mat2), row.names = rownames(mat2)), is.filter = FALSE)
sce3 = readsc(count = mat3, cell = coldata3, gene = data.frame(Symbol = rownames(mat3), row.names = rownames(mat3)), is.filter = FALSE)
sce4 = readsc(count = mat4, cell = coldata4, gene = data.frame(Symbol = rownames(mat4), row.names = rownames(mat4)), is.filter = FALSE)
sce5 = readsc(count = mat5, cell = coldata5, gene = data.frame(Symbol = rownames(mat5), row.names = rownames(mat5)), is.filter = FALSE)
sce6 = readsc(count = mat6, cell = coldata6, gene = data.frame(Symbol = rownames(mat6), row.names = rownames(mat6)), is.filter = FALSE)

## -----------------------------------------------------------------------------
process0 <- function(obj0){
  
  # filter cells
  obj0 = scFilter(obj0, min.UMI = 0, max.UMI = Inf, min.gene = 10, min.cell = 3, is.filter = FALSE)
  
  # normalize data
  obj0 = scNormalize(obj0, method = 'robust')
  
  # find highly variable genes
  obj0 = scDisperse(obj0)
  
  # here replace highly variable genes by all the genes for integraton
  obj0@vargene = rownames(sce1@rowdata)
  
  return(obj0)
  
}

sce1 = process0(sce1)
sce2 = process0(sce2)
sce3 = process0(sce3)
sce4 = process0(sce4)
sce5 = process0(sce5)
sce6 = process0(sce6)

## -----------------------------------------------------------------------------
set.seed(1)
var.genes = rownames(sce1@assay$count)
pcr0 = list(sce1, sce2, sce3, sce4, sce5, sce6)
pcr0 = scMultiIntegrate(pcr0, eigens = 9, var.gene = var.genes, align = 'OLS', npc = 15)
# pcr0 = scLargeIntegrate(pcr0, var.gene = var.genes, align = 'Predict', npc = 8)
pcr0 =scUMAP(pcr0, npc = 9, use = 'PLS', dist = 0.001, neighbors = 15)

## -----------------------------------------------------------------------------
pcr0@coldata$Group = factor(pcr0@coldata$Group0, levels = c('Group1', 'Group2', 'Group2*', 'Group3', 'Group3*'), labels = c("a", "b", "b'", "c", "c'"))
pcr0@coldata$Set0 = factor(pcr0@coldata$Set, levels = c('Set1', 'Set2', 'Set3', 'Set4', 'Set5', 'Set6'), labels = c('Set1 rep.1', 'Set2 rep.1', 'Set3 rep.1', 'Set1 rep.2', 'Set2 rep.2', 'Set3 rep.2'))
pcr0 = scCluster(pcr0, slot = "cell.umap", k = 4, method = "density", dc = 0.3)

## ----fig.show="hold", out.width="48", fig.dim=c(7, 5)-------------------------
DimPlot(pcr0, colFactor = 'Set0', size = 2)
DimPlot(pcr0, colFactor = 'Group', size = 2, Colors = brewer.pal(5, "Set1"))
DimPlot(pcr0, colFactor = 'Cluster', size = 2, Colors = brewer.pal(6, "Dark2"))

## -----------------------------------------------------------------------------
sessionInfo()

