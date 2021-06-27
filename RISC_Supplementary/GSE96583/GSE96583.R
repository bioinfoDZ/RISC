library(RISC)
library(irlba)
library(umap)
library(Matrix)
library(matrixStats)
library(ggplot2)
library(RColorBrewer)


PATH = "/Path to the data/GSE96583_PBMC"

###################################################################################
### Preparing RISC Objects ###
###################################################################################
# Input Data
data1 = read10Xgenomics(paste0(PATH, "/PBMC_Raw_Counts/Control"))
data2 = read10Xgenomics(paste0(PATH, "/PBMC_Raw_Counts/Stimulate"))

# Processing
process0 <- function(obj0){
  obj0 = scFilter(obj0, min.UMI = 200, min.gene = 200, max.UMI = 8000, min.cell = 3)
  obj0 = scNormalize(obj0)
  obj0 = scDisperse(obj0)
  print(length(obj0@vargene))
  return(obj0)
}

data1 = process0(data1)
data2 = process0(data2)


###################################################################################
### Uncorrected Data ###
###################################################################################
## Plot
Interplot <- function(m0, color0, shape0, group0 = "Type"){
  if(group0 == "Type"){
    g0 = ggplot(m0, aes(UMAP1, UMAP2)) + 
      geom_point(aes(color = CellType, shape = Group), size = 0.5, alpha = 1) + 
      scale_color_manual(values = color0) + 
      scale_shape_manual(values = shape0) + 
      theme_bw(base_line_size = 0) + 
      labs(color = "Cell-Type", shape = "Batch") + 
      theme(plot.title = element_text(hjust = 0.5, size = 16)) + 
      guides(color = guide_legend(override.aes = list(size = 6)), shape = guide_legend(override.aes = list(size = 5)))
  } else {
    g0 = ggplot(m0, aes(UMAP1, UMAP2)) + 
      geom_point(aes(color = Group, shape = Group), size = 0.5, alpha = 1) + 
      scale_color_manual(values = c("#FB8072", "#80B1D3")) + 
      scale_shape_manual(values = shape0) + 
      theme_bw(base_line_size = 0) + 
      labs(color = "Cell-Type", shape = "Batch") + 
      theme(plot.title = element_text(hjust = 0.5, size = 16)) + 
      guides(color = guide_legend(override.aes = list(size = 6)), shape = guide_legend(override.aes = list(size = 5)))
  }
  print(g0)
}

# Raw PCs
logmat1 = as.matrix(data1@assay$logcount)
logmat2 = as.matrix(data2@assay$logcount)
gene0 = intersect(rownames(logmat1), rownames(logmat2))
logmat0 = cbind(logmat1[gene0,], logmat2[gene0,])
pca0 = irlba(logmat0, nv = 50, center = T)$v
umap0 = umap(pca0[,1:16])$layout

# Plot
m0 = data.frame(UMAP1 = umap0[,1], UMAP2 = umap0[,2], Group = c(rep("IFN-", ncol(logmat1)), rep("IFN+", ncol(logmat2))))
shape0 = c(1, 8)
Interplot(m0, color0, shape0, group0 = "Set")


###################################################################################
### Data Integration with Full Cells ###
###################################################################################
## Integration
var0 = read.table(paste0(PATH, "/Var0_Ori.tsv"), sep = "\t", header = F, stringsAsFactors = F)
var0 = var0$V1
data3 = list(data1, data2)
InPlot(data3, var.gene = var0, nPC = 20)
data3 = scMultiIntegrate(data3, eigens = 15, var.gene = var0, add.Id = c('IFN-', 'IFN+'), adjust = FALSE, ncore = 4)
data3 = scUMAP(data3, npc = 18, use = 'PLS')

DimPlot(data3, slot = "cell.umap", colFactor = "Set", Colors = c("firebrick2", "skyblue"), size = 0.2)


###################################################################################
### Subset the Cells ###
###################################################################################
Anno0 = read.table(paste0(PATH, "/Anno0.tsv"), sep = "\t", header = T, stringsAsFactors = F)
cell1 = Anno0$scBarcode[Anno0$Set == "IFN-"]
cell2 = Anno0$scBarcode[Anno0$Set == "IFN+"]

data1 = SubSet(data1, cells = cell1)
data2 = SubSet(data2, cells = cell2)
data1 = process0(data1)
data2 = process0(data2)

FilterPlot(data1)
FilterPlot(data2)


###################################################################################
### Data Integration with Selected Cells ###
###################################################################################
## Integration
var0 = read.table(paste0(PATH, "/Var0.tsv"), sep = "\t", header = F, stringsAsFactors = F)
var0 = var0$V1
data3 = list(data1, data2)
InPlot(data3, var.gene = var0, nPC = 20)
data3 = scMultiIntegrate(data3, eigens = 16, var.gene = var0, add.Id = c('IFN-', 'IFN+'), adjust = TRUE, ncore = 4)
data3 = scUMAP(data3, npc = 16, use = 'PLS')
DimPlot(data3, slot = "cell.umap", colFactor = "Set", Colors = c("firebrick2", "skyblue"), size = 0.2)

# Cell Type
Anno0 = read.table(paste0(PATH, "/Anno0.tsv"), sep = "\t", header = T, stringsAsFactors = F)
data3@coldata$Type0 = Anno0$Type0
DimPlot(data3, slot = "cell.umap", colFactor = "Type0", size = 0.2, label = T)

# Marker Genes
DimPlot(data3, slot = "cell.umap", genes = "CD8A", size = 0.2)
DimPlot(data3, slot = "cell.umap", genes = "FYN", size = 0.2)
DimPlot(data3, slot = "cell.umap", genes = "CACYBP", size = 0.2)
DimPlot(data3, slot = "cell.umap", genes = "IL1B", size = 0.2)
DimPlot(data3, slot = "cell.umap", genes = "HES4", size = 0.2)
DimPlot(data3, slot = "cell.umap", genes = "GNLY", size = 0.2)
DimPlot(data3, slot = "cell.umap", genes = "AES", size = 0.2)
DimPlot(data3, slot = "cell.umap", genes = "CD79A", size = 0.2)
DimPlot(data3, slot = "cell.umap", genes = "CD79B", size = 0.2)
DimPlot(data3, slot = "cell.umap", genes = "SMPDL3A", size = 0.2)
DimPlot(data3, slot = "cell.umap", genes = "MIR155HG", size = 0.2)
DimPlot(data3, slot = "cell.umap", genes = "CKB", size = 0.2)
DimPlot(data3, slot = "cell.umap", genes = "PPBP", size = 0.2)


