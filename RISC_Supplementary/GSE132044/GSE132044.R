library(RISC)
library(Matrix)
library(irlba)
library(umap)
library(ggplot2)
library(RColorBrewer)


PATH = "/Path to the data/GSE132044"


###################################################################################
### Prepare Data ###
###################################################################################
## All Cells
mat0 = readMM(file = paste0(PATH, "/Raw_Data/counts.umi.mtx"))
cell0 = read.table(file = paste0(PATH, "/Raw_Data/cells.umi.txt"), sep = "\t", header = F, stringsAsFactors = F)
colnames(cell0) = "Mix"
cell0$Sort = sapply(cell0$Mix, function(x){strsplit(x, "_", fixed = T)[[1]][1]})
cell0$Platform = sapply(cell0$Mix, function(x){strsplit(x, "_", fixed = T)[[1]][2]})
gene0 = read.table(file = paste0(PATH, "/Raw_Data/genes.umi.txt"), sep = "\t", header = F, stringsAsFactors = F)
colnames(gene0) = "Mix"
gene0$Ensembl = sapply(gene0$Mix, function(x){strsplit(x, "_", fixed = T)[[1]][1]})
gene0$Symbol = sapply(gene0$Mix, function(x){strsplit(x, "_", fixed = T)[[1]][2]})
keep = !duplicated(gene0$Symbol)
gene0 = gene0[keep,]
mat0 = mat0[keep,]
colnames(mat0) = cell0$Mix
rownames(mat0) = gene0$Symbol


## PBMC ("10xChromiumv2A", "Drop", "inDrops")
# 10X Chromium V2
keep = cell0$Sort == "pbmc2" & cell0$Platform == "10xChromiumv2"
cell1 = cell0[keep,]
cell1$Barcode = sapply(cell1$Mix, function(x){strsplit(x, "_", fixed = T)[[1]][3]})
rownames(cell1) = cell1$Barcode
mat1 = mat0[,keep]
keep = rowSums(as.matrix(mat1) > 0) > 0
mat1 = mat1[keep,]
colnames(mat1) = rownames(cell1)
gene0 = data.frame(Symbol = rownames(mat1), row.names = rownames(mat1))
dat1 = readscdata(mat1, cell1, gene0, is.filter = F)

# Drop-seq
# keep = cell0$Sort == "pbmc2" & cell0$Platform == "Drop" & cell0$Mix %in% Ann0$Mix
keep = cell0$Sort == "pbmc2" & cell0$Platform == "Drop"
cell1 = cell0[keep,]
cell1$Barcode = sapply(cell1$Mix, function(x){strsplit(x, "_", fixed = T)[[1]][4]})
rownames(cell1) = cell1$Barcode
mat1 = mat0[,keep]
keep = rowSums(as.matrix(mat1) > 0) > 0
mat1 = mat1[keep,]
colnames(mat1) = rownames(cell1)
gene0 = data.frame(Symbol = rownames(mat1), row.names = rownames(mat1))
dat2 = readscdata(mat1, cell1, gene0, is.filter = F)

# inDrops-seq
# keep = cell0$Sort == "pbmc2" & cell0$Platform == "inDrops" & cell0$Mix %in% Ann0$Mix
keep = cell0$Sort == "pbmc2" & cell0$Platform == "inDrops"
cell1 = cell0[keep,]
cell1$Barcode = sapply(cell1$Mix, function(x){strsplit(x, "_", fixed = T)[[1]][4]})
rownames(cell1) = cell1$Barcode
mat1 = mat0[,keep]
keep = rowSums(as.matrix(mat1) > 0) > 0
mat1 = mat1[keep,]
colnames(mat1) = rownames(cell1)
gene0 = data.frame(Symbol = rownames(mat1), row.names = rownames(mat1))
dat3 = readscdata(mat1, cell1, gene0, is.filter = F)


###################################################################################
### RISC Objects ###
###################################################################################
process0 <- function(obj0){
  obj0 = scFilter(obj0, min.UMI = 100, max.UMI = 10000, min.gene = 100, min.cell = 3)
  obj0 = scNormalize(obj0, ncore = 4)
  obj0 = scDisperse(obj0)
  print(length(obj0@vargene))
  return(obj0)
}

dat1 = process0(dat1)
dat2 = process0(dat2)
dat3 = process0(dat3)


###################################################################################
### Uncorrect Data ###
###################################################################################
set.seed(123)
var0 = Reduce(intersect, list(dat1@rowdata$Symbol, dat2@rowdata$Symbol, dat3@rowdata$Symbol))
logmat1 = as.matrix(dat1@assay$logcount)[var0,]
logmat2 = as.matrix(dat2@assay$logcount)[var0,]
logmat3 = as.matrix(dat3@assay$logcount)[var0,]
logmat0 = cbind(logmat3, logmat1, logmat2)
pca0 = irlba(logmat0, nv = 50)$v
umap0 = umap(pca0[,1:15])$layout

ann0 = read.table(file = paste0(PATH, "/GSE132044_Anno.tsv"), sep = "\t", header = T, stringsAsFactors = F)
m0 = data.frame(UMAP1 = umap0[,1], UMAP2 = umap0[,2], Plat = ann0$Set, Type = ann0$Type)
m0$Platform = factor(m0$Plat, levels = c("inDrops", "10X-Genomics", "Drop-seq"), labels = c("inDrops", "10xGenomics", "Drop-seq"))
m0$Type = factor(m0$Type, levels = c("B", "Platelet", "pDC", "CD16+ Mono", "CD4 Memory T", "secretory B", "CD4 Naive T", "CD8 T", "DC", "NK", "CD14+ Mono"))
color0 = brewer.pal(11, "Spectral")
color1 = c("#D9D9D9", "#FDB462", "#BC80BD")

ggplot(m0, aes(UMAP1, UMAP2)) + 
  geom_point(aes(color = Platform, shape = Platform), size = 1, alpha = 1) + 
  scale_color_manual(values = color1) + 
  scale_shape_manual(values = c(1, 2, 6)) + 
  theme_bw(base_line_size = 0) + 
  labs(x = 'UMAP-1', y = 'UMAP-2', color = 'Platforms', shape = 'Set') + 
  theme(plot.title = element_text(hjust = 0.5, size = 16)) + 
  guides(colour = guide_legend(override.aes = list(size = 5), ncol = 1), shape = guide_legend(override.aes = list(size = 5), ncol = 1))

ggplot(m0, aes(UMAP1, UMAP2)) + 
  geom_point(aes(color = Type, shape = Platform), size = 1, alpha = 1) + 
  scale_color_manual(values = color0) + 
  scale_shape_manual(values = c(1, 2, 6)) + 
  theme_bw(base_line_size = 0) + 
  labs(x = 'UMAP-1', y = 'UMAP-2', color = 'Platforms', shape = 'Set') + 
  theme(plot.title = element_text(hjust = 0.5, size = 16)) + 
  guides(colour = guide_legend(override.aes = list(size = 5), ncol = 1), shape = guide_legend(override.aes = list(size = 5), ncol = 1))


###################################################################################
### RISC Integration ###
###################################################################################
## PBMC2
set.seed(123)
data0 = list(dat3, dat1, dat2)
var0 = read.table(file = paste0(PATH, "/var.tsv"), sep = "\t", header = F, stringsAsFactors = F)
var0 = var0$V1
InPlot(data0, var.gene = var0, ncore = 4, minPC = 11, nPC = 15)
data0 = scMultiIntegrate(data0, eigens = 14, var.gene = var0, add.Id = c("inDrops", "10X-Genomics", "Drop-seq"), ncore = 4)
data0 = scUMAP(data0, npc = 15, use = "PLS")
ann0 = read.table(file = paste0(PATH, "/GSE132044_Anno.tsv"), sep = "\t", header = T, stringsAsFactors = F)
data0@coldata$Type = ann0$Type

DimPlot(data0, slot = "cell.umap", colFactor = "Set", size = 0.5)
DimPlot(data0, slot = "cell.umap", colFactor = "Type", size = 0.5, label = T)
DimPlot(data0, slot = "cell.umap", genes = "MS4A1", size = 0.2)
DimPlot(data0, slot = "cell.umap", genes = "GNLY", size = 0.2)
DimPlot(data0, slot = "cell.umap", genes = "RCAN3", size = 0.2)
DimPlot(data0, slot = "cell.umap", genes = "CD8A", size = 0.2)
DimPlot(data0, slot = "cell.umap", genes = "GZMK", size = 0.2)
DimPlot(data0, slot = "cell.umap", genes = "CD14", size = 0.2)
DimPlot(data0, slot = "cell.umap", genes = "FCGR3A", size = 0.2)
DimPlot(data0, slot = "cell.umap", genes = "FCER1A", size = 0.2)
DimPlot(data0, slot = "cell.umap", genes = "PPBP",  size = 0.2)
DimPlot(data0, slot = "cell.umap", genes = "SERPINF1", size = 0.2)
DimPlot(data0, slot = "cell.umap", genes = "JCHAIN", size = 0.2)


