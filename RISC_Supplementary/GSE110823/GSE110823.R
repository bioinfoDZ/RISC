library(RISC)
library(Matrix)
library(R.matlab)


PATH = "/Path to the data/GSE110823"

###################################################################################
### Prepare Data ###
###################################################################################
## All Cells
dat0 = readMat(paste0(PATH, "/GSM3017261_150000_CNS_nuclei.mat"))
mat0 = dat0$DGE
coldata0 = data.frame(Barcode = paste0("Cell-", dat0$barcodes[1,]), Organ = dat0$sample.type[,1], Type = dat0$cluster.assignment[,1])
coldata0$Type = sapply(coldata0$Type, function(x){gsub(" ", "-", x, fixed = T)})
gene0 = data.frame(Symbol = dat0$genes[,1])
gene0$Symbol = sapply(gene0$Symbol, function(x){gsub(" ", "", x, fixed = T)})
colnames(mat0) = rownames(gene0) = gene0$Symbol
rownames(mat0) = rownames(coldata0) = coldata0$Barcode

# P2 Brain
keep = coldata0$Organ == "p2_brain " & !coldata0$Type %in% c("53-Unresolved------------------", "54-Unresolved-Kcng1------------")
coldata1 = coldata0[keep,]
mat1 = t(mat0[keep,])
keep = Matrix::rowSums(mat1 > 0) > 0
mat1 = mat1[keep,]
rowdata1 = data.frame(Symbol = rownames(mat1), row.names = rownames(mat1))
dat1 = readscdata(mat1, coldata1, rowdata1, is.filter = F)

# P11 Brain
keep = coldata0$Organ == "p11_brain" & !coldata0$Type %in% c("53-Unresolved------------------", "54-Unresolved-Kcng1------------")
coldata1 = coldata0[keep,]
mat1 = t(mat0[keep,])
keep = Matrix::rowSums(mat1 > 0) > 0
mat1 = mat1[keep,]
rowdata1 = data.frame(Symbol = rownames(mat1), row.names = rownames(mat1))
dat2 = readscdata(mat1, coldata1, rowdata1, is.filter = F)


###################################################################################
### RISC Objects ###
###################################################################################
process0 <- function(obj0){
  obj0 = scFilter(obj0, min.UMI = 500, max.UMI = 5000, min.gene = 200, min.cell = 3)
  obj0 = scNormalize(obj0)
  obj0 = scDisperse(obj0)
  print(length(obj0@vargene))
  return(obj0)
}

dat1 = process0(dat1)
dat2 = process0(dat2)


###################################################################################
### Uncorrect ###
###################################################################################
library(irlba)
library(Rtsne)
library(ggplot2)
library(RColorBrewer)

set.seed(1987)
gene0 = intersect(dat1@rowdata$Symbol, dat2@rowdata$Symbol)
ann0 = read.table(file = paste0(PATH, "/GSE110823_ann.tsv"), sep = "\t", header = T, stringsAsFactors = F)
logmat0 = cbind(dat2@assay$logcount[gene0, ann0$Barcode[ann0$Set == "P11"]], dat1@assay$logcount[gene0, ann0$Barcode[ann0$Set == "P2"]])
pca0 = irlba(logmat0, nv = 50)$v
tsne0 = Rtsne(pca0)$Y
m0 = data.frame(tSNE1 = tsne0[,1], tSNE2 = tsne0[,2], Set = ann0$Set, CellType = ann0$CellType)
m0$Set = factor(m0$Set, levels = c("P2", "P11"))
m0$CellType = factor(m0$CellType, levels = paste0("C", 1:59))
color0 = c("#1B9E77", "#288E96", "#367EB6", "#419486", "#4BAC50", "#579055", "#636C63", "#667F49", "#669E26", "#6B9254", "#72789C", "#8064AD", "#9154A5", "#98649F", "#98899B", "#9C867A", "#A26643", "#A65D25", "#A66D20", "#B07117", "#C9650A", "#DB5206", "#E03013", "#E42F18", "#E5760B", "#E69B12", "#E65C54", "#E8318E", "#F05BA8", "#F780B3", "#FB7F56", "#FF8201", "#FFC01A", "#FFFF33", "#9E0142", "#AF1446", "#C0274A", "#D13A4E", "#DC494C", "#E65848", "#F06744", "#F57948", "#F88D51", "#FBA15B", "#FDB466", "#FDC373", "#FDD380", "#FEE18E", "#FEEB9E", "#FEF5AE", "#FFFFBF", "#F7FBB2", "#EFF8A6", "#E7F59A", "#D7EF9B", "#C4E79E", "#B2E0A2", "#9ED7A4", "#88CFA4", "#72C7A4", "#5FBAA8", "#4FA8AF", "#3F96B7", "#3484BB", "#4272B2", "#5060AA", "#5E4FA2")

ggplot(m0, aes(tSNE1, tSNE2)) + 
  geom_point(aes(color = Set), size = 0.5) + 
  scale_color_manual(values = c("#FB8072", "#80B1D3")) + 
  theme_bw(base_size = 12, base_line_size = 0) + 
  labs(color = "Set") + 
  guides(color = guide_legend(override.aes = list(size = 8), ncol = 1))
ggplot(m0, aes(tSNE1, tSNE2)) + 
  geom_point(aes(color = CellType), size = 0.5) + 
  scale_color_manual(values = color0) + 
  theme_bw(base_size = 12, base_line_size = 0) + 
  labs(color = "Cell Type") + 
  theme(legend.text = element_text(size = 16), legend.title = element_text(size = 20)) + 
  guides(color = guide_legend(override.aes = list(size = 8), ncol = 4))


###################################################################################
### Integration Data ###
###################################################################################
## Integration
var0 = read.table(file = paste0(PATH, "/var.tsv"), sep = "\t", header = F, stringsAsFactors = F)
var0 = var0$V1
data0 = list(dat2, dat1)
InPlot(data0, var.gene = var0, ncore = 6)
data0 = scMultiIntegrate(data0, eigens = 40, var.gene = var0, add.Id = c("P11", "P2"), adjust = F, ncore = 1)
data0 = scTSNE(data0, npc = 42, use = "PLS")
ann0 = read.table(file = paste0(PATH, "/GSE110823_ann.tsv"), sep = "\t", header = T, stringsAsFactors = F)
ann0$scBarcode = paste0(ann0$Set, "_", ann0$Barcode)
cell0 = intersect(rownames(data0@coldata), ann0$scBarcode)
data0 = SubSet(data0, cells = cell0)
data0@coldata$CellType = factor(ann0$CellType, levels = paste0("C", 1:59))

color0 = c("#1B9E77", "#288E96", "#367EB6", "#419486", "#4BAC50", "#579055", "#636C63", "#667F49", "#669E26", "#6B9254", "#72789C", "#8064AD", "#9154A5", "#98649F", "#98899B", "#9C867A", "#A26643", "#A65D25", "#A66D20", "#B07117", "#C9650A", "#DB5206", "#E03013", "#E42F18", "#E5760B", "#E69B12", "#E65C54", "#E8318E", "#F05BA8", "#F780B3", "#FB7F56", "#FF8201", "#FFC01A", "#FFFF33", "#9E0142", "#AF1446", "#C0274A", "#D13A4E", "#DC494C", "#E65848", "#F06744", "#F57948", "#F88D51", "#FBA15B", "#FDB466", "#FDC373", "#FDD380", "#FEE18E", "#FEEB9E", "#FEF5AE", "#FFFFBF", "#F7FBB2", "#EFF8A6", "#E7F59A", "#D7EF9B", "#C4E79E", "#B2E0A2", "#9ED7A4", "#88CFA4", "#72C7A4", "#5FBAA8", "#4FA8AF", "#3F96B7", "#3484BB", "#4272B2", "#5060AA", "#5E4FA2")
DimPlot(data0, slot = "cell.tsne", colFactor = "Set")
DimPlot(data0, slot = "cell.tsne", colFactor = "CellType", Colors = color0)

DimPlot(data0, slot = "cell.tsne", genes = "Ebf3", size = 0.2)
DimPlot(data0, slot = "cell.tsne", genes = "Fat2", size = 0.2)


