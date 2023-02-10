library(RISC)
library(irlba)
library(umap)
library(ggplot2)
library(RColorBrewer)
library(biomaRt)


###################################################################################
### Prepare Data ###
###################################################################################
PATH = "/Path to the data/GSE84133"
Anno = read.table(file = paste0(PATH, "/GSE84133_Filter.tsv"), sep = "\t", header = T, stringsAsFactors = F)
ann0 = read.table(file = paste0(PATH, "/GSE84133_Anno.tsv"), sep = "\t", header = T, stringsAsFactors = F)

dat1 = read.table(file = paste0(PATH, '/Raw_Counts/GSM2230757_human1_umifm_counts.csv'), sep = ',', header = T, stringsAsFactors = F, check.names = F)
mat1 = as.matrix(t(dat1[,-c(1:3)]))
cel1 = dat1[,c(1, 3)]
colnames(cel1) = c('Barcode', 'Origi')
cel1$Set = 'Donor1'
colnames(mat1) = rownames(cel1) = cel1$Barcode
cel1 = cel1[cel1$Barcode %in% Anno$Barcode[Anno$Set == "Donor1"],]
mat1 = mat1[,colnames(mat1) %in% rownames(cel1)]

dat2 = read.table(file = paste0(PATH, '/Raw_Counts/GSM2230758_human2_umifm_counts.csv'), sep = ',', header = T, stringsAsFactors = F, check.names = F)
mat2 = as.matrix(t(dat2[,-c(1:3)]))
cel2 = dat2[,c(1, 3)]
colnames(cel2) = c('Barcode', 'Origi')
cel2$Set = 'Donor2'
colnames(mat2) = rownames(cel2) = cel2$Barcode
cel2 = cel2[cel2$Barcode %in% Anno$Barcode[Anno$Set == "Donor2"],]
mat2 = mat2[,colnames(mat2) %in% rownames(cel2)]

dat3 = read.table(file = paste0(PATH, '/Raw_Counts/GSM2230759_human3_umifm_counts.csv'), sep = ',', header = T, stringsAsFactors = F, check.names = F)
mat3 = as.matrix(t(dat3[,-c(1:3)]))
cel3 = dat3[,c(1, 3)]
colnames(cel3) = c('Barcode', 'Origi')
cel3$Set = 'Donor3'
colnames(mat3) = rownames(cel3) = cel3$Barcode
cel3 = cel3[cel3$Barcode %in% Anno$Barcode[Anno$Set == "Donor3"],]
mat3 = mat3[,colnames(mat3) %in% rownames(cel3)]

dat4 = read.table(file = paste0(PATH, '/Raw_Counts/GSM2230760_human4_umifm_counts.csv'), sep = ',', header = T, stringsAsFactors = F, check.names = F)
mat4 = as.matrix(t(dat4[,-c(1:3)]))
cel4 = dat4[,c(1, 3)]
colnames(cel4) = c('Barcode', 'Origi')
cel4$Set = 'Donor4'
colnames(mat4) = rownames(cel4) = cel4$Barcode
cel4 = cel4[cel4$Barcode %in% Anno$Barcode[Anno$Set == "Donor4"],]
mat4 = mat4[,colnames(mat4) %in% rownames(cel4)]

matrix1 = cbind(mat1, mat2, mat3, mat4)
keep = rowSums(matrix1 > 0) >= 0
matrix1 = matrix1[keep,]
gene0 = rownames(matrix1)
gen1 = gen2 = gen3 = gen4 = data.frame(Symbol = gene0, row.names = gene0, stringsAsFactors = F)

dat1 = readscdata(mat1[gene0,], cel1, gen1, is.filter = F)
dat2 = readscdata(mat2[gene0,], cel2, gen2, is.filter = F)
dat3 = readscdata(mat3[gene0,], cel3, gen3, is.filter = F)
dat4 = readscdata(mat4[gene0,], cel4, gen4, is.filter = F)

mouse = useMart('ensembl', dataset = 'mmusculus_gene_ensembl')
human = useMart('ensembl', dataset = 'hsapiens_gene_ensembl')
symbol = getLDS(attributes = 'hgnc_symbol', filters = 'hgnc_symbol', values = gene0, mart = human, attributesL = 'mgi_symbol', martL = mouse, uniqueRows = T)
colnames(symbol) = c('HGNC', 'MGI')
symbol0 = symbol[!duplicated(symbol$MGI, na.rm = T),]
symbol0 = symbol0[!duplicated(symbol0$HGNC, na.rm = T),]

dat5 = read.table(file = paste0(PATH, '/Raw_Counts/GSM2230761_mouse1_umifm_counts.csv'), sep = ',', header = T, stringsAsFactors = F, check.names = F)
mat5 = as.matrix(t(dat5[,-c(1:3)]))
cel5 = dat5[,c(1, 3)]
colnames(cel5) = c('Barcode', 'Origi')
cel5$Set = 'Mouse5'
colnames(mat5) = rownames(cel5) = cel5$Barcode
gene0 = intersect(symbol0$MGI, rownames(mat5))
symbol1 = symbol0[symbol0$MGI %in% gene0,]
mat5 = mat5[symbol1$MGI,]
rownames(mat5) = symbol1$HGNC
cel5 = cel5[cel5$Barcode %in% Anno$Barcode[Anno$Set == "Mouse1"],]
mat5 = mat5[,colnames(mat5) %in% rownames(cel5)]

dat6 = read.table(file = paste0(PATH, '/Raw_Counts/GSM2230762_mouse2_umifm_counts.csv'), sep = ',', header = T, stringsAsFactors = F, check.names = F)
mat6 = as.matrix(t(dat6[,-c(1:3)]))
cel6 = dat6[,c(1, 3)]
colnames(cel6) = c('Barcode', 'Origi')
cel6$Set = 'Mouse6'
colnames(mat6) = rownames(cel6) = cel6$Barcode
gene0 = intersect(symbol0$MGI, rownames(mat6))
symbol2 = symbol0[symbol0$MGI %in% gene0,]
mat6 = mat6[symbol2$MGI,]
rownames(mat6) = symbol2$HGNC
cel6 = cel6[cel6$Barcode %in% Anno$Barcode[Anno$Set == "Mouse2"],]
mat6 = mat6[,colnames(mat6) %in% rownames(cel6)]

matrix2 = cbind(mat5, mat6)
keep = rowSums(matrix2 > 0) >= 0
matrix2 = matrix2[keep,]
gene0 = rownames(matrix2)
gen5 = gen6 = data.frame(Symbol = gene0, row.names = gene0, stringsAsFactors = F)

dat5 = readscdata(mat5[gene0,], cel5, gen5, is.filter = F)
dat6 = readscdata(mat6[gene0,], cel6, gen6, is.filter = F)


###################################################################################
### RISC Objects ###
###################################################################################
process0 <- function(obj0){
  obj0 = scFilter(obj0, min.UMI = 1000, max.UMI = 20000, min.gene = 200, min.cell = 0, is.filter = F)
  obj0 = scNormalize(obj0)
  obj0 = scDisperse(obj0)
  print(length(obj0@vargene))
  return(obj0)
}

dat1 = process0(dat1)
dat2 = process0(dat2)
dat3 = process0(dat3)
dat4 = process0(dat4)
dat5 = process0(dat5)
dat6 = process0(dat6)


###################################################################################
### Uncorrect Data ###
###################################################################################
Color0 = data.frame(
  Type = c("activated_stellate", "acinar", "alpha", "beta", "delta", "ductal", "epsilon", "immune_other", "b_cell", "mast", "schwann", "t_cell", "gamma", "macrophage", "endothelial", "quiescent_stellate"), 
  Color = colorRampPalette(brewer.pal(11, 'Spectral'))(16), stringsAsFactors = F
)

coldata0 = rbind.data.frame(dat1@coldata, dat2@coldata, dat3@coldata, dat4@coldata)
var0 = Reduce(intersect, list(rownames(dat1@assay$logcount), rownames(dat2@assay$logcount), rownames(dat3@assay$logcount), rownames(dat4@assay$logcount)))
logmat1 = cbind(as.matrix(dat1@assay$logcount)[var0,], as.matrix(dat2@assay$logcount)[var0,], as.matrix(dat3@assay$logcount)[var0,], as.matrix(dat4@assay$logcount)[var0,])
pca0 = irlba(logmat1, nv = 10, center = T)$v
umap0 = umap(pca0)$layout
m0 = data.frame(UMAP1 = umap0[,1], UMAP2 = umap0[,2], Set = coldata0$Set, Type = coldata0$Origi)
m0$Type = factor(m0$Type, levels = c("activated_stellate", "acinar", "alpha", "beta", "delta", "ductal", "epsilon", "mast", "schwann", "t_cell", "gamma", "macrophage", "endothelial", "quiescent_stellate"))
color1 = Color0$Color[Color0$Type %in% as.character(m0$Type)]

ggplot(m0, aes(UMAP1, UMAP2)) + 
  geom_point(aes(color = Set, shape = Set), size = 1, alpha = 1) + 
  scale_color_manual(values = brewer.pal(4, 'Dark2')) + 
  scale_shape_manual(values = c(1, 2, 6, 8)) + 
  theme_bw(base_line_size = 0) + 
  labs(x = 'UMAP-1', y = 'UMAP-2', color = 'Set', shape = 'Source') + 
  theme(plot.title = element_text(hjust = 0.5, size = 16)) + 
  guides(colour = guide_legend(override.aes = list(size = 5), ncol = 1), shape = guide_legend(override.aes = list(size = 5), ncol = 1))
ggplot(m0, aes(UMAP1, UMAP2)) + 
  geom_point(aes(color = Type, shape = Set), size = 1, alpha = 1) + 
  scale_color_manual(values = color1) + 
  scale_shape_manual(values = c(1, 2, 6, 8)) + 
  theme_bw(base_line_size = 0) + 
  labs(x = 'UMAP-1', y = 'UMAP-2', color = 'Cell Type', shape = 'Source') + 
  theme(plot.title = element_text(hjust = 0.5, size = 16)) + 
  guides(colour = guide_legend(override.aes = list(size = 5), ncol = 1), shape = guide_legend(override.aes = list(size = 5), ncol = 1))

coldata0 = rbind.data.frame(dat5@coldata, dat6@coldata)
var0 = Reduce(intersect, list(rownames(dat5@assay$logcount), rownames(dat6@assay$logcount)))
logmat2 = cbind(as.matrix(dat5@assay$logcount)[var0,], as.matrix(dat6@assay$logcount)[var0,])
pca0 = irlba(logmat2, nv = 8, center = T)$v
umap0 = umap(pca0)$layout
m0 = data.frame(UMAP1 = umap0[,1], UMAP2 = umap0[,2], Set = coldata0$Set, Type = coldata0$Origi)
m0$Set = factor(m0$Set, levels = c("Mouse5", "Mouse6"), labels = c("Mouse1", "Mouse2"))
m0$Type = factor(
  m0$Type, 
  levels = c("activated_stellate", "beta", "ductal", "immune_other", "B_cell", "schwann", "T_cell", "gamma", "macrophage", "endothelial", "quiescent_stellate"), 
  labels = c("activated_stellate", "beta", "ductal", "immune_other", "b_cell", "schwann", "t_cell", "gamma", "macrophage", "endothelial", "quiescent_stellate")
)
color2 = Color0$Color[Color0$Type %in% m0$Type]

ggplot(m0, aes(UMAP1, UMAP2)) + 
  geom_point(aes(color = Set, shape = Set), size = 1, alpha = 1) + 
  scale_color_manual(values = brewer.pal(4, 'Dark2')) + 
  scale_shape_manual(values = c(4, 5)) + 
  theme_bw(base_line_size = 0) + 
  labs(x = 'UMAP-1', y = 'UMAP-2', color = 'Set', shape = 'Source') + 
  theme(plot.title = element_text(hjust = 0.5, size = 16)) + 
  guides(colour = guide_legend(override.aes = list(size = 5), ncol = 1), shape = guide_legend(override.aes = list(size = 5), ncol = 1))
ggplot(m0, aes(UMAP1, UMAP2)) + 
  geom_point(aes(color = Type, shape = Set), size = 1, alpha = 1) + 
  scale_color_manual(values = color2) + 
  scale_shape_manual(values = c(4, 5)) + 
  theme_bw(base_line_size = 0) + 
  labs(x = 'UMAP-1', y = 'UMAP-2', color = 'Cell Type', shape = 'Source') + 
  theme(plot.title = element_text(hjust = 0.5, size = 16)) + 
  guides(colour = guide_legend(override.aes = list(size = 5), ncol = 1), shape = guide_legend(override.aes = list(size = 5), ncol = 1))


###################################################################################
### RISC Integration ###
###################################################################################
## Human
data1 = list(dat1, dat2, dat3, dat4)
var0 = Reduce(intersect, list(dat1@vargene, dat2@vargene, dat3@vargene, dat4@vargene))
InPlot(data1, var.gene = var0)
data1 = scMultiIntegrate(data1, eigens = 14, var.gene = var0, ncore = 4, add.Id = c("Donor1", "Donor2", "Donor3", "Donor4"))
data1 = scUMAP(data1, npc = 14, use = 'PLS')
data1@coldata$Type = as.factor(ann0$Type[ann0$Set %in% c("Donor1", "Donor2", "Donor3", "Donor4")])

DimPlot(data1, slot = "cell.umap", colFactor = "Set", Colors = brewer.pal(4, "Dark2"))
DimPlot(data1, slot = "cell.umap", colFactor = "Type")

## Mouse
data2 = list(dat6, dat5)
var0 = Reduce(intersect, list(dat5@vargene, dat6@vargene))
InPlot(data2, var.gene = var0)
data2 = scMultiIntegrate(data2, eigens = 14, var.gene = var0, ncore = 4, add.Id = c("Mouse2", "Mouse1"))
data2 = scUMAP(data2, npc = 14, use = 'PLS')
data2@coldata$Type = as.factor(ann0$Type[ann0$Set %in% c("Mouse2", "Mouse1")])

DimPlot(data2, slot = "cell.umap", colFactor = "Set", Colors = brewer.pal(4, "Dark2"))
DimPlot(data2, slot = "cell.umap", colFactor = "Type")

## Both species
data3 = list(dat1, dat2, dat3, dat4, dat5, dat6)
var0 = read.table(file = paste0(PATH, "/var.tsv"), sep = "\t", header = F, stringsAsFactors = F)
var0 = var0$V1
InPlot(data3, var.gene = var0)
data3 = scMultiIntegrate(data3, eigens = 15, var.gene = var0, ncore = 4, add.Id = c("Donor1", "Donor2", "Donor3", "Donor4", "Mouse1", "Mouse2"))
data3 = scUMAP(data3, npc = 15, use = 'PLS')
data3@coldata$Type = as.factor(ann0$Type)

DimPlot(data3, slot = "cell.umap", colFactor = "Set", Colors = brewer.pal(6, "Dark2"))
DimPlot(data3, slot = "cell.umap", colFactor = "Type")

DimPlot(data3, slot = "cell.umap", genes = "GCG", size = 0.5)
DimPlot(data3, slot = "cell.umap", genes = "INS", size = 0.5, exp.range = c(6, 10))
DimPlot(data3, slot = "cell.umap", genes = "RESP18", size = 0.5)
DimPlot(data3, slot = "cell.umap", genes = "KRT19", size = 0.5)
DimPlot(data3, slot = "cell.umap", genes = "COL6A3", size = 0.5)
DimPlot(data3, slot = "cell.umap", genes = "SST", size = 0.5, exp.range = c(3, 11))
DimPlot(data3, slot = "cell.umap", genes = "RGS5", size = 0.5)
DimPlot(data3, slot = "cell.umap", genes = "PRSS1", size = 0.5)
DimPlot(data3, slot = "cell.umap", genes = "FLT1", size = 0.5)


