library(RISC)
library(irlba)
library(umap)
library(ggplot2)
library(RColorBrewer)


###################################################################################
### Prepare Data ###
###################################################################################
PATH = "/Path to the data/GSE85241_GSE81076_GSE83139_EMTAB_5061"

## GSE85241 CEL-Seq2
gse85241.df <- read.table(file = paste0(PATH, "/Raw_Counts/GSE85241_cellsystems_dataset_4donors_updated.csv"), sep = '\t', h = TRUE, row.names = 1, stringsAsFactors = FALSE)
dim(gse85241.df)
donor.names <- sub("^(D[0-9]+).*", "\\1", colnames(gse85241.df))
table(donor.names)
plate.id <- sub("^D[0-9]+\\.([0-9]+)_.*", "\\1", colnames(gse85241.df))
table(plate.id)
gene.symb <- gsub("__chr.*$", "", rownames(gse85241.df))
is.spike <- grepl("^ERCC-", gene.symb)
table(is.spike)
gse85241.df$Symbol = gene.symb
gse85241.df = gse85241.df[!duplicated(gse85241.df$Symbol),]
rownames(gse85241.df) = gse85241.df$Symbol
gse85241.df = gse85241.df[,-3073]

count0 = gse85241.df
cell = data.frame(Donor = donor.names, Plate = plate.id)
rownames(cell) = colnames(count0)
gene = data.frame(Symbol = rownames(count0))
rownames(gene) = rownames(count0)

data1 = readscdata(count = count0, cell = cell, gene = gene)
rm(cell, count0, donor.names, gene, gene.symb, gse85241.df, is.spike, plate.id)


## GSE81076 CEL-Seq
gse81076.df <- read.table(file = paste0(PATH, "/Raw_Counts/GSE81076_D2_3_7_10_17.txt"), sep='\t', header=TRUE, stringsAsFactors=FALSE, row.names=1)
donor.names <- sub("^(D[0-9]+).*", "\\1", colnames(gse81076.df))
table(donor.names)
plate.id <- sub("^D[0-9]+(.*)_.*", "\\1", colnames(gse81076.df))
table(plate.id)
gene.symb <- gsub("__chr.*$", "", rownames(gse81076.df))
is.spike <- grepl("^ERCC-", gene.symb)
table(is.spike)
gse81076.df$Symbol = gene.symb
gse81076.df = gse81076.df[!duplicated(gse81076.df$Symbol),]
rownames(gse81076.df) = gse81076.df$Symbol
gse81076.df = gse81076.df[,-1729]

count0 = gse81076.df
cell = data.frame(Donor = donor.names, Plate = plate.id)
rownames(cell) = colnames(count0)
gene = data.frame(Symbol = rownames(count0))
rownames(gene) = rownames(count0)

data2 = readscdata(count = count0, cell = cell, gene = gene)
rm(cell, count0, donor.names, gene, gene.symb, gse81076.df, is.spike, plate.id)


## E-MTAB-5061 Smart-Seq2
header <- read.table(file = paste0(PATH, "/Raw_Counts/pancreas_refseq_rpkms_counts_3514sc.txt"), nrow = 1, sep = "\t", comment.char = "", stringsAsFactors = FALSE)
ncells <- ncol(header) - 1L
col.types <- vector("list", ncells*2 + 2)
col.types[1:2] <- "character"
col.types[2+ncells + seq_len(ncells)] <- "integer"
e5601.df <- read.table(file = paste0(PATH, "/Raw_Counts/pancreas_refseq_rpkms_counts_3514sc.txt"), sep = "\t", colClasses = col.types)
gene.data <- e5601.df[,1:2]
e5601.df <- e5601.df[,-(1:2)]
colnames(e5601.df) <- as.character(header[1,-1])
dim(e5601.df)
is.spike <- grepl("^ERCC-", gene.data[,2])
table(is.spike)
e5601.df = data.frame(e5601.df, Symbol = gene.data$V1)
e5601.df = e5601.df[!duplicated(e5601.df$Symbol),]
rownames(e5601.df) = e5601.df$Symbol
e5601.df = e5601.df[,-3515]
col.types[[1]] = sapply(colnames(e5601.df), function(x){strsplit(x, '_', fixed = T)[[1]][1]})
col.types[[2]] = sapply(colnames(e5601.df), function(x){strsplit(x, '_', fixed = T)[[1]][2]})

count0 = e5601.df
cell = data.frame(Donor = col.types[[1]], Plate = col.types[[2]])
rownames(cell) = colnames(count0)
gene = data.frame(Symbol = rownames(count0))
rownames(gene) = rownames(count0)

data3 = readscdata(count = count0, cell = cell, gene = gene)
rm(cell, count0, gene, e5601.df, is.spike, col.types, gene.data, header, ncells)


## GSE83139 SMARTer
gse83139.df <- read.table(file = paste0(PATH, "/Raw_Counts/GSE83139_tbx-v-f-norm-ntv-cpms.csv"), sep = "\t", header = T, stringsAsFactors = F)
gse83139.df <- gse83139.df[!duplicated(gse83139.df$gene),]
gse83139.df <- gse83139.df[!is.na(gse83139.df$gene),]
count0 = as.matrix(gse83139.df[,-c(1:7)])
rownames(count0) = as.character(gse83139.df$gene)
cell = data.frame(Sample = colnames(count0), row.names = colnames(count0), stringsAsFactors = F)
gene = data.frame(ID = rownames(count0), row.names = rownames(count0), stringsAsFactors = F)

data4 = readscdata(count = count0, cell = cell, gene = gene)
rm(count0, gse83139.df, cell, gene)


###################################################################################
### RISC Data ###
###################################################################################
process1 <- function(obj0){
  obj0 <- scFilter(obj0, min.UMI = 3000, max.UMI = 100000, min.gene = 1000, min.cell = 5)
  obj0 <- scNormalize(obj0)
  obj0 <- scDisperse(obj0)
  print(length(obj0@vargene))
  obj0 <- scPCA(obj0)
  return(obj0)
}

process2 <- function(obj0){
  obj0 <- scFilter(obj0)
  obj0 <- scNormalize(obj0)
  obj0 <- scDisperse(obj0)
  print(length(obj0@vargene))
  obj0 <- scPCA(obj0)
  return(obj0)
}

data1 = process1(data1)
data2 = process1(data2)
data3 = process1(data3)
data4 = process2(data4)


###################################################################################
### Uncorrect Data ###
###################################################################################
var0 = read.table(file = paste0(PATH, "/var.tsv"), sep = "\t", header = F, stringsAsFactors = F)
var0 = var0$V1
logmat0 = cbind(
  as.matrix(data1@assay$logcount)[var0,], 
  as.matrix(data2@assay$logcount)[var0,], 
  as.matrix(data3@assay$logcount)[var0,], 
  as.matrix(data4@assay$logcount)[var0,]
)
pca0 = irlba(logmat0, nv = 12, center = T)$v
umap0 = umap(pca0)$layout
ann0 = read.table(file = paste0(PATH, "/Pancreas_Annotation.tsv"), sep = "\t", header = T, stringsAsFactors = F)
m0 = data.frame(UMAP1 = umap0[,1], UMAP2 = umap0[,2], Set = ann0$Set, Type = ann0$Type)
m0$Set = factor(m0$Set, levels = c("E-MTAB-5061", "GSE81076", "GSE85241", "GSE83139"))

ggplot(m0, aes(UMAP1, UMAP2)) + 
  geom_point(aes(color = Set, shape = Set), size = 2, alpha = 1) + 
  scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A")) + 
  scale_shape_manual(values = c(1, 2, 6, 8)) + 
  theme_bw(base_line_size = 0) + 
  labs(x = 'UMAP-1', y = 'UMAP-2', color = 'Cell Type', shape = 'Set') + 
  theme(plot.title = element_text(hjust = 0.5, size = 16)) + 
  guides(colour = guide_legend(override.aes = list(size = 5), ncol = 1), shape = guide_legend(override.aes = list(size = 5), ncol = 1))
ggplot(m0, aes(UMAP1, UMAP2)) + 
  geom_point(aes(color = Type, shape = Set), size = 2, alpha = 1) + 
  scale_color_manual(values = brewer.pal(8, 'Spectral')) + 
  scale_shape_manual(values = c(1, 2, 6, 8)) + 
  theme_bw(base_line_size = 0) + 
  labs(x = 'UMAP-1', y = 'UMAP-2', color = 'Cell Type', shape = 'Set') + 
  theme(plot.title = element_text(hjust = 0.5, size = 16)) + 
  guides(colour = guide_legend(override.aes = list(size = 5), ncol = 1), shape = guide_legend(override.aes = list(size = 5), ncol = 1))


###################################################################################
### Data Integration ###
###################################################################################
# Integrating
data0 = list(data1, data2, data3, data4)
var0 = read.table(file = paste0(PATH, "/var.tsv"), sep = "\t", header = F, stringsAsFactors = F)
var0 = var0$V1
InPlot(data0, var.gene = var0, minPC = 8, nPC = 12)
data0 = scMultiIntegrate(data0, eigens = 8, var.gene = var0, ncore = 4, adjust = F, add.Id = c("GSE85241", "GSE81076", "E-MTAB-5061", "GSE83139"))
data0 = scUMAP(data0, npc = 11, use = 'PLS')
ann0 = read.table(file = paste0(PATH, "/Pancreas_Annotation.tsv"), sep = "\t", header = T, stringsAsFactors = F)
data0@coldata$Platform = as.character(ann0$Platform)
data0@coldata$Type = as.character(ann0$Type)

DimPlot(data0, slot = "cell.umap", colFactor = "Set", Colors = brewer.pal(4, "Dark2"))
DimPlot(data0, slot = "cell.umap", colFactor = "Platform", Colors = brewer.pal(4, "Dark2"))
DimPlot(data0, slot = "cell.umap", colFactor = "Type", label = T)
DimPlot(data0, slot = "cell.umap", genes = "SST", size = 1)
DimPlot(data0, slot = "cell.umap", genes = "PPY", size = 1)
DimPlot(data0, slot = "cell.umap", genes = "GCG", size = 1)
DimPlot(data0, slot = "cell.umap", genes = "INS", size = 1)
DimPlot(data0, slot = "cell.umap", genes = "COL1A1", size = 1)
DimPlot(data0, slot = "cell.umap", genes = "PRSS1", size = 1)
DimPlot(data0, slot = "cell.umap", genes = "KRT19", size = 1)
DimPlot(data0, slot = "cell.umap", genes = "FLT1", size = 1)


