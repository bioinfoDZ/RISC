library(RISC)
library(RColorBrewer)
library(ggplot2)


colname <- function(x0){
  colname0 <- colnames(x0)[-ncol(x0)]
  x0 <- x0[,-ncol(x0)]
  colnames(x0) <- colname0
  x0 <- as.matrix(x0)
  rownames(x0) <- sapply(rownames(x0), function(x){strsplit(x, "__", fixed = T)[[1]][1]})
  return(x0)
}

PATH = "/Path to the data/GSE125688"

###################################################################################
### Prepare data ###
###################################################################################
Ctrl1 <- read.table(file = paste0(PATH, '/Raw_Counts/GSM3580367_bil-adult1.coutb.csv'), sep = '\t', header = T, stringsAsFactors = F)
Ctrl2 <- read.table(file = paste0(PATH, '/Raw_Counts/GSM3580368_bil-adult2.coutb.csv'), sep = '\t', header = T, stringsAsFactors = F)
Ctrl3 <- read.table(file = paste0(PATH, '/Raw_Counts/GSM3580369_bil-adult3.coutb.csv'), sep = '\t', header = T, stringsAsFactors = F)
DDC <- read.table(file = paste0(PATH, '/Raw_Counts/GSM3580370_bil-DDC1.coutb.csv'), sep = '\t', header = T, stringsAsFactors = F)

Ctrl1 <- colname(Ctrl1)
Ctrl2 <- colname(Ctrl2)
Ctrl3 <- colname(Ctrl3)
DDC <- colname(DDC)

Ctrl.dat1 <- readscdata(count = Ctrl1, cell = data.frame(Barcode = colnames(Ctrl1), row.names = colnames(Ctrl1)), gene = data.frame(Symbol = rownames(Ctrl1), row.names = rownames(Ctrl1)))
Ctrl.dat2 <- readscdata(count = Ctrl2, cell = data.frame(Barcode = colnames(Ctrl2), row.names = colnames(Ctrl2)), gene = data.frame(Symbol = rownames(Ctrl2), row.names = rownames(Ctrl2)))
Ctrl.dat3 <- readscdata(count = Ctrl3, cell = data.frame(Barcode = colnames(Ctrl3), row.names = colnames(Ctrl3)), gene = data.frame(Symbol = rownames(Ctrl3), row.names = rownames(Ctrl3)))
DDC.dat <- readscdata(count = DDC, cell = data.frame(Barcode = colnames(DDC), row.names = colnames(DDC)), gene = data.frame(Symbol = rownames(DDC), row.names = rownames(DDC)))


###################################################################################
### Process data ###
###################################################################################
Ctrl.dat1 = scFilter(Ctrl.dat1, min.UMI = 500, max.UMI = 8000, min.gene = 500, min.cell = 5)
Ctrl.dat2 = scFilter(Ctrl.dat2, min.UMI = 500, max.UMI = 8000, min.gene = 500, min.cell = 5)
Ctrl.dat3 = scFilter(Ctrl.dat3, min.UMI = 500, max.UMI = 8000, min.gene = 500, min.cell = 5)
DDC.dat = scFilter(DDC.dat, min.UMI = 500, max.UMI = 8000, min.gene = 500, min.cell = 5)

Ctrl.dat1 = scNormalize(Ctrl.dat1)
Ctrl.dat2 = scNormalize(Ctrl.dat2)
Ctrl.dat3 = scNormalize(Ctrl.dat3)
DDC.dat = scNormalize(DDC.dat)

Ctrl.dat1 = scDisperse(Ctrl.dat1)
Ctrl.dat2 = scDisperse(Ctrl.dat2)
Ctrl.dat3 = scDisperse(Ctrl.dat3)
DDC.dat = scDisperse(DDC.dat)

length(Ctrl.dat1@vargene)
length(Ctrl.dat2@vargene)
length(Ctrl.dat3@vargene)
length(DDC.dat@vargene)

Ctrl.dat1 = scPCA(Ctrl.dat1)
Ctrl.dat2 = scPCA(Ctrl.dat2)
Ctrl.dat3 = scPCA(Ctrl.dat3)
DDC.dat = scPCA(DDC.dat)


###################################################################################
### Uncorrected data ###
###################################################################################
## Uncorrect
library(irlba)
library(umap)

gene0 = Reduce(intersect, list(Ctrl.dat1@rowdata$Symbol, Ctrl.dat2@rowdata$Symbol, Ctrl.dat3@rowdata$Symbol, DDC.dat@rowdata$Symbol))
logmat0 = cbind(Ctrl.dat1@assay$logcount[gene0,], Ctrl.dat2@assay$logcount[gene0,], Ctrl.dat3@assay$logcount[gene0,], DDC.dat@assay$logcount[gene0,])
pca0 = irlba(logmat0, nv = 50)$v
umap0 = umap(pca0[,1:9])$layout
m0 = data.frame(UMAP1 = umap0[,1], UMAP2 = umap0[,2], Set = c(rep("BEC1", nrow(Ctrl.dat1@coldata)), rep("BEC2", nrow(Ctrl.dat2@coldata)), rep("BEC3", nrow(Ctrl.dat3@coldata)), rep("DDC", nrow(DDC.dat@coldata))))
ggplot(m0, aes(UMAP1, UMAP2)) + 
  geom_point(aes(color = Set), shape = 19, size = 1, alpha = 1) + 
  scale_color_manual(values = c("#9ECAE1", "#4292C6", "#08519C", "#EF3B2C")) + 
  theme_bw(base_line_size = 0, base_size = 12) + 
  labs(x = 'UMAP1', y = 'UMAP2', color = 'Set') + 
  theme(plot.title = element_text(hjust = 0.5, size = 16)) + 
  guides(colour = guide_legend(override.aes = list(size = 5), ncol = 1))


###################################################################################
### Integration data ###
###################################################################################
## Integration
BEC = list(Ctrl.dat3, Ctrl.dat2, Ctrl.dat1, DDC.dat)
var0 = Reduce(intersect, list(
  rownames(Ctrl.dat1@rowdata), rownames(Ctrl.dat2@rowdata), 
  rownames(Ctrl.dat3@rowdata), rownames(DDC.dat@rowdata)
))
BEC = scMultiIntegrate(BEC, eigens = 12, var.gene = var0, ncore = 4, add.Id = c('Ctrl3', 'Ctrl2', 'Ctrl1', 'DDC'))
BEC = scUMAP(BEC, npc = 12, use = 'PLS')

DimPlot(BEC, slot = "cell.umap", colFactor = "Set", size = 2, Colors = c("#9ECAE1", "#4292C6", "#08519C", "#EF3B2C"))


