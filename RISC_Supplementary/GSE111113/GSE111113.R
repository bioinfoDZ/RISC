library(RISC)
library(ggplot2)
library(RColorBrewer)
library(irlba)
library(umap)


PATH = "/Path to the data/GSE111113"

###################################################################################
### RISC raw ###
###################################################################################
mat0 = read.table(file = paste0(PATH, '/GSE111113_Table_S1_FilterNormal10xExpMatrix.txt'), sep = '\t', header = T, stringsAsFactors = F)
mat1 = as.matrix(mat0[,-c(1:3)])
rownames(mat1) = mat0$gene_id
Group = sapply(colnames(mat1), function(x){strsplit(x, "_")[[1]][1]})
Symbol0 = mat0[,c(1, 3)]
colnames(Symbol0) = c('Ensembl', 'Symbol')


###################################################################################
### Uncorrected data ###
###################################################################################
dat0 = readscdata(count = mat1, cell = data.frame(Time = Group, row.names = colnames(mat1)), gene = data.frame(Symbol = rownames(mat1), row.names = rownames(mat1)))
dat0 = scFilter(dat0, min.UMI = 500, max.UMI = Inf, min.gene = 200, min.cell = 5)
cell0 = rownames(dat0@coldata)[dat0@coldata$Time %in% c('E16', 'E18', 'P4', 'Adu1', 'Adu2')]
dat0 = SubSet(dat0, cells = cell0)
dat0 = scNormalize(dat0)
dat0 = scDisperse(dat0)
length(dat0@vargene)
dat0 = scPCA(dat0)
dat0 = scUMAP(dat0)

UMAPlot(dat0, colFactor = 'Time', Colors = brewer.pal(5, 'Spectral'))


###################################################################################
### RISC integration ###
###################################################################################
cell0 = rownames(dat0@coldata)[dat0@coldata$Time == 'E16']
dat1 = SubSet(dat0, cells = cell0)
cell0 = rownames(dat0@coldata)[dat0@coldata$Time == 'E18']
dat2 = SubSet(dat0, cells = cell0)
cell0 = rownames(dat0@coldata)[dat0@coldata$Time == 'P4']
dat3 = SubSet(dat0, cells = cell0)
cell0 = rownames(dat0@coldata)[dat0@coldata$Time == 'Adu1']
dat4 = SubSet(dat0, cells = cell0)
cell0 = rownames(dat0@coldata)[dat0@coldata$Time == 'Adu2']
dat5 = SubSet(dat0, cells = cell0)

dat1 = scDisperse(dat1)
length(dat1@vargene)
dat2= scDisperse(dat2)
length(dat2@vargene)
dat3= scDisperse(dat3)
length(dat3@vargene)
dat4= scDisperse(dat4)
length(dat4@vargene)
dat5= scDisperse(dat5)
length(dat5@vargene)


###################################################################################
### RISC integration ###
###################################################################################
### Integration All
var0 = read.table(file = paste0(PATH, "/var.tsv"), sep = "\t", header = F, stringsAsFactors = F)
var0 = var0$V1
dat.all = list(dat4, dat5, dat3, dat2, dat1)
InPlot(dat.all, var.gene = var0)
dat.all = scMultiIntegrate(dat.all, eigens = 15, var.gene = var0, add.Id = c('Adult1', 'Adult2', 'P4', 'E18', 'E16'), ncore = 4)
dat.all = scUMAP(dat.all, npc = 15, use = 'PLS')
dat.all@coldata$Set = factor(dat.all@coldata$Set, levels = c('E16', 'E18', 'P4', 'Adult1', 'Adult2'))

DimPlot(dat.all, slot = "cell.umap", colFactor = 'Set')
UMAPlot(dat.all, genes = Symbol0$Ensembl[Symbol0$Symbol == 'Wfdc18'])
UMAPlot(dat.all, genes = Symbol0$Ensembl[Symbol0$Symbol == 'Sostdc1'])

## Integration Patial
dat.par = list(dat4, dat5, dat2, dat1)
dat.par = scMultiIntegrate(dat.par, eigens = 20, var.gene = var0, add.Id = c('Adult1', 'Adult2', 'E18', 'E16'), ncore = 4)
dat.par = scUMAP(dat.par, npc = 20, use = 'PLS')
dat.par@coldata$Set = factor(dat.par@coldata$Set, levels = c('E16', 'E18', 'Adult1', 'Adult2'))

DimPlot(dat.par, slot = "cell.umap", colFactor = 'Set')
UMAPlot(dat.par, genes = Symbol0$Ensembl[Symbol0$Symbol == 'Wfdc18'])
UMAPlot(dat.par, genes = Symbol0$Ensembl[Symbol0$Symbol == 'Sostdc1'])


