library(RISC)
library(irlba)
library(umap)
library(ggplot2)
library(RColorBrewer)


PATH = "/Path to the data/GSE131181"

###################################################################################
### Prepare data ###
###################################################################################
library(data.table)

E10 = fread(file = paste0(PATH, '/Raw_Counts/GSE131181_e10.5.raw.data.csv'), sep = ',')
E13 = fread(file = paste0(PATH, '/Raw_Counts/GSE131181_e13.5.raw.data.csv'), sep = ',')
dE10 = E10[,-1]
gene1 = E10$V1
dE13 = E13[,-1]
gene2 = E13$V1
rm(E10, E13)

mE10 = fread(file = paste0(PATH, '/Raw_Counts/GSE131181_e10.5.meta.data.csv'), sep = ',')
mE13 = fread(file = paste0(PATH, '/Raw_Counts/GSE131181_e13.5.meta.data.csv'), sep = ',')
mdE10 = mE10[,c(1, 4, 9, 10)]
mdE13 = mE13[,c(1, 4, 10, 8)]
colnames(mdE10) = colnames(mdE13) = c('Sample', 'Experiment', 'Group', 'Set')
rm(mE10, mE13)

dE10 = dE10[, mdE10$Sample, with = F]
dE13 = dE13[, mdE13$Sample, with = F]

mdE10 = as.data.frame(mdE10)
mE10wt = mdE10[mdE10$Set == 'Control', , drop = F]
mE10pd = mdE10[mdE10$Set == 'PhD', , drop = F]
mdE13 = as.data.frame(mdE13)
mE13wt = mdE13[mdE13$Set == 'Control', , drop = F]
mE13pd = mdE13[mdE13$Set == 'PhD', , drop = F]

E10wt = as.matrix(dE10[, mE10wt$Sample, with = F])
E10pd = as.matrix(dE10[, mE10pd$Sample, with = F])
E13wt = as.matrix(dE13[, mE13wt$Sample, with = F])
E13pd = as.matrix(dE13[, mE13pd$Sample, with = F])
rm(dE10, dE13, mdE10, mdE13)
rownames(E10wt) = rownames(E10pd) = gene1
rownames(E13wt) = rownames(E13pd) = gene2
rownames(mE10wt) = mE10wt$Sample
rownames(mE10pd) = mE10pd$Sample
rownames(mE13wt) = mE13wt$Sample
rownames(mE13pd) = mE13pd$Sample


###################################################################################
### Individual datasets ###
###################################################################################
dat1 = readscdata(count = E10wt, cell = mE10wt, gene = data.frame(Symbol = gene1, row.names = gene1))
dat2 = readscdata(count = E10pd, cell = mE10pd, gene = data.frame(Symbol = gene1, row.names = gene1))
dat3 = readscdata(count = E13wt, cell = mE13wt, gene = data.frame(Symbol = gene2, row.names = gene2))
dat4 = readscdata(count = E13pd, cell = mE13pd, gene = data.frame(Symbol = gene2, row.names = gene2))
rm(E10wt, E10pd, E13wt, E13pd, mE10wt, mE10pd, mE13wt, mE13pd, gene1, gene2)

process0 <- function(obj0){
  obj0 = scFilter(obj0, min.UMI = 1000, max.UMI = Inf, min.gene = 200, min.cell = 5)
  obj0 = scNormalize(obj0)
  obj0 = scDisperse(obj0)
  length(obj0@vargene)
  return(obj0)
}

dat1 = process0(dat1)
dat2 = process0(dat2)
dat3 = process0(dat3)
dat4 = process0(dat4)

FilterPlot(dat1)
FilterPlot(dat2)
FilterPlot(dat3)
FilterPlot(dat4)


###################################################################################
### Uncorrect Data ###
###################################################################################
# Integration
ann0 = read.table(file = paste0(PATH, "/GSE131181_Ann.tsv"), sep = "\t", header = T, stringsAsFactors = F)
var0 = read.table(file = paste0(PATH, "/GSE131181_var.tsv"), sep = "\t", header = F, stringsAsFactors = F)
var0 = var0$V1
logmat0 = cbind(dat3@assay$logcount[var0,], dat4@assay$logcount[var0,], dat1@assay$logcount[var0,], dat2@assay$logcount[var0,])
pca0 = irlba(as.matrix(logmat0), nv = 50)$v
umap0 = umap(pca0[,1:25])$layout

m0 = data.frame(UMAP1 = umap0[,1], UMAP2 = umap0[,2], Set = ann0$Set, Type = ann0$Type)
ggplot(m0, aes(UMAP1, UMAP2)) + 
  geom_point(aes(color = Set, shape = Set), size = 0.1, alpha = 0.85) + 
  scale_color_manual(values = c('#FB8072', '#80B1D3', '#8DD3C7', '#BC80BD')) + 
  scale_shape_manual(values = c(3, 4, 2, 6)) + 
  theme_bw(base_line_size = 0, base_size = 24) + 
  labs(x = 'UMAP-1', y = 'UMAP-2', color = 'Set', shape = 'Set') + 
  theme(plot.title = element_text(hjust = 0.5, size = 16)) + 
  guides(colour = guide_legend(override.aes = list(size = 6), ncol = 1))

ggplot(m0, aes(UMAP1, UMAP2)) + 
  geom_point(aes(color = Type, shape = Set), size = 0.1, alpha = 0.9) + 
  scale_color_manual(values = colorRampPalette(brewer.pal(11, 'Spectral'))(18)) + 
  scale_shape_manual(values = c(3, 4, 2, 6)) + 
  theme_bw(base_line_size = 0, base_size = 24) + 
  labs(x = 'UMAP-1', y = 'UMAP-2', color = 'Cell Type', shape = 'Set') + 
  theme(plot.title = element_text(hjust = 0.5, size = 16)) + 
  guides(colour = guide_legend(override.aes = list(size = 6), ncol = 1), shape = guide_legend(override.aes = list(size = 5), ncol = 1))


###################################################################################
### Data Integration ###
###################################################################################
# Integration
var0 = read.table(file = paste0(PATH, "/GSE131181_var.tsv"), sep = "\t", header = F, stringsAsFactors = F)
var0 = var0$V1
dat.all = list(dat3, dat4, dat1, dat2)
InPlot(dat.all, var.gene = var0, nPC = 30, ncore = 4)
dat.all = scMultiIntegrate(dat.all, eigens = 20, var.gene = var0, ncore = 2, adjust = FALSE, add.Id = c('E13wt', 'E13pd', 'E10wt', 'E10pd'))
dat.all = scUMAP(dat.all, npc = 25, use = 'PLS')
dat.all = scTSNE(dat.all, npc = 25, use = 'PLS')

ann0 = read.table(file = paste0(PATH, "/GSE131181_Ann.tsv"), sep = "\t", header = T, stringsAsFactors = F)
dat.all@coldata$Type = ann0$Type

DimPlot(dat.all, slot = "cell.umap", colFactor = "Set", size = 0.2)
DimPlot(dat.all, slot = "cell.umap", colFactor = "Type", size = 0.2)


