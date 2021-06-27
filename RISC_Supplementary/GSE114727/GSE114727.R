library(RISC)
library(irlba)
library(umap)
library(ggplot2)
library(RColorBrewer)


PATH = "/Path to the data/GSE114727"

###################################################################################
### Raw Data ###
###################################################################################
# BC09 rep1
cell1 = read.csv(file = paste0(PATH, "/10X_Genomics/GSM3148580_BC09_TUMOR1_filtered_contig_annotations.csv"), sep = ",", header = T, stringsAsFactors = F)
cell1 = cell1[cell1$full_length == 'True' & cell1$productive == 'True',]
data1 = read10Xgenomics(data.path = paste0(PATH, "/10X_Genomics/BC09_Tech1"))
data1 = SubSet(data1, cells = cell1$barcode)

# BC09 rep2
cell2 = read.csv(file = paste0(PATH, "/10X_Genomics/GSM3148581_BC09_TUMOR2_filtered_contig_annotations.csv"), sep = ",", header = T, stringsAsFactors = F)
cell2 = cell2[cell2$full_length == 'True' & cell2$productive == 'True',]
data2 = read10Xgenomics(data.path = paste0(PATH, "/10X_Genomics/BC09_Tech2"))
data2 = SubSet(data2, cells = cell2$barcode)

# BC10
cell3 = read.csv(file = paste0(PATH, "/10X_Genomics/GSM3148582_BC10_TUMOR1_filtered_contig_annotations.csv"), sep = ",", header = T, stringsAsFactors = F)
cell3 = cell3[cell3$full_length == 'True' & cell3$productive == 'True',]
data3 = read10Xgenomics(data.path = paste0(PATH, "/10X_Genomics/BC10_Tech1"))
data3 = SubSet(data3, cells = cell3$barcode)

# BC11 rep1
cell4 = read.csv(file = paste0(PATH, "/10X_Genomics/GSM3148583_BC11_TUMOR1_filtered_contig_annotations.csv"), sep = ",", header = T, stringsAsFactors = F)
cell4 = cell4[cell4$full_length == 'True' & cell4$productive == 'True',]
data4 = read10Xgenomics(data.path = paste0(PATH, "/10X_Genomics/BC11_Tech1"))
data4 = SubSet(data4, cells = cell4$barcode)

# BC11 rep2
cell5 = read.csv(file = paste0(PATH, "/10X_Genomics/GSM3148584_BC11_TUMOR2_filtered_contig_annotations.csv"), sep = ",", header = T, stringsAsFactors = F)
cell5 = cell5[cell5$full_length == 'True' & cell5$productive == 'True',]
data5 = read10Xgenomics(data.path = paste0(PATH, "/10X_Genomics/BC11_Tech2"))
data5 = SubSet(data5, cells = cell5$barcode)

# GSE110686
data6 = read10Xgenomics(data.path = paste0(PATH, "/10X_Genomics/GSE110686"))


###################################################################################
### Prepare Data ###
###################################################################################
process0 <- function(obj0){
  obj0 = scFilter(obj0, min.UMI = 1000, max.UMI = 20000, min.gene = 500, min.cell = 5)
  obj0 = scNormalize(obj0, ncore = 4)
  obj0 = scDisperse(obj0)
  print(length(obj0@vargene))
  return(obj0)
}

data1 = process0(data1)
data2 = process0(data2)
data3 = process0(data3)
data4 = process0(data4)
data5 = process0(data5)
data6 = process0(data6)

FilterPlot(data1)
FilterPlot(data2)
FilterPlot(data3)
FilterPlot(data4)
FilterPlot(data5)
FilterPlot(data6)


###################################################################################
### Uncorrect Data ###
###################################################################################
var0 = Reduce(intersect, list(
  rownames(data1@assay$logcount), rownames(data2@assay$logcount), rownames(data3@assay$logcount), 
  rownames(data4@assay$logcount), rownames(data5@assay$logcount), rownames(data6@assay$logcount)
))
logmat0 = cbind(
  as.matrix(data2@assay$logcount)[var0,], as.matrix(data1@assay$logcount)[var0,], 
  as.matrix(data3@assay$logcount)[var0,], as.matrix(data4@assay$logcount)[var0,], 
  as.matrix(data5@assay$logcount)[var0,], as.matrix(data6@assay$logcount)[var0,]
)
pca0 = irlba(logmat0, nv = 12)$v
umap0 = umap(pca0)$layout
ann0 = read.table(file = paste0(PATH, "/GSE114727_Anno_All.tsv"), sep = "\t", header = T, stringsAsFactors = F)
m0 = data.frame(UMAP1 = umap0[,1], UMAP2 = umap0[,2], Set = ann0$Set, Type = ann0$Type)
m0$Type = factor(m0$Type, levels = c(
  "CD4+ sub1", "CD4+ sub1 TN", "CD4+ sub2", "CD4+ sub3", "CD4+ sub4", 
  "CD4+ Treg", "CD8+ sub1", "CD8+ sub2", "CD8+ sub3", "CD8+ Trm"
))
m0$Set0 = factor(
  m0$Set, 
  levels = c("BC_ER", "BC1_ER_PR", "BC2_ER_PR", "BC1_Her2", "BC2_Her2", "BC_TN"), 
  labels = c("BC ER+", "BC ER+PR+", "BC ER+PR+", "BC Her2+", "BC Her2+", "BC TN")
)

ggplot(m0, aes(UMAP1, UMAP2)) + 
  geom_point(aes(color = Set0, shape = Set), size = 1, alpha = 1) + 
  scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A")) + 
  scale_shape_manual(values = c(1, 2, 4, 5, 6, 8)) + 
  theme_bw(base_line_size = 0) + 
  labs(x = 'UMAP-1', y = 'UMAP-2', color = 'Set', shape = 'Source') + 
  theme(plot.title = element_text(hjust = 0.5, size = 16)) + 
  guides(colour = guide_legend(override.aes = list(size = 5), ncol = 1), shape = guide_legend(override.aes = list(size = 5), ncol = 1))

ggplot(m0, aes(UMAP1, UMAP2)) + 
  geom_point(aes(color = Type, shape = Set), size = 0.5, alpha = 1) + 
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A", "#6A3D9A", "#FF7F00", "#A65628", "#F781BF", "#984EA3", "#999999", "#E7298A")) + 
  scale_shape_manual(values = c(1, 2, 4, 5, 6, 8)) + 
  theme_bw(base_line_size = 0) + 
  labs(x = 'UMAP-1', y = 'UMAP-2', color = 'Set', shape = 'Source') + 
  theme(plot.title = element_text(hjust = 0.5, size = 16)) + 
  guides(colour = guide_legend(override.aes = list(size = 5), ncol = 1), shape = guide_legend(override.aes = list(size = 5), ncol = 1))


###################################################################################
### Data Integration ###
###################################################################################
## BC positive
dat.pos = list(data2, data1, data3, data4, data5)
var0 = read.table(file = paste0(PATH, "/var_pos.tsv"), sep = "\t", header = F, stringsAsFactors = F)
var0 = var0$V1
InPlot(dat.pos, var.gene = var0)
dat.pos = scMultiIntegrate(dat.pos, eigens = 12, var.gene = var0, add.Id = c("BC2_ER_PR", "BC1_ER_PR", "BC_ER", "BC1_Her2", "BC2_Her2"), ncore = 4)
dat.pos = scUMAP(dat.pos, npc = 12, use = "PLS")

ann0 = read.table(file = paste0(PATH, "/GSE114727_Anno_All.tsv"), sep = "\t", header = T, stringsAsFactors = F)
ann0 = ann0[ann0$Set != "BC_TN",]
dat.pos@coldata$Type = as.character(ann0$Type)

DimPlot(dat.pos, slot = "cell.umap", colFactor = "Set", size = 0.5)
DimPlot(dat.pos, slot = "cell.umap", colFactor = "Type", size = 0.5)
DimPlot(dat.pos, slot = "cell.umap", genes = c("HAVCR2", "CD8A", "CD4", "FOXP3"), size = 0.2)


## BC positive and triple negative
dat.all = list(data2, data1, data3, data4, data5, data6)
var0 = read.table(file = paste0(PATH, "/var_all.tsv"), sep = "\t", header = F, stringsAsFactors = F)
var0 = var0$V1
InPlot(dat.all, var.gene = var0)
dat.all = scMultiIntegrate(dat.all, eigens = 12, var.gene = var0, add.Id = c("BC2_ER_PR", "BC1_ER_PR", "BC_ER", "BC1_Her2", "BC2_Her2", "BC_TN"), ncore = 4)
dat.all = scUMAP(dat.all, npc = 12, use = "PLS")

ann0 = read.table(file = paste0(PATH, "/GSE114727_Anno_All.tsv"), sep = "\t", header = T, stringsAsFactors = F)
dat.all@coldata$Type = as.character(ann0$Type)

DimPlot(dat.all, slot = "cell.umap", colFactor = "Set", size = 0.5)
DimPlot(dat.all, slot = "cell.umap", colFactor = "Type", size = 0.5)

UMAPlot(dat.all, genes = "HAVCR2", size = 0.5, exp.col = "firebrick2")
UMAPlot(dat.all, genes = "CD8A", size = 0.5, exp.col = "firebrick2")
UMAPlot(dat.all, genes = "CD4", size = 0.5, exp.col = "firebrick2")
UMAPlot(dat.all, genes = "FOXP3", size = 0.5, exp.col = "firebrick2")
UMAPlot(dat.all, genes = "CD40LG", size = 0.5, exp.col = "firebrick2")
UMAPlot(dat.all, genes = "DPP4", size = 0.5, exp.col = "firebrick2")
UMAPlot(dat.all, genes = "CHN1", size = 0.5, exp.col = "firebrick2")
UMAPlot(dat.all, genes = "KRT86", size = 0.5, exp.col = "firebrick2")
UMAPlot(dat.all, genes = "LINC00402", size = 0.5, exp.col = "firebrick2")
UMAPlot(dat.all, genes = "IKZF2", size = 0.5, exp.col = "firebrick2")
UMAPlot(dat.all, genes = "ZNF683", size = 0.5, exp.col = "firebrick2")
UMAPlot(dat.all, genes = "AIF1", size = 0.5, exp.col = "firebrick2")
UMAPlot(dat.all, genes = "PLEK", size = 0.5, exp.col = "firebrick2")


