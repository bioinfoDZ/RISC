####################################################################################
#' Find Cluster Markers
####################################################################################
#' 
#' This is the basic function in RISC, it can identify the cluster markers by
#' comparing samples in the selected cluster to the samples of the rest clusters.
#' Therefore, it is possible one gene labeled as a marker in more than one clusters.
#' Two methods are employed in RISC, one is based on Negative Binomial model 
#' while the other using QuasiPoisson model.
#' 
#' Because log2 cannot handle counts with value 0, we use log1p to calculate average 
#' values of counts and log2 to format fold-change.
#' 
#' @rdname Cluster-Marker
#' @param object RISC object: a framework dataset.
#' @param cluster Select the cluster that we want to detect cluster marker genes.
#' @param positive Whether only output the cluster markers with positive log2FC.
#' @param frac A fraction cutoff, the cluster marekr genes expressed at least a 
#' cutoff fraction of the cluster cells.
#' @param log2FC The cutoff of log2 Fold-change for differentially expressed marker 
#' genes.
#' @param Padj The cutoff of the adjusted P-value.
#' @param latent.factor The latent factor from coldata, which represents number 
#' values or factors, and only one latent factor can be inputed.
#' @param method Which method is used to identify cluster markers, three options: 'NB' 
#' for Negative Binomial model, 'QP' for QuasiPoisson model, and 'Wilcox' for Wilcoxon 
#' Rank Sum and Signed Rank model.
#' @param min.cells The minimum cells for each cluster to calculate marker genes.
#' @param ncore The multiple cores for parallel calculating.
#' @param fast In the fast way, making cell bins to calculate marker genes.
#' @param bin_bg When run fast, how many cell bins will use in background cells.
#' @param bin_cl When run fast, how many cell bins will use in cluster cells.
#' @importFrom Matrix.utils aggregate.Matrix
#' @importFrom MASS glm.nb
#' @importFrom stats var wilcox.test
#' @references Paternoster et al., Criminology (1997)
#' @references Berk et al., Journal of Quantitative Criminology (2008)
#' @references Liu et al., Nature Biotech. (2021)
#' @name scMarker
#' @export
#' @examples 
#' # RISC object
#' obj0 = raw.mat[[3]]
#' obj0 = scPCA(obj0, npc = 10)
#' obj0 = scUMAP(obj0, npc = 3)
#' obj0 = scCluster(obj0, slot = "cell.umap", k = 3, method = 'density')
#' DimPlot(obj0, slot = "cell.umap", colFactor = 'Cluster', size = 2)
#' marker1 = scMarker(obj0, cluster = 1, method = 'QP', min.cells = 3)

scMarker <- function(
  object, 
  cluster = 1, 
  positive = TRUE, 
  frac = 0.25, 
  log2FC = 0.5, 
  Padj = 0.05, 
  latent.factor = NULL, 
  method = 'QP', 
  min.cells = 10, 
  ncore = 1, 
  fast = FALSE, 
  bin_bg = 500, 
  bin_cl = 100
  ) {
  
  pboptions(type = "txt", style = 3)
  count = object@assay$logcount
  coldata = as.data.frame(object@coldata)
  latent.factor0 = intersect(colnames(coldata), latent.factor[1])
  ncore = as.integer(ncore)
  registerDoParallel(ncore)
  Padj = as.numeric(Padj)
  log2FC = as.numeric(log2FC)
  min.cells = as.integer(min.cells)
  bin_bg = as.integer(bin_bg)
  bin_cl = as.integer(bin_cl)
  
  cell1 = rownames(coldata)[coldata$Cluster %in% cluster]
  cell2 = rownames(coldata)[!coldata$Cluster %in% cluster]
  
  if(length(cell1) < min.cells | length(cell2) < min.cells){
    stop("The number of cells in either control or sample cluster is too low")
  } else {
    cell1 = cell1
    cell2 = cell2
  }
  
  
  if('Integration' %in% names(object@metadata)){
    
    coldata$Set = factor(coldata$Set)
    set0 = levels(coldata$Set)
    rowsum0 = lapply(count, FUN = function(y){Matrix::rowSums(y[, colnames(y) %in% cell1, drop = FALSE] > 0)})
    rowsum0 = do.call(cbind, rowsum0)
    percent0 = round(Matrix::rowSums(rowsum0) / length(cell1), digits = 2)
    keep = Matrix::rowSums(rowsum0) > floor(length(cell1) * frac)
    gene0 = rownames(object@rowdata)[keep]
    percent0 = percent0[keep]
    count = lapply(count, FUN = function(y){y[gene0, , drop = FALSE]})
    
    cell1 = data.frame(id = cell1, seq = 1L:length(cell1))
    rownames(cell1) = cell1$id
    cell1$set = as.character(coldata[rownames(cell1), "Set"])
    set_bin = data.frame(table(cell1$set))
    colnames(set_bin) = c("Label", "Value")
    set_bin$Label = as.character(set_bin$Label)
    set_bin$Value = as.integer(set_bin$Value)
    set_bin = set_bin[set_bin$Value != 0,]
    
    ave0 = lapply(1L:nrow(set_bin), FUN = function(i){spr_count(count0 = count, i = i, id0 = set_bin$Label, cell0 = cell1$id)})
    ave0 = lapply(ave0, FUN = function(y){Matrix::rowMeans(y) * ncol(y)})
    ave0 = do.call(cbind, ave0)
    ave0 = Matrix::rowSums(ave0) / nrow(cell1)
    
    if(fast & length(cell2) > bin_bg){
      
      if(nrow(cell1) <= bin_cl){
        
        count1 = lapply(1L:nrow(set_bin), FUN = function(i){spr_count(count0 = count, i = i, id0 = set_bin$Label, cell0 = cell1$id)})
        count1 = do.call(cbind, count1)
        colnames(count1) = cell1 = paste0("sam_", 1L:ncol(count1))
        
      } else {
        
        bin1 = lapply(1L:nrow(set_bin), FUN = function(i){spr_bin(cell0 = cell1, id0 = set_bin, i = i, bins = bin_cl)})
        names(bin1) = set_bin$Label
        count1 = lapply(1L:nrow(set_bin), FUN = function(i){aggregate.Matrix(x = Matrix::t(spr_count(count0 = count, i = i, id0 = set_bin$Label, cell0 = cell1$id)), groupings = bin1[[i]], fun = "mean")})
        names(count1) = set_bin$Label
        count1 = lapply(count1, FUN = function(x){spr_na_t(x)})
        count1 = do.call(cbind, count1)
        rownames(count1) = gene0
        colnames(count1) = cell1 = paste0("sam_", 1L:ncol(count1))
        
      }
      
      cell2 = data.frame(id = cell2, seq = 1L:length(cell2))
      rownames(cell2) = cell2$id
      cell2$set = as.character(coldata[rownames(cell2), "Set"])
      set_bin = data.frame(table(cell2$set))
      colnames(set_bin) = c("Label", "Value")
      set_bin$Label = as.character(set_bin$Label)
      set_bin$Value = as.integer(set_bin$Value)
      set_bin = set_bin[set_bin$Value != 0,]
      
      bin2 = lapply(1L:nrow(set_bin), FUN = function(i){spr_bin(cell0 = cell2, id0 = set_bin, i = i, bins = bin_bg)})
      names(bin2) = set_bin$Label
      count2 = lapply(1L:nrow(set_bin), FUN = function(i){aggregate.Matrix(x = Matrix::t(spr_count(count0 = count, i = i, id0 = set_bin$Label, cell0 = cell2$id)), groupings = bin2[[i]], fun = "mean")})
      names(count2) = set_bin$Label
      count2 = lapply(count2, FUN = function(x){spr_na_t(x)})
      count2 = do.call(cbind, count2)
      rownames(count2) = gene0
      colnames(count2) = cell2 = paste0("ctrl_", 1L:ncol(count2))
      
      count = cbind(count2, count1)
      group1 = rep('sam', length(cell1))
      group2 = rep('ctrl', length(cell2))
      group = factor(c(group2, group1), levels = c('ctrl', 'sam'))
      group = as.integer(group)
      latent.factor0 = NULL
      
    } else {
      
      group = rep('ctrl', nrow(coldata))
      group[rownames(coldata) %in% rownames(cell1)] = 'sam'
      group = factor(group, levels = c('ctrl', 'sam'))
      group = as.integer(group)
      fast = FALSE
      
    }
    
    if(length(coldata$Cluster) == 0){
      stop('Please cluster cells first.')
    } else if(length(latent.factor0) == 0 & fast) {
      
      if(method == 'NB'){
        result = pblapply(1L:length(gene0), FUN = function(x){NB.detect(group, count[x,])}, cl = ncore)
      } else if(method == 'QP') {
        result = pblapply(1L:length(gene0), FUN = function(x){QP.detect(group, count[x,])}, cl = ncore)
      } else {
        result = pblapply(1L:length(gene0), FUN = function(x){wil.detect(group, count[x,])}, cl = ncore)
      }
      
    } else if(length(latent.factor0) == 0){
      
      if(method == 'NB'){
        result = pblapply(1L:length(gene0), FUN = function(x){NB.detect(group, unlist(lapply(count, FUN = function(y){y[x,]})))}, cl = ncore)
      } else if(method == 'QP') {
        result = pblapply(1L:length(gene0), FUN = function(x){QP.detect(group, unlist(lapply(count, FUN = function(y){y[x,]})))}, cl = ncore)
      } else {
        result = pblapply(1L:length(gene0), FUN = function(x){wil.detect(group, unlist(lapply(count, FUN = function(y){y[x,]})))}, cl = ncore)
      }
      
    } else if(length(latent.factor0) == 1){
      
      coldata0 = data.frame(coldata[,latent.factor0])
      
      if(method == 'NB'){
        result = pblapply(1L:length(gene0), FUN = function(x){NB.latent.detect(unlist(lapply(count, FUN = function(y){y[x,]})), group, coldata0)}, cl = ncore)
      } else {
        result = pblapply(1L:length(gene0), FUN = function(x){QP.latent.detect(unlist(lapply(count, FUN = function(y){y[x,]})), group, coldata0)}, cl = ncore)
      }
      
    } else {
      stop('Check input')
      }
    
  } else {
    
    rowsum0 = Matrix::rowSums(count[, colnames(count) %in% cell1, drop = FALSE] > 0)
    percent0 = round(rowsum0 / length(cell1), digits = 2)
    keep = rowsum0 > floor(length(cell1) * frac)
    percent0 = percent0[keep]
    count = count[keep, , drop = FALSE]
    gene0 = rownames(object@rowdata)[keep]
    ave0 = Matrix::rowMeans(count[, colnames(count) %in% cell1, drop = FALSE])
    
    if(fast & length(cell2) > bin_bg){
      
      cell1 = data.frame(id = cell1, seq = 1L:length(cell1))
      rownames(cell1) = cell1$id
      
      if(nrow(cell1) <= bin_cl){
        count1 = count[, colnames(count) %in% cell1$id, drop = FALSE]
        colnames(count1) = cell1 = paste0("sam_", 1L:ncol(count1))
      } else {
        bin1 = cut(cell1$seq, breaks = bin_cl, labels = FALSE)
        cell1$bin = as.factor(bin1)
        count1 = aggregate.Matrix(x = Matrix::t(count[, colnames(count) %in% cell1$id, drop = FALSE]), groupings = cell1$bin, fun = "mean")
        count1 = Matrix::t(count1)
        rownames(count1) = rownames(count)
        colnames(count1) = cell1 = paste0("sam_", 1L:ncol(count1))
      }
      
      cell2 = data.frame(id = cell2, seq = 1L:length(cell2))
      rownames(cell2) = cell2$id
      
      bin2 = cut(cell2$seq, breaks = bin_bg, labels = FALSE)
      cell2$bin = as.factor(bin2)
      count2 = aggregate.Matrix(x = Matrix::t(count[, colnames(count) %in% cell2$id, drop = FALSE]), groupings = cell2$bin, fun = "mean")
      count2 = Matrix::t(count2)
      rownames(count2) = rownames(count)
      colnames(count2) = cell2 = paste0("ctrl_", 1L:ncol(count2))
      
      count = cbind(count2, count1)
      group1 = rep('sam', length(cell1))
      group2 = rep('ctrl', length(cell2))
      group = factor(c(group2, group1), levels = c('ctrl', 'sam'))
      group = as.integer(group)
      latent.factor0 = NULL
      
    } else {
      
      group = rep('ctrl', nrow(coldata))
      group[rownames(coldata) %in% cell1] = 'sam'
      group = factor(group, levels = c('ctrl', 'sam'))
      group = as.integer(group)
      fast = FALSE
      
    }
    
    if(length(coldata$Cluster) == 0){
      stop('Please cluster cells first.')
    } else if(length(latent.factor0) == 0){
      
      if(method == 'NB'){
        result = pblapply(1L:length(gene0), FUN = function(x){NB.detect(group, count[x,])}, cl = ncore)
      } else if(method == 'QP') {
        result = pblapply(1L:length(gene0), FUN = function(x){QP.detect(group, count[x,])}, cl = ncore)
      } else {
        result = pblapply(1L:length(gene0), FUN = function(x){wil.detect(group, count[x,])}, cl = ncore)
      }
      
    } else if(length(latent.factor0) == 1){
      
      coldata0 = data.frame(coldata[,latent.factor0])
      
      if(method == 'NB'){
        result = pblapply(1L:length(gene0), FUN = function(x){NB.latent.detect(count[x,], group, coldata0)}, cl = ncore)
      } else {
        result = pblapply(1L:length(gene0), FUN = function(x){QP.latent.detect(count[x,], group, coldata0)}, cl = ncore)
      }
      
    } else {
      stop('Check input')
    }
    
  }
  
  result = do.call(rbind, result)
  result = as.matrix(result)
  
  colnames(result) = c('log2FC', 'Pvalue')
  result = data.frame(Symbol = gene0, Percent = percent0, logAve = ave0, result)
  result$Padj = p.adjust(result$Pvalue, method = 'BH')
  result = result[order(result$Pvalue, decreasing = FALSE),]
  # result$logAve = sapply(result$logAve, FUN = function(x){log1p_2(x)})
  
  if(is.null(Padj)){
    result0 = result[result$Pvalue < 0.05,]
  } else {
    if(Padj == 1) {
      result0 = result
    } else {
      result0 = result[result$Padj < Padj,]
    }
  }
  
  if(positive) {
    result0 = result0[result0$log2FC > log2FC,]
  } else {
    result0 = result0[abs(result0$log2FC) > log2FC,]
  }
  
  return(result0)
  
}



####################################################################################
#' Find All Cluster Markers
####################################################################################
#' 
#' This function depends on "scMarker" function by using the same criteria, and 
#' generates markers for all the clusters. Here, if the cell number of any cluster 
#' is less than 10 , RISC will skip this cluster and not detect its cluster markers 
#' in the default parameters.
#' 
#' Because log2 cannot handle counts with value 0, we use log1p to calculate average 
#' values of counts and log2 to format fold-change. 
#' 
#' @rdname All-Cluster-Marker
#' @param object RISC object: a framework dataset.
#' @param positive Whether only output the cluster markers with positive log2FC.
#' @param frac A fraction cutoff, the cluster marker genes expressed at least a 
#' cutoff fraction of the cluster cells.
#' @param log2FC The cutoff of log2 Fold-change for differentially expressed marker 
#' genes.
#' @param Padj The cutoff of the adjusted P-value.
#' @param latent.factor The latent factor from coldata, which represents number 
#' values or factors, and only one latent factor can be inputed.
#' @param min.cells The threshold for the cell number of valid clusters.
#' @param method Which method is used to identify cluster markers, three options: 'NB' 
#' for Negative Binomial model, 'QP' for QuasiPoisson model, and 'Wilcox' for Wilcoxon 
#' Rank Sum and Signed Rank model.
#' @param ncore The multiple cores for parallel calculating.
#' @param fast In the fast way, making cell bins to calculate marker genes.
#' @param bin_bg When run fast, how many cell bins will use in background cells.
#' @param bin_cl When run fast, how many cell bins will use in cluster cells.
#' @importFrom MASS glm.nb
#' @name AllMarker
#' @export

AllMarker <- function(
  object, 
  positive = TRUE, 
  frac = 0.25, 
  log2FC = 0.5, 
  Padj = 0.05, 
  latent.factor = NULL, 
  min.cells = 25L, 
  method = 'QP', 
  ncore = 1, 
  fast = FALSE, 
  bin_bg = 500, 
  bin_cl = 100
  ) {
  
  if(length(object@cluster) == 0) {
    stop('Please cluster cells first')
  } else {
    
    Padj = as.numeric(Padj)
    log2FC = as.numeric(log2FC)
    ncore = as.numeric(ncore)
    
    coldata = as.data.frame(object@coldata)
    cluster0 = data.frame(table(coldata$Cluster))
    colnames(cluster0) = c('Cluster', 'No')
    
    if(length(bin_bg) == nrow(cluster0)){
      bin_bg = as.integer(bin_bg)
    } else {
      bin_bg = rep(bin_bg[1], nrow(cluster0))
    }
    
    if(length(bin_cl) == nrow(cluster0)){
      bin_cl = as.integer(bin_cl)
    } else {
      bin_cl = rep(bin_cl[1], nrow(cluster0))
    }
    
    keep = (cluster0$No >= min.cells)
    bin_bg = bin_bg[keep]
    bin_cl = bin_cl[keep]
    cluster0 = cluster0[keep,]
    result = data.frame(Symbol = 'Gene', Percent = 0, logAve = 0, log2FC = 0, Pvalue = 1, Padj = 1, Cluster = 0)
    
    if(length(bin_bg) == nrow(cluster0))
    
    i = 1L
    while(i <= nrow(cluster0)) {
      
      clusteri = levels(as.factor(cluster0$Cluster))[i]
      cell1 = rownames(coldata)[coldata$Cluster %in% clusteri]
      cell2 = rownames(coldata)[!coldata$Cluster %in% clusteri]
      
      if(length(cell1) < min.cells | length(cell2) < min.cells){
        
        result = result
        i = i + 1L
        
      } else {
        
        resulti = scMarker(object = object, cluster = clusteri, positive = positive, frac = frac, log2FC = log2FC, Padj = Padj, latent.factor = latent.factor, method = method, ncore = ncore, min.cells = min.cells, fast = fast, bin_bg = bin_bg[i], bin_cl = bin_cl[i])
        
        if(nrow(resulti) > 0) {
          resulti$Cluster = clusteri
          result = rbind(result, resulti)
          i = i + 1L
        } else {
          result = result
          i = i + 1L
        }
        
      }
      
    }
    
    result0 = result[-1L,]
    
  }
  
  return(result0)
  
}



####################################################################################
#' Find Differentially Expressed Genes between Clusters
####################################################################################
#' 
#' This is the basic function in RISC, it can identify the differentially expressed 
#' genes (DEGs) by comparing samples between the selected clusters. The criteria 
#' used for the cluster markers are also appropriate to DEGs.
#' 
#' Here RISC provides two algorithms to detect DEGs, the primary one is a model 
#' "Quasi-Poisson" which has advantage to identify DEGs from the cluster with 
#' a small number of cells. Meanwhile, RISC also has alternative algorithm: 
#' "Negative Binomial" model.
#' 
#' Because log2 cannot handle counts with value 0, we use log1p to calculate average 
#' values of counts and log2 to format fold-change.
#' 
#' @rdname Pairwise-DEGs
#' @param object RISC object: a framework dataset.
#' @param cell.ctrl Select the cells as the reference cells for detecting DEGs.
#' @param cell.sam Select the cells as the sample cells for detecting DEGs.
#' @param frac A fraction cutoff, the cluster marker genes expressed at least a 
#' cutoff fraction of the cluster cells.
#' @param log2FC The cutoff of log2 Fold-change for differentially expressed marker 
#' genes.
#' @param Padj The cutoff of the adjusted P-value. If Padj is NULL, use p-value 
#' < 0.05 as a threshold. Set Padj as 1, without any cutoff.
#' @param latent.factor The latent factor from coldata, which represents number 
#' values or factors, and only one latent factor can be inputed.
#' @param method Which method is used to identify cluster markers, two options: 
#' 'NB' for Negative Binomial model, 'QP' for QuasiPoisson model, and 'wil' for
#' Wilcoxon Rank-Sum model.
#' @param min.cells The minimum cells for each cluster to calculate marker genes.
#' @param ncore The multiple cores for parallel calculating.
#' @param fast In the fast way, making cell bins to calculate marker genes.
#' @param bin_ctrl When run fast, how many cell bins will use in control cells.
#' @param bin_sam When run fast, how many cell bins will use in sample cells.
#' @references Paternoster et al., Criminology (1997)
#' @references Berk et al., Journal of Quantitative Criminology (2008)
#' @references Liu et al., Nature Biotech. (2021)
#' @name scDEG
#' @export
#' @examples
#' # RISC object
#' obj0 = raw.mat[[4]]
#' obj0 = scPCA(obj0, npc = 10)
#' obj0 = scUMAP(obj0, npc = 3)
#' obj0 = scCluster(obj0, slot = "cell.umap", k = 3, method = 'density')
#' DimPlot(obj0, slot = "cell.umap", colFactor = 'Cluster', size = 2)
#' cell.ctrl = rownames(obj0@coldata)[obj0@coldata$Cluster == 1]
#' cell.sam = rownames(obj0@coldata)[obj0@coldata$Cluster == 3]
#' DEG0 = scDEG(obj0, cell.ctrl = cell.ctrl, cell.sam = cell.sam, 
#'              min.cells = 3, method = 'QP')

scDEG <- function(
  object, 
  cell.ctrl = NULL, 
  cell.sam = NULL, 
  frac = 0.10, 
  log2FC = 1, 
  Padj = 0.01, 
  latent.factor = NULL, 
  method = 'NB', 
  min.cells = 10, 
  ncore = 1, 
  fast = FALSE, 
  bin_ctrl = 500, 
  bin_sam = 100
  ) {
  
  pboptions(type = "txt", style = 3)
  count = object@assay$logcount
  coldata = as.data.frame(object@coldata)
  latent.factor0 = intersect(colnames(coldata), latent.factor)
  ncore = as.numeric(ncore)
  registerDoParallel(ncore)
  log2FC = as.numeric(log2FC)
  min.cells = as.numeric(min.cells)
  cell.ctrl = intersect(cell.ctrl, rownames(coldata))
  cell.sam = intersect(cell.sam, rownames(coldata))
  bin_ctrl = as.integer(bin_ctrl)
  bin_sam = as.integer(bin_sam)
  
  if(length(cell.ctrl) < min.cells | length(cell.sam) < min.cells){
    stop("The number of cells in either control or sample is too low")
  } else {
    cell1 = cell.sam
    cell2 = cell.ctrl
  }
  
  coldata = coldata[rownames(coldata) %in% unique(c(cell1, cell2)),]
  
  if('Integration' %in% names(object@metadata)){
    coldata$Set = as.character(coldata$Set)
    if(length(table(coldata$Set)) > 1) {
      coldata$Set = factor(coldata$Set)
      set0 = levels(coldata$Set)
      Inte = TRUE
    } else {
      set0 = unique(coldata$Set)
      count = count[[set0]]
      Inte = FALSE
    }
  } else {
    Inte = FALSE
  }
  
  if(Inte){
    
    count = lapply(count, FUN = function(y){y[, colnames(y) %in% unique(c(cell1, cell2)), drop = FALSE]})
    rowsum0 = lapply(count, FUN = function(y){Matrix::rowSums(y > 0)})
    rowsum0 = do.call(cbind, rowsum0)
    percent0 = round(Matrix::rowSums(rowsum0) / (length(cell1) + length(cell2)), digits = 2)
    keep = Matrix::rowSums(rowsum0) > floor(min(length(cell1), length(cell2)) * frac)
    percent0 = percent0[keep]
    gene0 = rownames(object@rowdata)[keep]
    count = lapply(count, FUN = function(y){y[gene0, , drop = FALSE]})
    
    cell1 = data.frame(id = cell1, seq = 1L:length(cell1))
    rownames(cell1) = cell1$id
    cell1$set = as.character(coldata[rownames(cell1), "Set"])
    set_bin1 = data.frame(table(cell1$set))
    colnames(set_bin1) = c("Label", "Value")
    set_bin1$Label = as.character(set_bin1$Label)
    set_bin1$Value = as.integer(set_bin1$Value)
    set_bin1 = set_bin1[set_bin1$Value != 0,]
    
    ave1 = lapply(1L:nrow(set_bin1), FUN = function(i){spr_count(count0 = count, i = i, id0 = set_bin1$Label, cell0 = cell1$id)})
    ave1 = lapply(ave1, FUN = function(y){Matrix::rowMeans(y) * ncol(y)})
    ave1 = do.call(cbind, ave1)
    ave1 = Matrix::rowSums(ave1) / nrow(cell1)
    
    cell2 = data.frame(id = cell2, seq = 1L:length(cell2))
    rownames(cell2) = cell2$id
    cell2$set = as.character(coldata[rownames(cell2), "Set"])
    set_bin2 = data.frame(table(cell2$set))
    colnames(set_bin2) = c("Label", "Value")
    set_bin2$Label = as.character(set_bin2$Label)
    set_bin2$Value = as.integer(set_bin2$Value)
    set_bin2 = set_bin2[set_bin2$Value != 0,]
    
    ave2 = lapply(1L:nrow(set_bin2), FUN = function(i){spr_count(count0 = count, i = i, id0 = set_bin2$Label, cell0 = cell2$id)})
    ave2 = lapply(ave2, FUN = function(y){Matrix::rowMeans(y) * ncol(y)})
    ave2 = do.call(cbind, ave2)
    ave2 = Matrix::rowSums(ave2) / nrow(cell2)
    
    if(fast & length(cell2) > bin_ctrl){
      
      if(length(cell1) <= bin_sam){
        
        count1 = lapply(1L:nrow(set_bin1), FUN = function(i){spr_count(count0 = count, i = i, id0 = set_bin1$Label, cell0 = cell1$id)})
        count1 = do.call(cbind, count1)
        colnames(count1) = cell1 = paste0("sam_", 1L:ncol(count1))
        
      } else {
        
        bin1 = lapply(1L:nrow(set_bin1), FUN = function(i){spr_bin(cell0 = cell1, id0 = set_bin1, i = i, bins = bin_sam)})
        names(bin1) = set_bin1$Label
        count1 = lapply(1L:nrow(set_bin1), FUN = function(i){aggregate.Matrix(x = Matrix::t(spr_count(count0 = count, i = i, id0 = set_bin1$Label, cell0 = cell1$id)), groupings = bin1[[i]], fun = "mean")})
        names(count1) = set_bin1$Label
        count1 = lapply(count1, FUN = function(x){spr_na_t(x)})
        count1 = do.call(cbind, count1)
        rownames(count1) = gene0
        colnames(count1) = cell1 = paste0("sam_", 1L:ncol(count1))
        
      }
      
      bin2 = lapply(1L:nrow(set_bin2), FUN = function(i){spr_bin(cell0 = cell2, id0 = set_bin2, i = i, bins = bin_ctrl)})
      names(bin2) = set_bin2$Label
      count2 = lapply(1L:nrow(set_bin2), FUN = function(i){aggregate.Matrix(x = Matrix::t(spr_count(count0 = count, i = i, id0 = set_bin2$Label, cell0 = cell2$id)), groupings = bin2[[i]], fun = "mean")})
      names(count2) = set_bin2$Label
      count2 = lapply(count2, FUN = function(x){spr_na_t(x)})
      count2 = do.call(cbind, count2)
      rownames(count2) = gene0
      colnames(count2) = cell2 = paste0("ctrl_", 1L:ncol(count2))
      
      count = cbind(count2, count1)
      group1 = rep('sam', length(cell1))
      group2 = rep('ctrl', length(cell2))
      group = factor(c(group2, group1), levels = c('ctrl', 'sam'))
      group = as.integer(group)
      latent.factor0 = NULL
      
    } else {
      
      group = rep('ctrl', nrow(coldata))
      group[rownames(coldata) %in% rownames(cell1)] = 'sam'
      group = factor(group, levels = c('ctrl', 'sam'))
      group = as.integer(group)
      fast = FALSE
      
    }
    
    if(length(latent.factor0) == 0 & fast) {
      
      if(method == 'NB'){
        result = pblapply(1L:length(gene0), FUN = function(x){NB.detect(group, count[x,])}, cl = ncore)
      } else if(method == 'QP') {
        result = pblapply(1L:length(gene0), FUN = function(x){QP.detect(group, count[x,])}, cl = ncore)
      } else {
        result = pblapply(1L:length(gene0), FUN = function(x){wil.detect(group, count[x,])}, cl = ncore)
      }
      
    } else if(length(latent.factor0) == 0){
      
      if(method == 'NB'){
        result = pblapply(1L:length(gene0), FUN = function(x){NB.detect(group, unlist(lapply(count, FUN = function(y){y[x,]})))}, cl = ncore)
      } else if(method == 'QP') {
        result = pblapply(1L:length(gene0), FUN = function(x){QP.detect(group, unlist(lapply(count, FUN = function(y){y[x,]})))}, cl = ncore)
      } else {
        result = pblapply(1L:length(gene0), FUN = function(x){wil.detect(group, unlist(lapply(count, FUN = function(y){y[x,]})))}, cl = ncore)
      }
      
    } else if(length(latent.factor0) == 1){
      
      coldata0 = data.frame(coldata[,latent.factor0])
      
      if(method == 'NB'){
        result = pblapply(1L:length(gene0), FUN = function(x){NB.latent.detect(unlist(lapply(count, FUN = function(y){y[x,]})), group, coldata0)}, cl = ncore)
      } else {
        result = pblapply(1L:length(gene0), FUN = function(x){QP.latent.detect(unlist(lapply(count, FUN = function(y){y[x,]})), group, coldata0)}, cl = ncore)
      }
      
    } else {
      stop('Check input')
    }
    
  } else {
    
    count = count[, colnames(count) %in% unique(c(cell1, cell2)), drop = FALSE]
    rowsum0 = Matrix::rowSums(count > 0)
    percent0 = round(rowsum0 / (length(cell1) + length(cell2)), digits = 2)
    keep = rowsum0 > floor(min(length(cell1), length(cell2)) * frac)
    percent0 = percent0[keep]
    count = count[keep, , drop = FALSE]
    gene0 = rownames(object@rowdata)[keep]
    ave1 = Matrix::rowMeans(count[, colnames(count) %in% cell1, drop = FALSE])
    ave2 = Matrix::rowMeans(count[, colnames(count) %in% cell2, drop = FALSE])
    
    if(fast & length(cell2) > bin_ctrl){
      
      cell1 = data.frame(id = cell1, seq = 1L:length(cell1))
      rownames(cell1) = cell1$id
      
      if(length(cell1) <= bin_sam){
        count1 = count[, colnames(count) %in% cell1$id, drop = FALSE]
        colnames(count1) = cell1 = paste0("sam_", 1L:ncol(count1))
      } else {
        bin1 = cut(cell1$seq, breaks = bin_sam, labels = FALSE)
        cell1$bin = as.factor(bin1)
        count1 = aggregate.Matrix(x = Matrix::t(count[, colnames(count) %in% cell1$id, drop = FALSE]), groupings = cell1$bin, fun = "mean")
        count1 = Matrix::t(count1)
        rownames(count1) = rownames(count)
        colnames(count1) = cell1 = paste0("sam_", 1L:ncol(count1))
      }
      
      cell2 = data.frame(id = cell2, seq = 1L:length(cell2))
      rownames(cell2) = cell2$id
      
      bin2 = cut(cell2$seq, breaks = bin_ctrl, labels = FALSE)
      cell2$bin = as.factor(bin2)
      count2 = aggregate.Matrix(x = Matrix::t(count[, colnames(count) %in% cell2$id, drop = FALSE]), groupings = cell2$bin, fun = "mean")
      count2 = Matrix::t(count2)
      rownames(count2) = rownames(count)
      colnames(count2) = cell2 = paste0("ctrl_", 1L:ncol(count2))
      
      count = cbind(count2, count1)
      group1 = rep('sam', length(cell1))
      group2 = rep('ctrl', length(cell2))
      group = factor(c(group2, group1), levels = c('ctrl', 'sam'))
      group = as.integer(group)
      latent.factor0 = NULL
      
    } else {
      
      group = rep('ctrl', nrow(coldata))
      group[rownames(coldata) %in% cell1] = 'sam'
      group = factor(group, levels = c('ctrl', 'sam'))
      group = as.integer(group)
      fast = FALSE
      
    }
    
    if(length(latent.factor0) == 0){
      
      if(method == 'NB'){
        result = pblapply(1L:length(gene0), FUN = function(x){NB.detect(group, count[x,])}, cl = ncore)
      } else if(method == 'QP') {
        result = pblapply(1L:length(gene0), FUN = function(x){QP.detect(group, count[x,])}, cl = ncore)
      } else {
        result = pblapply(1L:length(gene0), FUN = function(x){wil.detect(group, count[x,])}, cl = ncore)
      }
      
    } else if(length(latent.factor0) == 1){
      
      coldata0 = data.frame(coldata[,latent.factor0])
      
      if(method == 'NB'){
        result = pblapply(1L:length(gene0), FUN = function(x){NB.latent.detect(count[x,], group, coldata0)}, cl = ncore)
      } else {
        result = pblapply(1L:length(gene0), FUN = function(x){QP.latent.detect(count[x,], group, coldata0)}, cl = ncore)
      }
      
    } else {
      stop('Check input')
    }
    
  }
  
  result = do.call(rbind, result)
  result = as.matrix(result)
  
  colnames(result) = c('log2FC', 'Pvalue')
  result = data.frame(Symbol = gene0, Percent = percent0, avelogCtl = ave2, avelogSam = ave1, result)
  result$Padj = p.adjust(result$Pvalue, method = 'BH')
  result = result[order(result$Pvalue, decreasing = FALSE),]
  # result$avelogCtl = sapply(result$avelogCtl, FUN = function(x){log1p_2(x)})
  # result$avelogSam = sapply(result$avelogSam, FUN = function(x){log1p_2(x)})
  
  if(is.null(Padj)){
    result0 = result[result$Pvalue < 0.05 & abs(result$log2FC) >= log2FC,]
  } else {
    Padj = as.numeric(Padj)
    
    if(Padj == 1) {
      result0 = result
    } else {
      result0 = result[result$Padj < Padj & abs(result$log2FC) >= log2FC,]
    }
    
  }
  
  return(result0)
  
}



####################################################################################
####################################################################################
## Negative Binomial
# NB.core <- function(x0, y0){
#   l0 = glm(y0 ~ x0, family = MASS::negative.binomial(1))
#   l1 = c(coef(summary(l0))[2, 1], coef(summary(l0))[2, 4])
#   return(l1)
# }

NB.core <- function(x0, y0){
  logfc0 = aggregate(y0, list(x0), mean)
  logfc0[1, 2] = log1p_2(logfc0[1, 2])
  logfc0[2, 2] = log1p_2(logfc0[2, 2])
  logfc0 = logfc0[2, 2] - logfc0[1, 2]
  l0 = suppressWarnings(glm.nb(y0 ~ x0))
  res0 = c(logfc0, summary(l0)$coefficients[2, 4])
  return(res0)
}

NB.latent.core <- function(y0, x0, z0){
  logfc0 = aggregate(y0, list(x0), mean)
  logfc0[1, 2] = log1p_2(logfc0[1, 2])
  logfc0[2, 2] = log1p_2(logfc0[2, 2])
  logfc0 = logfc0[2, 2] - logfc0[1, 2]
  l0 = suppressWarnings(glm.nb(y0 ~ z0 + x0))
  res0 = c(logfc0, summary(l0)$coefficients[3, 4])
  return(res0)
}

NB.detect <- function(x0, y0){
  return(
    tryCatch(
      NB.core(x0, y0), 
      error = function(e){c(0, 1)}
    )
  )
}

NB.latent.detect <- function(y0, x0, z0){
  return(
    tryCatch(
      NB.latent.core(y0, x0, z0), 
      error = function(e){c(0, 1)}
    )
  )
}


## QuasiPoisson
# QP.core <- function(x0, y0){
#   l0 = fastglm(x = x0, y = y0, family = quasipoisson(), method = 2)
#   l1 = c(coef(summary(l0))[2, 1], coef(summary(l0))[2, 4])
#   return(l1)
# }

QP.core <- function(x0, y0){
  logfc0 = aggregate(y0, list(x0), mean)
  logfc0[1, 2] = log1p_2(logfc0[1, 2])
  logfc0[2, 2] = log1p_2(logfc0[2, 2])
  logfc0 = logfc0[2, 2] - logfc0[1, 2]
  l0 = suppressWarnings(glm(y0 ~ x0, family = quasipoisson))
  res0 = c(logfc0, summary(l0)$coefficients[2, 4])
  return(res0)
}

QP.latent.core <- function(y0, x0, z0){
  logfc0 = aggregate(y0, list(x0), mean)
  logfc0[1, 2] = log1p_2(logfc0[1, 2])
  logfc0[2, 2] = log1p_2(logfc0[2, 2])
  logfc0 = logfc0[2, 2] - logfc0[1, 2]
  l0 = suppressWarnings(glm(y0 ~ z0 + x0, family = quasipoisson))
  res0 = c(logfc0, summary(l0)$coefficients[3, 4])
  return(res0)
}

QP.detect <- function(x0, y0){
  return(
    tryCatch(
      QP.core(x0, y0), 
      error = function(e){c(0, 1)}
    )
  )
}

QP.latent.detect <- function(y0, x0, z0){
  return(
    tryCatch(
      QP.latent.core(y0, x0, z0), 
      error = function(e){c(0, 1)}
    )
  )
}


## Wilcoxon Rank Sum and Signed Rank
wil.core <- function(x0, y0){
  logfc0 = aggregate(y0, list(x0), mean)
  logfc0[1, 2] = log1p_2(logfc0[1, 2])
  logfc0[2, 2] = log1p_2(logfc0[2, 2])
  logfc0 = logfc0[2, 2] - logfc0[1, 2]
  l0 = suppressWarnings(wilcox.test(y0 ~ x0))
  res0 = c(logfc0, l0$p.value)
  return(res0)
}

wil.detect <- function(x0, y0){
  return(
    tryCatch(
      wil.core(x0, y0), 
      error = function(e){c(0, 1)}
    )
  )
}

spr_na_t <- function(obj0){
  obj0[is.na(obj0)] = 0
  obj0 = as(obj0, "CsparseMatrix")
  obj0 = Matrix::t(obj0)
  return(obj0)
}

spr_count <- function(count0, i, id0, cell0){
  counti = count0[[id0[i]]]
  counti = counti[, colnames(counti) %in% cell0, drop = FALSE]
  return(counti)
}

spr_bin <- function(cell0, id0, i, bins){
  break0 = ceiling(bins * (id0$Value[i] / nrow(cell0)))
  if(break0 > 1){
    bin0 = cut(cell0$seq[cell0$set == id0$Label[i]], breaks = break0, labels = FALSE)
  } else {
    bin0 = rep(1, id0$Value[i])
  }
  return(bin0)
}

log1p_2 <- function(x0){
  if(x0 <= 0){
    x0 = 0
  } else {
    x0 = log2(expm1(x0))
  }
  return(x0)
}
