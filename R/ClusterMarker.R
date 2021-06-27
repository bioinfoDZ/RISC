####################################################################################
#' Find Cluster Markers
####################################################################################
#' 
#' This is the basic function in RISC, it can identify the cluster markers by
#' comparing samples in the selected cluster to the samples of the rest clusters.
#' Therefore, it is possible one gene labeled as a marker in more than one clusters.
#' Two methods are employed in RISC, one is based on Negative Binomial model while 
#' the other using QuasiPoisson model.
#' 
#' @rdname Cluster-Marker
#' @param object RISC object: a framework dataset.
#' @param cluster Select the cluster that we want to detect cluster marker genes.
#' @param positive Whether only output the cluster markers with positive logFC.
#' @param frac A fraction cutoff, the cluster marekr genes expressed at least a 
#' cutoff fraction of the cluster cells.
#' @param logFC The cutoff of log Fold-change for differentially expressed marker 
#' genes.
#' @param Padj The cutoff of the adjusted P-value.
#' @param latent.factor The latent factor from coldata, which represents number 
#' values or factors, and only one latent factor can be inputed.
#' @param method Which method is used to identify cluster markers, two options: 'NB' 
#' for Negative Binomial model and 'QP' for QuasiPoisson model.
#' @param ncore The multiple cores for parallel calculating.
#' @param cell.lim The minimum cells for each cluster to calculate marker genes.
#' @importFrom MASS glm.nb
#' @importFrom pbmcapply pbmclapply
#' @importFrom stats var
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
#' marker1 = scMarker(obj0, cluster = 1, method = 'QP', cell.lim = 3)

scMarker <- function(
  object, 
  cluster = 1, 
  positive = TRUE, 
  frac = 0.25, 
  logFC = 0.25, 
  Padj = 0.05, 
  latent.factor = NULL, 
  method = 'NB', 
  ncore = 1, 
  cell.lim = 10
  ) {
  
  # count = as.matrix(object@assay$logcount)
  count = object@assay$logcount
  coldata = as.data.frame(object@coldata)
  latent.factor0 = intersect(colnames(coldata), latent.factor[1])
  ncore = as.integer(ncore)
  Padj = as.numeric(Padj)
  logFC = as.numeric(logFC)
  cell.lim = as.integer(cell.lim)
  
  if(length(coldata$Cluster) == 0){
    stop('Please cluster cells first.')
  } else if(length(latent.factor0) > 0){
    
    cell1 = rownames(coldata)[coldata$Cluster %in% cluster]
    cell2 = rownames(coldata)[!coldata$Cluster %in% cluster]
    
    if(length(cell1) < cell.lim | length(cell2) < cell.lim){
      stop("The number of cells in either control or sample cluster is too low")
    } else {
      cell1 = cell1
      cell2 = cell2
    }
    
    cell0 = intersect(colnames(count), c(cell1, cell2))
    count1 = count[,colnames(count) %in% cell1]
    count2 = count[,colnames(count) %in% cell2]
    count0 = cbind(count2, count1)
    keep = Matrix::rowSums(count0 > 0) > floor(min(length(cell1), length(cell2)) * frac)
    count0 = count0[keep,]
    
    coldata0 = data.frame(coldata[,latent.factor0])
    coldata0 = coldata0[rownames(coldata0) %in% cell0,]
    group = factor(c(rep('ctrl', length(cell2)), rep('sam', length(cell1))), levels = c('ctrl', 'sam'))
    ave0 = Matrix::rowMeans(count0[,colnames(count0) %in% cell1])
    
    if(method == 'NB'){
      
      out = pbmcmapply(FUN = function(x){NB.latent.detect(count0[x,], coldata0, group)}, 1L:nrow(count0), mc.cores = ncore)
      # out = apply(count0, 1, FUN = function(county){glm.nb(formula = formula0, data = data.frame(coldata0, county, group))})
      
    } else {
      
      out = pbmcmapply(FUN = function(x){QP.latent.detect(count0[x,], coldata0, group)}, 1L:nrow(count0), mc.cores = ncore)
      # out = apply(count0, 1, FUN = function(county){glm(formula = formula0, data = data.frame(coldata0, county, group), family = quasipoisson)})
      
    }
    
    result = t(as.matrix(out))
    colnames(result) = c('logFC', 'Pvalue')
    result = data.frame(Symbol = rownames(count0), Avelog2 = ave0, result)
    result$Padj = p.adjust(result$Pvalue, method = 'BH')
    result = result[order(result$Padj, decreasing = FALSE),]
    
  } else {
    
    cell1 = rownames(coldata)[coldata$Cluster %in% cluster]
    cell2 = rownames(coldata)[!coldata$Cluster %in% cluster]
    
    if(length(cell1) < cell.lim | length(cell2) < cell.lim){
      stop("The number of cells in either control or sample cluster is too low")
    } else {
      cell1 = cell1
      cell2 = cell2
    }
    
    count1 = count[,colnames(count) %in% cell1]
    count2 = count[,colnames(count) %in% cell2]
    count0 = cbind(count2, count1)
    keep = Matrix::rowSums(count0 > 0) > floor(min(length(cell1), length(cell2)) * frac)
    count0 = count0[keep,]
    
    group = factor(c(rep('ctrl', length(cell2)), rep('sam', length(cell1))), levels = c('ctrl', 'sam'))
    ave0 = Matrix::rowMeans(count0[,colnames(count0) %in% cell1])
    
    if(method == 'NB'){
      
      out = pbmcmapply(FUN = function(x){NB.detect(count0[x,], group)}, 1L:nrow(count0), mc.cores = ncore)
      # out = apply(count0, 1, FUN = function(county){glm.nb(county ~ group)})
      
    } else {
      
      out = pbmcmapply(FUN = function(x){QP.detect(count0[x,], group)}, 1L:nrow(count0), mc.cores = ncore)
      # out = apply(count0, 1, FUN = function(county){glm(county ~ group, family = quasipoisson)})
      
    }
    
    result = t(as.matrix(out))
    colnames(result) = c('logFC', 'Pvalue')
    result = data.frame(Symbol = rownames(count0), Avelog2 = ave0, result)
    result$Padj = p.adjust(result$Pvalue, method = 'BH')
    result = result[order(result$Pvalue, decreasing = FALSE),]
    
  }
  
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
    result0 = result0[result0$logFC > logFC,]
  } else {
    result0 = result0[abs(result0$logFC) > logFC,]
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
#' @rdname All-Cluster-Marker
#' @param object RISC object: a framework dataset.
#' @param positive Whether only output the cluster markers with positive logFC.
#' @param frac A fraction cutoff, the cluster marker genes expressed at least a 
#' cutoff fraction of the cluster cells.
#' @param logFC The cutoff of log Fold-change for differentially expressed marker 
#' genes.
#' @param Padj The cutoff of the adjusted P-value.
#' @param latent.factor The latent factor from coldata, which represents number 
#' values or factors, and only one latent factor can be inputed.
#' @param min.cells The threshold for the cell number of valid clusters.
#' @param method Which method is used to identify cluster markers, two options: 
#' 'NB' for Negative Binomial model and 'QP' for QuasiPoisson model.
#' @param ncore The multiple cores for parallel calculating.
#' @importFrom MASS glm.nb
#' @name AllMarker
#' @export

AllMarker <- function(
  object, 
  positive = TRUE, 
  frac = 0.25, 
  logFC = 0.5, 
  Padj = 0.05, 
  latent.factor = NULL, 
  min.cells = 25L, 
  method = 'NB', 
  ncore = 1
  ) {
  
  if(length(object@cluster) == 0) {
    stop('Please cluster cells first')
  } else {
    
    Padj = as.numeric(Padj)
    logFC = as.numeric(logFC)
    ncore = as.numeric(ncore)
    coldata = as.data.frame(object@coldata)
    cluster0 = data.frame(table(coldata$Cluster))
    colnames(cluster0) = c('Cluster', 'No')
    keep = (cluster0$No >= min.cells)
    cluster0 = cluster0[keep,]
    result = data.frame(Symbol = 'Gene', Avelog2 = 0, logFC = 0, Pvalue = 1, Padj = 1, Cluster = 0)
    
    i = 1L
    while(i <= nrow(cluster0)) {
      
      clusteri = levels(as.factor(cluster0$Cluster))[i]
      cell1 = rownames(coldata)[coldata$Cluster %in% clusteri]
      cell2 = rownames(coldata)[!coldata$Cluster %in% clusteri]
      
      if(length(cell1) < min.cells | length(cell2) < min.cells){
        
        result = result
        i = i + 1L
        
      } else {
        
        resulti = scMarker(object = object, cluster = clusteri, positive = positive, frac = frac, logFC = logFC, Padj = Padj, latent.factor = latent.factor, method = method, ncore = ncore, cell.lim = min.cells)
        
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
#' @rdname Cluster-DEGs
#' @param object RISC object: a framework dataset.
#' @param cluster.ctrl Select the clusters that we want to detect DEGs, it/they as 
#' the control.
#' @param cluster.sam Select the clusters that we want to detect DEGs, it/they as 
#' the sample.
#' @param cell.ctrl Select the cells as the reference cells for detecting DEGs.
#' @param cell.sam Select the cells as the sample cells for detecting DEGs.
#' @param frac A fraction cutoff, the cluster marker genes expressed at least a 
#' cutoff fraction of the cluster cells.
#' @param logFC The cutoff of log Fold-change for differentially expressed marker 
#' genes.
#' @param Padj The cutoff of the adjusted P-value. If Padj is NULL, use p-value 
#' < 0.05 as a threshold. Set Padj as 1, without any cutoff.
#' @param latent.factor The latent factor from coldata, which represents number 
#' values or factors, and only one latent factor can be inputed.
#' @param method Which method is used to identify cluster markers, two options: 
#' 'NB' for Negative Binomial model and 'QP' for QuasiPoisson model.
#' @param ncore The multiple cores for parallel calculating.
#' @param cell.lim The minimum cells for each cluster to calculate marker genes.
#' @importFrom MASS glm.nb
#' @importFrom pbmcapply pbmclapply
#' @references Paternoster et al., Criminology (1997)
#' @references Berk et al., Journal of Quantitative Criminology (2008)
#' @references Liu et al., Nature Biotech. (2021)
#' @name scDEG
#' @export
#' @examples 
#' # RISC object
#' obj0 = raw.mat[[3]]
#' obj0 = scPCA(obj0, npc = 10)
#' obj0 = scUMAP(obj0, npc = 3)
#' obj0 = scCluster(obj0, slot = "cell.umap", k = 3, method = 'density')
#' DimPlot(obj0, slot = "cell.umap", colFactor = 'Cluster', size = 2)
#' DEG0 = scDEG(obj0, cluster.ctrl = 1, cluster.sam = 3, 
#'        cell.lim = 3, method = 'QP')

scDEG <- function(
  object, 
  cluster.ctrl = NULL, 
  cluster.sam = NULL, 
  cell.ctrl = NULL, 
  cell.sam = NULL, 
  frac = 0.10, 
  logFC = 0.25, 
  Padj = 0.05, 
  latent.factor = NULL, 
  method = 'NB', 
  ncore = 1, 
  cell.lim = 10
  ) {
  
  # count = as.matrix(object@assay$logcount)
  count = object@assay$logcount
  coldata = as.data.frame(object@coldata)
  latent.factor0 = intersect(colnames(coldata), latent.factor)
  ncore = as.numeric(ncore)
  logFC = as.numeric(logFC)
  cell.lim = as.numeric(cell.lim)
  
  if(length(object@cluster) == 0 & length(cell.ctrl) == 0) {
    stop('Please cluster cells first or input cell.ctrl and cell.sam.')
  } else if(length(latent.factor0) > 0) {
      
    if(length(cell.ctrl) > 0 & length(cell.sam) > 0){
      cell1 = cell.ctrl
      cell2 = cell.sam
    } else if(length(object@cluster) > 0 & length(cluster.ctrl) > 0 & length(cluster.sam) > 0){
      cluster0 = object@cluster
      cell1 = names(cluster0)[cluster0 %in% cluster.ctrl]
      cell2 = names(cluster0)[cluster0 %in% cluster.sam]
    }
    
    if(length(cell1) < cell.lim | length(cell2) < cell.lim){
      stop("The number of cells in either control or sample cluster is too low")
    } else {
      cell1 = cell1
      cell2 = cell2
    }
    
    cell0 = intersect(colnames(count), c(cell1, cell2))
    count1 = count[,colnames(count) %in% cell1]
    count2 = count[,colnames(count) %in% cell2]
    count0 = cbind(count1, count2)
    keep = Matrix::rowSums(count0 > 0) > floor(min(length(cell1), length(cell2)) * frac)
    count0 = count0[keep,]
    
    coldata0 = data.frame(coldata[,latent.factor0])
    coldata0 = coldata0[rownames(coldata0) %in% cell0,]
    group = factor(c(rep('ctrl', length(cell1)), rep('sam', length(cell2))), levels = c('ctrl', 'sam'))
    ave1 = Matrix::rowMeans(count0[,colnames(count0) %in% cell1])
    ave2 = Matrix::rowMeans(count0[,colnames(count0) %in% cell2])
    
    if(method == 'NB'){
      
      out = pbmcmapply(FUN = function(x){NB.latent.detect(count0[x,], coldata0, group)}, 1L:nrow(count0), mc.cores = ncore)
      # out = apply(count0, 1, FUN = function(county){glm.nb(formula = formula0, data = data.frame(coldata0, county, group))})
      
    } else {
      
      out = pbmcmapply(FUN = function(x){QP.latent.detect(count0[x,], coldata0, group)}, 1L:nrow(count0), mc.cores = ncore)
      # out = apply(count0, 1, FUN = function(county){glm(formula = formula0, data = data.frame(coldata0, county, group), family = quasipoisson)})
      
    }
    
    result = t(as.matrix(out))
    colnames(result) = c('logFC', 'Pvalue')
    result = data.frame(Symbol = rownames(count0), Avelog2ctrl = ave1, Avelog2sam = ave2, result)
    result$Padj = p.adjust(result$Pvalue, method = 'BH')
    result = result[order(result$Padj, decreasing = FALSE),]
    
  } else {
    
    if(length(cell.ctrl) > 0 & length(cell.sam) > 0){
      cell1 = cell.ctrl
      cell2 = cell.sam
    } else if(length(object@cluster) > 0 & length(cluster.ctrl) > 0 & length(cluster.sam) > 0){
      cluster0 = object@cluster
      cell1 = names(cluster0)[cluster0 %in% cluster.ctrl]
      cell2 = names(cluster0)[cluster0 %in% cluster.sam]
    }
    
    if(length(cell1) < cell.lim | length(cell2) < cell.lim){
      stop("The number of either control or sample cells is too low")
    } else {
      cell1 = cell1
      cell2 = cell2
    }
    
    count1 = count[,colnames(count) %in% cell1]
    count2 = count[,colnames(count) %in% cell2]
    count0 = cbind(count1, count2)
    keep = Matrix::rowSums(count0 > 0) > floor(min(length(cell1), length(cell2)) * frac)
    count0 = count0[keep,]
    
    group = factor(c(rep('ctrl', length(cell1)), rep('sam', length(cell2))), levels = c('ctrl', 'sam'))
    ave1 = Matrix::rowMeans(count0[,colnames(count0) %in% cell1])
    ave2 = Matrix::rowMeans(count0[,colnames(count0) %in% cell2])
    
    if(method == 'NB'){
      
      out = pbmcmapply(FUN = function(x){NB.detect(count0[x,], group)}, 1L:nrow(count0), mc.cores = ncore)
      # out = apply(count0, 1, FUN = function(county){glm.nb(county ~ group)})
      
    } else {
      
      out = pbmcmapply(FUN = function(x){QP.detect(count0[x,], group)}, 1L:nrow(count0), mc.cores = ncore)
      # out = apply(count0, 1, FUN = function(county){glm(county ~ group, family = quasipoisson)})
      
    }
    
    result = t(as.matrix(out))
    colnames(result) = c('logFC', 'Pvalue')
    result = data.frame(Symbol = rownames(count0), Avelog2ctrl = ave1, Avelog2sam = ave2, result)
    result$Padj = p.adjust(result$Pvalue, method = 'BH')
    result = result[order(result$Pvalue, decreasing = FALSE),]
    
  }
  
  if(is.null(Padj)){
    result0 = result[result$Pvalue < 0.05 & abs(result$logFC) >= logFC,]
  } else {
    Padj = as.numeric(Padj)
    
    if(Padj == 1) {
      result0 = result
    } else {
      result0 = result[result$Padj < Padj & abs(result$logFC) >= logFC,]
    }
    
  }
  
  return(result0)
  
}



####################################################################################
####################################################################################
NB.core <- function(x0, y0){
   l0 = suppressWarnings(glm.nb(x0 ~ y0))
   l1 = c(l0$coefficients[length(l0$coefficients)][[1]], summary(l0)$coef[nrow(summary(l0)$coef), 4])
   return(l1)
}

NB.latent.core <- function(x0, y0, z0){
  l0 = suppressWarnings(glm.nb(x0 ~ z0 + y0))
  l1 = c(l0$coefficients[length(l0$coefficients)][[1]], summary(l0)$coef[nrow(summary(l0)$coef), 4])
  return(l1)
}

NB.detect <- function(x0, y0){
  return(
    tryCatch(
      NB.core(x0, y0), 
      error = function(e){c(0, 1)}
    )
  )
}

NB.latent.detect <- function(x0, y0, z0){
  return(
    tryCatch(
      NB.latent.core(x0, y0, z0), 
      error = function(e){c(0, 1)}
    )
  )
}


## QuasiPoisson
QP.core <- function(x0, y0){
  l0 = glm(x0 ~ y0, family = quasipoisson)
  l1 = c(l0$coefficients[length(l0$coefficients)][[1]], summary(l0)$coef[nrow(summary(l0)$coef), 4])
  return(l1)
}

QP.latent.core <- function(x0, y0, z0){
  l0 = glm(x0 ~ z0 + y0, family = quasipoisson)
  l1 = c(l0$coefficients[length(l0$coefficients)][[1]], summary(l0)$coef[nrow(summary(l0)$coef), 4])
  return(l1)
}

QP.detect <- function(x0, y0){
  return(
    tryCatch(
      QP.core(x0, y0), 
      error = function(e){c(0, 1)}
    )
  )
}

QP.latent.detect <- function(x0, y0, z0){
  return(
    tryCatch(
      QP.latent.core(x0, y0), 
      error = function(e){c(0, 1)}
    )
  )
}


