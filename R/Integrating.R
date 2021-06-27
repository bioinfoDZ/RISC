####################################################################################
#' Integrating Multiple Datasets
####################################################################################
#' 
#' The "scMultiIntegrate" function can be used for data integration of multiple 
#' datasets, it is basically based on our new approach RPCI (reference principal 
#' components integration), which decomposes all the target datasets based on the 
#' reference data. The output of this function is RISC object, including the 
#' integrated eigenvectors and aligned gene expression values.
#' 
#' @rdname Multiple-Integrating
#' @param objects The list of multiple RISC objects: 
#' list{object1, object2, object3, ...}. The first set is the reference to generate 
#' gene-eigenvectors.
#' @param eigens The number of eigenvectors used for data integration.
#' @param add.Id Add a vector of Id to label different datasets, a character vector.
#' @param var.gene Define the variable genes manually. Here input a vector of gene 
#' names as variable genes
#' @param method Two methods are available, RPCI algorithm and SIMPLS algorithm, 
#' the default is RPCI. 
#' @param align The method for alignment of gene expression values: "Optimal" for 
#' alignment by experience, "Predict" for alignment by RPCI prediction, and "OLS" 
#' for alignment by the ordinary linear regression.
#' @param npc The number of the PCs returns from "scMultiIntegrate" function, 
#' they are usually used for the subsequent analyses, like cell embedding and 
#' cell clustering.
#' @param adjust Whether adjust the number of eigenvectors.
#' @param ncore The number of multiple cores for data integration.
#' @param do.fast With enough PC memory, we could do.fast, otherwise, 
#' let it as the default. Three options: AUTO (default), TRUE, and FALSE.
#' @param seed The random seed to keep consistent result.
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
#' @importFrom Matrix rowMeans colMeans rowSums colSums
#' @importFrom irlba irlba
#' @importFrom stats model.matrix contr.sum contrasts<-
#' @importFrom pbmcapply pbmcmapply pbmclapply
#' @references Liu et al., Nature Biotech. (2021)
#' @name scMultiIntegrate
#' @export
#' @examples
#' data("raw.mat")
#' obj1 = raw.mat[[3]]
#' obj2 = raw.mat[[4]]
#' obj0 = list(obj1, obj2)
#' var0 = intersect(obj1@vargene, obj2@vargene)
#' obj0 = scMultiIntegrate(obj0, eigens = 10, var.gene = var0, 
#'                        npc = 15, add.Id = c("Set-1", "Set-2"))
#' obj0 = scUMAP(obj0, npc = 10, use = "PLS")
#' DimPlot(obj0, slot = "cell.umap", colFactor = "Set", size = 2)
#' DimPlot(obj0, slot = "cell.umap", colFactor = "Group", size = 2, label = TRUE)

scMultiIntegrate <- function(
  objects, 
  eigens = 10, 
  add.Id = NULL, 
  var.gene = NULL, 
  method = 'RPCI', 
  align = 'OLS', 
  npc = 50, 
  adjust = TRUE, 
  ncore = 1, 
  do.fast = 'AUTO', 
  seed = 123
  ) {
  
  ncore = as.numeric(ncore)
  registerDoParallel(ncore)
  nset = nrow(summary(objects))
  
  if(is.null(add.Id)){
    names(objects) = paste0('Set', 1L:nset)
    add.Id = names(objects)
  } else {
    names(objects) = add.Id
  }
  
  if(is.null(objects)){
    stop('Please keep objects valid')
  } else {
    
    # prepare individual datasets
    i = j = l = 1L
    count0 = logcount0 = Var0 = gene0 = type0 = coldata0 = rowdata0 = list()
    Group = design0 = coef0 = vector()
    
    while(i <= nset){
      
      object = objects[[i]]
      object@coldata$scUMI = Matrix::colSums(object@assay$count)
      Var0[[i]] = object@vargene
      count0[[i]] = object@assay$count
      logcount0[[i]] = object@assay$logcount
      coldata0[[i]] = object@coldata
      rowdata0[[i]] = object@rowdata
      gene0[[i]] = rownames(object@assay$logcount)
      type0[[i]] = colnames(object@coldata)
      i = i + 1L
      
    }
    
    names(count0) = names(logcount0) = names(Var0) = names(coldata0) = names(rowdata0) = names(objects)
    type0 = Reduce(intersect, type0)
    rm(objects, object)
    gene0 = Reduce(intersect, gene0)
    if(isTRUE(adjust)){
      eigen0 = as.integer(eigens) - ifelse((as.integer(eigens)/10) > 1, ceiling(as.integer(eigens)/10), 0)
    } else {
      eigen0 = as.integer(eigens)
    }
    
    # select highly variable genes
    if(is.null(var.gene)) {
      var.gene = Reduce(intersect, Var0)
      
      if(length(var.gene) >= 3000){
        var.gene = unique(unlist(Var0, use.names = FALSE))
        var.gene = var.gene[!is.na(var.gene)]
      } else {
        var.gene = gene0
      }
      
    } else {
      var.gene = var.gene
    }
    var.gene = intersect(var.gene, gene0)
    
    # rank datasets
    while(l <= length(add.Id)){
      countl = count0[[l]]
      logcountl = logcount0[[l]]
      coldatal = coldata0[[l]][,type0]
      coldatal$Set = add.Id[l]
      rowdatal = rowdata0[[l]]
      rownames(rowdatal) = rownames(logcountl) = rownames(countl)
      colnames(countl) = colnames(logcountl) = rownames(coldatal) = coldatal[,'Barcode'] = paste0(add.Id[l], "_", colnames(countl))
      count0[[l]] = countl[gene0,]
      logcount0[[l]] = logcountl[gene0,]
      coldata0[[l]] = coldatal
      rowdata0[[l]] = rowdatal[gene0,]
      rm(countl, logcountl, coldatal, rowdatal)
      l = l + 1L
    }
    
    
    ### core integration
    logcount1L = as.matrix(logcount0[[1L]])
    Var1L = logcount1L[var.gene,]
    seed = as.numeric(seed)
    set.seed(seed)
    
    if(method == 'RPCI'){
      
      # RPCI
      Var.pca1L = irlba(Var1L, nv = eigen0, center = TRUE)
      
      scale.beta = foreach(j = 1L:length(add.Id)) %dopar% {
        logcountl = as.matrix(logcount0[[j]])
        Varj = scale(logcountl, center = TRUE, scale = FALSE)
        Varj = Varj[var.gene,]
        PCA1Lu = multipleCpp(Var.pca1L$u, diag(Var.pca1L$d, nrow = eigen0))
        PCAjv = crossprodCpp(PCA1Lu, Varj) / (Var.pca1L$d)^2
        beta0 = multipleCpp(Var.pca1L$v, PCAjv)
        rownames(beta0) = colnames(Var1L)
        colnames(beta0) = colnames(Varj)
        return(beta0)
      }
      
      rm(Var.pca1L)
      
    } else {
      
      # SIMPLS 
      scale.beta = foreach(j = 1L:length(add.Id)) %dopar% {
        logcountl = as.matrix(logcount0[[j]])
        Varj = logcountl[var.gene,]
        Varj = scale(Varj, center = TRUE, scale = FALSE)
        beta0 = simplsCpp(Var1L, Varj, ncomp = eigen0)$B
        rownames(beta0) = colnames(Var1L)
        colnames(beta0) = colnames(Varj)
        return(beta0)
      }
      
    }
    
    scale.beta0 = do.call(cbind, scale.beta)
    beta.pca = irlba(scale.beta0, nv = npc)
    beta.pca0 = as.matrix(beta.pca$v)
    colnames(beta.pca0) = paste0('PC', 1L:npc)
    rownames(beta.pca0) = colnames(scale.beta0)
    rm(Var1L, beta.pca, scale.beta)
    
    
    ### logcount alignment
    if(align == 'Predict'){
      
      logcount1L = as.matrix(logcount0[[1L]])
      logcountp = matrix()
      
      # RPCI prediction
      logcountp = multipleCpp(logcount1L, scale.beta0)
      colnames(logcountp) = colnames(scale.beta0)
      rownames(logcountp) = rownames(logcount1L)
      logcountp = as(logcountp, 'dgCMatrix')
      logcount0 = do.call(cbind, logcount0)
      keep = (logcount0 == 0)
      logcountp[logcountp <= 0] = 0
      logcountp[keep] = 0
      logcount.integrate = as(logcountp, 'dgCMatrix')
      rm(logcount1L, logcount0, keep, logcountp)
      
    } else if(align == 'Optimal'){
      
      logcount1L = logcount0[[1L]]
      ratio1L = pbmcmapply(function(x){quantile(logcount1L[x,], probs = 0.975)[[1]]}, 1L:nrow(logcount1L), mc.cores = ncore)
      # ratio1L = apply(logcount1L, 1, FUN = function(x){quantile(x, probs = 0.995)[[1]]})
      rowmean1L = Matrix::rowMeans(logcount1L)
      logcountj = list()
      j = as.integer()
      
      # Scale by experience
      logcountj = foreach(j = 2L:length(add.Id)) %dopar% {
        
        logcountl = logcount0[[j]]
        ratioj = pbmcmapply(function(x){quantile(logcountl[x,], probs = 0.975)[[1]]}, 1L:nrow(logcountl), mc.cores = ncore)
        # ratioj = apply(logcountl, 1, FUN = function(x){quantile(x, probs = 0.995)[[1]]})
        ratiom = ratio1L / ratioj
        ratiom[is.na(ratiom)] = ratiom[is.infinite(ratiom)] = ratiom[ratiom == 0] = 1
        rowmean0 = Matrix::rowMeans(logcountl)
        names(rowmean0) = rownames(logcountl)
        keep = (rowmean0 > 0)
        ratio.correct = glm(log2(ratiom[keep] + 1) ~ log2(rowmean0[keep] + 1), family = 'quasipoisson')$fitted.values
        max.correct = 2^max(ratio.correct) - 1
        min.correct = 2^min(ratio.correct) - 1
        ratiom[ratiom > max.correct] = max.correct
        ratiom[ratiom < min.correct] = min.correct
        logcountl0 = as(logcountl * ratiom, 'dgCMatrix')
        return(logcountl0)
        
      }
      
      logcount.integrate = cbind(as(logcount1L, 'dgCMatrix'), do.call(cbind, logcountj))
      rm(logcount1L, logcountj)
      
    } else if(align == 'OLS') {
      
      coldata.OLS = do.call(rbind, coldata0)
      Group = factor(coldata.OLS$Set, levels = add.Id)
      contrasts(Group) = contr.sum(levels(Group))
      design0 = model.matrix(~Group)
      Group = design0[, -1, drop = F]
      logcount.integrate = do.call(cbind, logcount0)
      
      if(do.fast == 'AUTO') {
        
        if(length(logcount.integrate) < 1e+08){
          coef0 = apply(logcount.integrate, 1, FUN = function(x){lmCpp(design0, x)})
        } else {
          coef0 = pbmcmapply(FUN = function(x){lmCpp(design0, logcount.integrate[x,])}, 1L:nrow(logcount.integrate), mc.cores = ncore)
        }
        
      } else if (do.fast) {
        coef0 = apply(logcount.integrate, 1, FUN = function(x){lmCpp(design0, x)})
      } else {
        coef0 = pbmcmapply(FUN = function(x){lmCpp(design0, logcount.integrate[x,])}, 1L:nrow(logcount.integrate), mc.cores = ncore)
      }
      
      colnames(coef0) = rownames(logcount.integrate)
      coef0 = t(as.matrix(coef0))
      coef0 = coef0[, -1L, drop = FALSE]
      coef0[is.na(coef0)] = 0
      logcount.integrate = logcount.integrate - multipleCpp(coef0, t(Group))
      logcount.integrate[logcount.integrate <= 0] = 0
      count.OLS = do.call(cbind, count0)
      keep = (count.OLS == 0)
      logcount.integrate[keep] = 0
      logcount.integrate = as(logcount.integrate, 'dgCMatrix')
      rm(coldata.OLS, Group, design0, coef0, count.OLS, keep)
      
      
    } else {
      stop('Please input align, Optimal or Predict')
    }
    
    rm(Var0)
    
    
    ### combine data
    count.integrate = do.call(cbind, count0)
    count.integrate = as(count.integrate, 'dgCMatrix')
    rm(count0)
    
    coldata.integrate = do.call(rbind, coldata0)
    rownames(coldata.integrate) = coldata.integrate[,'Barcode']
    coldata.integrate$Set = as.factor(coldata.integrate$Set)
    rm(coldata0)
    
    rowdata.integrate = data.frame(Symbol = rownames(logcount.integrate), RNA = "Gene Expression", stringsAsFactors = F)
    rownames(rowdata.integrate) = rownames(logcount.integrate)
    
    # Output
    object = SingleCellData(assay = list(count = count.integrate, logcount = logcount.integrate), rowdata = rowdata.integrate, coldata = coldata.integrate)
    object@metadata[['normalise']] = 'Yes'
    object@metadata[['Beta']] = as.matrix(scale.beta0)
    object@metadata[['Integration']] = data.frame(eigens = eigen0, npc = npc, method = "RPCI", align = align, stringsAsFactors = FALSE)
    object@vargene = var.gene
    object@DimReduction[['cell.pls']] = as.matrix(beta.pca0)
    return(object)
    
  }
  
}



####################################################################################
#' Integrating Multiple Large Datasets
####################################################################################
#' 
#' The "scMultiIntegrate" function can be used for data integration of multiple 
#' datasets, it is basically based on our new algorithm: reference principal 
#' components integration (RPCI). RPCI decomposes all the target datasets based 
#' on the reference. The output of this function can be used for low dimension 
#' visualization.
#' 
#' @rdname Large-Integrating
#' @param objects The list of multiple RISC objects: 
#' list{object1, object2, object3, ...}.
#' @param K The number of PCs used for data integration, the optimal number of 
#' PCs is equal to the number of PCs used for dimension reduction of individual 
#' datasets.
#' @param npc The number of the PCs returns from "scMultiIntegrate" function, 
#' they are usually used for the subsequent analyses, like cell embedding.
#' @param ncore The number of multiple cores for data integration.
#' @param Add.Id Add a vector of Ids to label different datasets, a character vector.
#' @param var.gene Define the variable genes mannually. Here input a vector of gene 
#' names as variable genes
#' @param seed The random seed to keep consistent result.
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
#' @importFrom irlba irlba
#' @references Liu et al., Nature Biotech. (2021)
#' @name scLargeIntegrate


scLargeIntegrate <- function(
  objects, 
  K = 10, 
  npc = 50, 
  var.gene = NULL, 
  Add.Id = NULL, 
  ncore = 1, 
  seed = 123
) {
  
  ncore = as.numeric(ncore)
  registerDoParallel(ncore)
  nset = nrow(summary(objects))
  
  if(is.null(Add.Id)){
    Add.Id = names(objects) = paste0('Set', 1L:nset)
  } else {
    names(objects) = Add.Id
  }
  
  if(is.null(var.gene)){
    stop("Please input var.gene")
  } else {
    var.gene = as.vector(var.gene)
  }
  
  if(is.null(objects)){
    stop('Please keep objects valid')
  } else {
    
    # prepare individual datasets
    i = j = 1L
    Var0 = beta0 = list()
    logcounti = matrix()
    
    while(i <= nset){
      object = objects[[i]]
      objects[[i]] = 0
      logcounti = object@assay$logcount
      rm(object)
      var.gene = intersect(var.gene, rownames(logcounti))
      logcounti = logcounti[var.gene,]
      colnames(logcounti) = paste0(Add.Id[i], "_", colnames(logcounti))
      Var0[[i]] = logcounti
      i = i + 1L
    }
    rm(objects, logcounti)
    
    
    ### core integration
    Var1L = Var0[[1L]]
    seed = as.numeric(seed)
    set.seed(seed)
    
    Var1L = scale(as.matrix(Var1L), center = TRUE, scale = FALSE)
    Var.pca1L = irlba(Var1L, nv = K)
    
    scale.beta = foreach(j = 1L:nset) %dopar% {
      Varj = Var0[[j]]
      Varj = scale(as.matrix(Varj), center = TRUE, scale = FALSE)
      PCA1Lu = multipleCpp(Var.pca1L$u, diag(Var.pca1L$d, nrow = K))
      PCAjv = crossprodCpp(PCA1Lu, Varj) / (Var.pca1L$d)^2
      beta0 = multipleCpp(Var.pca1L$v, PCAjv)
      beta0 = as.matrix(beta0)
      rownames(beta0) = colnames(Var1L)
      colnames(beta0) = colnames(Varj)
      return(beta0)
    }
    
    scale.beta0 = do.call(cbind, scale.beta)
    set.seed(seed)
    beta.pca = irlba(scale.beta0, nv = npc)
    beta.pca0 = as.matrix(beta.pca$v)
    colnames(beta.pca0) = paste0('PC', 1L:npc)
    rownames(beta.pca0) = colnames(scale.beta0)
    rm(Var1L, Var.pca1L, beta.pca, scale.beta)
    beta0[['PLS']] = beta.pca0
    beta0[['Coef']] = scale.beta0
    
    return(beta0)
    
  }
  
}



####################################################################################
#' Integrating Multiple Datasets
####################################################################################
#' 
#' The "scMultiIntegrate" function can be used for data integration of multiple 
#' datasets, it is basically based on our new algorithm reference principal 
#' components integration (RPCI). RPCI decomposes all the target datasets based 
#' on the reference data. The output of this function aligns the gene expression 
#' of different datasets.
#' 
#' @rdname Align-Gene-Expression
#' @param objects The list of multiple RISC objects: 
#' list{object1, object2, object3, ...}. 
#' @param Align The method for alignment of gene expression values: "Optimal" for 
#' alignment by experience, "Predict" for alignment by RPCI prediction, and "OLS" 
#' for alignment by the ordinary linear regression.
#' @param ncore The number of multiple cores for data integration.
#' @param Add.Id Add a vector of Id to label different datasets, a character vector.
#' @param beta The coefficiency to generate predicted counts.
#' @param do.fast With enough PC memory, we could do.fast, otherwise, let it as the 
#' default. Three options: AUTO (default), TRUE, and FALSE.
#' @param seed The random seed to keep consistent result.
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
#' @importFrom Matrix rowMeans colMeans rowSums colSums
#' @importFrom irlba irlba
#' @importFrom stats model.matrix contr.sum contrasts<-
#' @importFrom pbmcapply pbmcmapply pbmclapply
#' @references Liu et al., Nature Biotech. (2021)
#' @name scAlignGene


scAlignGene <- function(
  objects, 
  Align = 'OLS', 
  ncore = 1, 
  Add.Id = NULL, 
  do.fast = 'AUTO', 
  beta = NULL, 
  seed = 123
) {
  
  ncore = as.numeric(ncore)
  registerDoParallel(ncore)
  nset = nrow(summary(objects))
  
  if(is.null(Add.Id)){
    Add.Id = names(objects) = paste0('Set', 1L:nset)
  } else {
    names(objects) = Add.Id
  }
  
  if(is.null(objects)){
    stop('Please keep objects valid')
  } else {
    
    # prepare individual datasets
    i = j = l = 1L
    logcount0 = coldata0 = gene0 = list()
    logcounti = matrix()
    
    while(i <= nset){
      object = objects[[i]]
      logcounti = object@assay$logcount
      coldatai = object@coldata
      rm(object)
      rownames(coldatai) = colnames(logcounti) = paste0(Add.Id[i], "_", colnames(logcounti))
      coldatai$Set = Add.Id[i]
      logcount0[[i]] = logcounti
      coldata0[[i]] = as.data.frame(coldatai)
      gene0[[i]] = rownames(logcounti)
      i = i + 1L
    }
    rm(objects, logcounti)
    logcountl = matrix()
    gene0 = Reduce(intersect, gene0)
    
    while(l <= nset){
      logcountl = logcount0[[l]]
      logcount0[[l]] = logcountl[gene0,]
      l = l + 1L
    }
    rm(logcountl)
    
    ### logcount alignment
    if(Align == 'Predict'){
      
      logcount1L = as.matrix(logcount0[[1L]])
      logcountp = matrix()
      
      if(is.null(scale.beta0)){
        stop('Input matrix coefficient')
      } else {
        scale.beta0 = as.matrix(scale.beta0)
      }
      
      # RPCI prediction
      logcountp = multipleCpp(logcount1L, scale.beta0)
      colnames(logcountp) = colnames(scale.beta0)
      rownames(logcountp) = rownames(logcount1L)
      # logcountp = cbind(logcount1L, logcountp)
      logcount0 = do.call(cbind, logcount0)
      keep = (as.matrix(logcount0) == 0)
      logcountp[logcountp < 0] = 0
      logcountp[keep] = 0
      logcount.integrate = as(logcountp, 'dgCMatrix')
      rm(logcount1L, logcount0, keep, logcountp)
      
    } else if(Align == 'Optimal'){
      
      logcount1L = logcount0[[1L]]
      ratio1L = pbmcmapply(function(x){quantile(logcount1L[x,], probs = 0.975)[[1]]}, 1L:nrow(logcount1L), mc.cores = ncore)
      # ratio1L = apply(logcount1L, 1, FUN = function(x){quantile(x, probs = 0.995)[[1]]})
      rowmean1L = rowMeans(as.matrix(logcount1L))
      logcountj = list()
      j = as.integer()
      
      # Scale by experience
      logcountj = foreach(j = 2L:nset) %dopar% {
        
        logcountl = logcount0[[j]]
        ratioj = pbmcmapply(function(x){quantile(logcountl[x,], probs = 0.975)[[1]]}, 1L:nrow(logcountl), mc.cores = ncore)
        # ratioj = apply(logcountl, 1, FUN = function(x){quantile(x, probs = 0.995)[[1]]})
        ratiom = ratio1L / ratioj
        ratiom[is.na(ratiom)] = ratiom[is.infinite(ratiom)] = ratiom[ratiom == 0] = 1
        rowmean0 = rowMeans(as.matrix(logcountl))
        names(rowmean0) = rownames(logcountl)
        keep = (rowmean0 > 0)
        ratio.correct = glm(log2(ratiom[keep] + 1) ~ log2(rowmean0[keep] + 1), family = 'quasipoisson')$fitted.values
        max.correct = 2^max(ratio.correct) - 1
        min.correct = 2^min(ratio.correct) - 1
        ratiom[ratiom > max.correct] = max.correct
        ratiom[ratiom < min.correct] = min.correct
        logcountl0 = as(logcountl * ratiom, 'dgCMatrix')
        return(logcountl0)
        
      }
      
      logcount0 = c(as(logcount1L, 'dgCMatrix'), logcountj)
      rm(logcount1L, logcountj)
      logcount.integrate = do.call(cbind, logcount0)
      
    } else if(Align == 'OLS') {
      
      coldata.OLS = do.call(rbind, coldata0)
      Group = factor(coldata.OLS$Set, levels = Add.Id)
      contrasts(Group) = contr.sum(levels(Group))
      design0 = model.matrix(~Group)
      Group = design0[, -1, drop = FALSE]
      logcount.integrate = do.call(cbind, logcount0)
      
      if(do.fast == 'AUTO') {
        
        if(length(logcount.integrate) < 1e+08){
          coef0 = apply(logcount.integrate, 1, FUN = function(x){lmCpp(design0, x)})
        } else {
          coef0 = pbmcmapply(FUN = function(x){lmCpp(design0, logcount.integrate[x,])}, 1L:nrow(logcount.integrate), mc.cores = ncore)
        }
        
      } else if (do.fast) {
        coef0 = apply(logcount.integrate, 1, FUN = function(x){lmCpp(design0, x)})
      } else {
        coef0 = pbmcmapply(FUN = function(x){lmCpp(design0, logcount.integrate[x,])}, 1L:nrow(logcount.integrate), mc.cores = ncore)
      }
      
      colnames(coef0) = rownames(logcount.integrate)
      coef0 = t(as.matrix(coef0))
      coef0 = coef0[, -1L, drop = FALSE]
      coef0[is.na(coef0)] = 0
      logcount.integrate = as.matrix(logcount.integrate) - multipleCpp(coef0, t(Group))
      logcount.integrate[logcount.integrate < 0] = 0
      logcount.integrate[keep] = 0
      logcount.integrate = as(logcount.integrate, 'dgCMatrix')
      rm(coldata.OLS, Group, design0, coef0, keep)
      
    } else {
      stop('Please input Align, Optimal or Predict')
    }
    
    return(logcount.integrate)
    
  }
  
}



####################################################################################
#' Integration Plot
####################################################################################
#' 
#' The "InPlot" function makes the plot to show how the PCs explain the variance 
#' for data integration. This plot helps the users to select the optimal reference 
#' and the PCs to perform data integration. 
#' 
#' @rdname InPlot
#' @param object A list of RISC objects.
#' @param var.gene The highly variable genes.
#' @param Colors The colors labeling for different data sets.
#' @param nPC The PCs will be calculated.
#' @param neighbor The nearest neighbors. 
#' @param method The method of cell clustering for individual datasets.
#' @param algorithm The algorithm for knn, the default is "kd_tree", all options: 
#' "kd_tree", "cover_tree", "CR", "brute".
#' @param ncore The number of multiple cores for testing.
#' @param minPC The minimal PCs for detecting cell clustering.
#' @param Std.cut The cutoff of standard deviation of the PCs.
#' @param bin The bin number for calculating cell clustering.
#' @importFrom stats ks.test reshape
#' @importFrom matrixStats rowMedians
#' @references Liu et al., Nature Biotech. (2021)
#' @name InPlot
#' @export

InPlot <- function(
  object = NULL, 
  var.gene = NULL, 
  Colors = NULL, 
  nPC = 20, 
  neighbor = 30, 
  method = "louvain", 
  algorithm = "kd_tree", 
  ncore = 1, 
  minPC = 11, 
  Std.cut = 0.95, 
  bin = 5
) {
  
  set.seed(123)
  if(is.null(objects)){
    stop('Please keep objects valid')
  } else {
    ncore = as.integer(ncore)
    registerDoParallel(ncore)
    npc = as.integer(nPC)
    neighbor = as.integer(neighbor)
    method = as.character(method)
    algorithm = as.character(algorithm)
    nset = dim(summary(object))[1]
    minpc = as.integer(minPC)
    Std.cut = as.numeric(Std.cut)
    bin = as.integer(bin)
    
    if(minpc >= npc){
      minpc = npc
    } else {
      minpc = minpc
    }
  }
  
  # choose the reference by Cells PCs
  gene0 = pc.gene = var0 = PC0 = list()
  vari = PCi = PC1 = PC2 = k0 = ki = ka = kb = kc = Std.Num = Std0 = Cluster = Num = PCs = Group = Group1 = Group2 = mean0 = NULL
  i = j = 1L
  
  while(i <= nset){
    gene0[[i]] = object[[i]]@rowdata$Symbol
    var0[[i]] = object[[i]]@vargene
    i = i + 1L
  }
  
  gene0 = Reduce(intersect, gene0)
  var0 = Reduce(intersect, var0)
  
  if(is.null(var.gene)) {
    var0 = intersect(gene0, var0)
  } else {
    var0 = intersect(gene0, var.gene)
  }
  
  PC0 = foreach(i = 1L:nset) %dopar% {
    vari = object[[i]]@assay$logcount[var0,]
    PCi = irlba(vari, nv = npc, center = T)
    return(PCi)
  }
  
  PC1 = foreach(i = 1L:nset) %dopar% {
    ka = apply(PC0[[i]]$u, 2, FUN = function(x){ks.test(x, "pnorm", mean=mean(x), sd=sd(x))$statistic})
    kb = glm(ka ~ c(1:npc), family = quasipoisson)$deviance
    kc = mean(ka)
    return(list(ka, kb, kc))
  }
  ka = unlist(lapply(PC1, `[[`, 2))
  ka = (1-ka)*(1/max(1-ka))
  kc = unlist(lapply(PC1, `[[`, 3))
  
  PC2 = foreach(i = 1L:nset) %dopar% {
    PCi = PC0[[i]]$d^2 / sum(PC0[[i]]$d^2)
    PCi = cumsum(PCi)
    Std0 = c(1:npc)[PCi > Std.cut][1]
    return(list(PCi, Std0))
  }
  kb = sapply(1L:nset, FUN = function(x){glm(PC2[[x]][[1]] ~ c(1L:npc), family = poisson)$coefficients[2]})
  #  kb = kb*(0.5/min(kb))
  #  kb = 1 - (kb - max(kb)) / (min(kb) - max(kb))
  kb = 1 - (1/npc - kb) / (1/npc)
  Std.Num = unlist(lapply(PC2, `[[`, 2))
  
  bin0 = as.integer(seq(minpc, npc, length.out = bin))
  if(method == "louvain") {
    
    clust0 = foreach(i = rep(1L:nset, length(bin0)), j = rep(bin0, each = nset)) %dopar% {
      k0 = get.knn(PC0[[i]]$v[,1L:j], k = neighbor, algorithm = algorithm)
      ki = data.frame(NodStar = rep(1L:nrow(PC0[[i]]$v[,1L:j]), neighbor), NodEnd = as.vector(k0$nn.index), stringsAsFactors = FALSE)
      ki = graph_from_data_frame(ki, directed = FALSE)
      ki = simplify(ki)
      ki = cluster_louvain(ki, weights = 1/(1 + as.vector(k0$nn.dist)))
      k0 = length(unique(ki$membership))
      return(as.integer(k0))
    }
    
  } else {stop('A new method later')}
  
  PC3 = data.frame(PC1 = unlist(lapply(PC1, `[[`, 1)), PC2 = unlist(lapply(PC2, `[[`, 1)), Group = rep(paste0("Set-", 1L:nset), each=npc), Num = rep(1L:npc, nset), stringsAsFactors = F)
  PC3$Group = factor(PC3$Group, levels = paste0("Set-", 1L:nset))
  PC3$mean0 = rep(round(kc, 2), each = npc)
  PC3$Group1 = factor(PC3$Group, labels = paste0(levels(PC3$Group), ": ", round(ka, 2)))
  PC3$Group2 = factor(PC3$Group, labels = paste0(levels(PC3$Group), ": ", round(as.numeric(kb), 2)))
  PC4 = data.frame(Group = paste0("Set-", 1L:nset), Cluster = matrix(unlist(clust0), nset, ), stringsAsFactors = F)
  colnames(PC4) = c("Set", paste0("PC-", bin0))
  cmean = rowMedians(matrix(unlist(clust0), nset, ))
  cmean = round(cmean / max(cmean), 2)
  PC4 = reshape(PC4, idvar = "Group", varying = list(2:ncol(PC4)), v.names = "Cluster", direction = "long")
  colnames(PC4) = c("Set", "PCs", "Cluster", "Group")
  PC4$Group = factor(PC4$Set, levels = paste0("Set-", 1L:nset), labels = paste0("Set-", 1L:nset, ": ", cmean))
  PC4$PCs = factor(PC4$PCs, levels = 1:length(bin0), labels = paste0("PC-", bin0))
  PC4$Cluster = as.integer(PC4$Cluster)
  
  if(length(levels(PC3$Group)) <= 35){
    
    if(is.null(Colors)){
      Cols0 = Color0(length(levels(PC3$Group)))
    } else {
      Cols0 = Colors
    }
    
    plot1 = ggplot(PC3, aes(x = Num, y = PC2, g = Group)) + 
      geom_point(aes(fill = Group2), size = 3, shape = 21, alpha = 0.8) + 
      geom_smooth(aes(color = Group), se=FALSE, method="loess", formula = "y~x", span=1, linetype="twodash", show.legend = FALSE) + 
      scale_fill_manual(values = Cols0) + 
      scale_color_manual(values = Cols0) + 
      labs(x = "PC Number", y = "Score (Stv by PCs)", fill = "Scores") + 
      theme_bw(base_line_size = 0, base_size = 12) + 
      theme(axis.title = element_text(size = 12, face = "bold")) + 
      guides(fill = guide_legend(override.aes = list(size = 4)))
    
    plot2 = ggplot(PC3, aes(x = Group, y = PC1)) + 
      geom_violin(aes(color = Group), show.legend = FALSE) + 
      geom_boxplot(aes(fill = Group1), width = 0.1, show.legend = FALSE) + 
      geom_jitter(aes(fill = Group1), size = 3, pch = 21, alpha = 0.8, width = 0.25) + 
      geom_label(aes(x = Group, y = mean0, label = mean0), color = "firebrick", nudge_x = 0.25) + 
      scale_fill_manual(values = Cols0) + 
      scale_color_manual(values = Cols0) + 
      labs(x = "", y = "Score (Kolmogorov-Smirnov)", fill = "Scores") + 
      theme_bw(base_line_size = 0, base_size = 12) + 
      theme(axis.title = element_text(size = 12, face = "bold")) + 
      guides(fill = guide_legend(override.aes = list(size = 4)))
    
    plot3 = ggplot(PC4, aes(x = PCs, y = Cluster, group = Group)) + 
      geom_bar(aes(fill = Group), stat = "identity", width = 0.8, position = "dodge") + 
      scale_fill_manual(values = Cols0) + 
      labs(x = "PCs", y = "Cluster Num.", fill = "Scores") + 
      theme_bw(base_line_size = 0, base_size = 12) + 
      theme(axis.title = element_text(size = 12, face = "bold")) + 
      guides(fill = guide_legend(override.aes = list(size = 4)))
    
  } else {
    
    plot1 = ggplot(PC3, aes(x = Num, y = PC2, g = Group)) + 
      geom_point(aes(fill = Group2), size = 3, shape = 21, alpha = 0.8) + 
      geom_smooth(aes(color = Group), se=FALSE, method="loess", formula = "y~x", span=1, linetype="twodash", show.legend = FALSE) + 
      labs(x = "PC Number", y = "Score (Stv by PCs)", fill = "Scores") + 
      theme_bw(base_line_size = 0, base_size = 12) + 
      theme(axis.title = element_text(size = 12, face = "bold")) + 
      guides(fill = guide_legend(override.aes = list(size = 4)))
    
    plot2 = ggplot(PC3, aes(x = Group, y = PC1)) + 
      geom_violin(aes(color = Group), show.legend = FALSE) + 
      geom_boxplot(aes(fill = Group1), width = 0.1, show.legend = FALSE) + 
      geom_jitter(aes(fill = Group1), size = 3, pch = 21, alpha = 0.8, width = 0.25) + 
      geom_label(aes(x = Group, y = mean0, label = mean0), color = "firebrick", nudge_x = 0.25) + 
      labs(x = "", y = "Score (Kolmogorov-Smirnov)", fill = "Scores") + 
      theme_bw(base_line_size = 0, base_size = 12) + 
      theme(axis.title = element_text(size = 12, face = "bold")) + 
      guides(fill = guide_legend(override.aes = list(size = 4)))
    
    plot3 = ggplot(PC4, aes(x = PCs, y = Cluster, group = Group)) + 
      geom_bar(aes(fill = Group), stat = "identity", width = 0.8, position = "dodge") + 
      labs(x = "PCs", y = "Cluster Num.", fill = "Scores") + 
      theme_bw(base_line_size = 0, base_size = 12) + 
      theme(axis.title = element_text(size = 12, face = "bold")) + 
      guides(fill = guide_legend(override.aes = list(size = 4)))
    
  }
  
  # print(paste0("The Number of PCs for ", Std.cut*100, "% Gene Expression Variance:"))
  # cat(paste0(levels(PC3$Group), ": ", as.integer(Std.Num)), sep = "\n")
  return(grid.arrange(plot3, plot1, plot2, ncol = 1))
  
}



####################################################################################
#' Integration Algorithm SIMPLS
####################################################################################
#' 
#' The partial least square (PLS) with SIMPLS algorithm is an extension of the 
#' multiple linear regression model and considered as bilinear factor models. 
#' Instead of embedding the reference and target matrices into a hyperplane 
#' of maximum variance, the PLS utilizes a linear regression to project the 
#' reference and target matrices into a new place. The SIMPLS algorithm provides 
#' the regularization procedure for PLS. The matrices need to be centered before 
#' SIMPLS integraton. 
#' 
#' @rdname Integration-Algorithm-SIMPLS
#' @param X The reference matrix, row for genes and column for cells.
#' @param Y The target matrix, row for genes and column for cells.
#' @param npcs The number of the PCs used for data integration.
#' @param seed The random seed to keep consistent result.
#' @importFrom Matrix colMeans
#' @importFrom matrixStats colSds
#' @references De-Jong et al. (1993)
#' @name SIMPLS
#' @export
#' @examples
#' ## SIMPLS with two matrices
#' mat0 = as.matrix(raw.mat[[1]])
#' coldata0 = as.data.frame(raw.mat[[2]])
#' coldata1 = coldata0[coldata0$Batch0 == 'Batch1',]
#' coldata2 = coldata0[coldata0$Batch0 == 'Batch4',]
#' mat1 = mat0[,rownames(coldata1)]
#' mat2 = mat0[,rownames(coldata2)]
#' SIM0 = SIMPLS(mat1, mat2, npcs = 4)

SIMPLS <- function(
  X, 
  Y, 
  npcs = 10, 
  seed = 123
  ) {
  
  row0 = nrow(X)
  xcol = ncol(X)
  ycol = ncol(Y)
  V = R = matrix(0, xcol, npcs)
  tQ = matrix(0, npcs, ycol)
  B = matrix(0, xcol, ycol)
  S = crossprodCpp(X, Y)
  
  i = 1L
  while(i <= npcs){
    
    if(ycol < xcol){
      set.seed(seed)
      q.a = eigen(selfcrossprodCpp(S), symmetric = TRUE)$vectors[,1]
    } else {
      q.a = c(crossprodCpp(S, eigen(selftcrossprodCpp(S), symmetric = TRUE)$vectors[,1]))
      q.a = q.a / sqrt(c(selfcrossprodCpp(q.a)))
    }
    
    r.a = multipleCpp(S, q.a)
    t.a = multipleCpp(X, r.a)
    t.a = t.a - mean(t.a)
    tnorm = sqrt(c(selfcrossprodCpp(t.a)))
    t.a = t.a / tnorm
    r.a = r.a / tnorm
    p.a = crossprodCpp(X, t.a)
    q.a = crossprodCpp(Y, t.a)
    v.a = p.a
    v.a = v.a - multipleCpp(V, crossprodCpp(V, p.a))
    v.a = v.a /sqrt(c(selfcrossprodCpp(v.a)))
    S = S - multipleCpp(v.a, crossprodCpp(v.a, S))
    R[,i] = r.a
    tQ[i,] = q.a
    V[,i] = v.a
    B = R[, 1:i, drop = FALSE] %*% tQ[1:i,, drop = FALSE]
    i = i + 1L
    
  }
  
  return(B)
  
}



####################################################################################
#' Alignment of gene expression values
####################################################################################
#' 
#' This funciton is not used in data integration but for adjusting the results, 
#' usually the users do not need this.
#' 
#' @rdname MSC
#' @param X The reference matrix, row for genes and column for cells.
#' @param Y The target matrix, row for genes and column for cells.
#' @references Mevik et al., JSS (2007)
#' @name MSC
#' @export
#' @examples
#' ## Calculate MSC with the inputs of two matrices
#' mat0 = as.matrix(raw.mat[[1]])
#' coldata0 = as.data.frame(raw.mat[[2]])
#' coldata1 = coldata0[coldata0$Batch0 == 'Batch1',]
#' coldata2 = coldata0[coldata0$Batch0 == 'Batch4',]
#' mat1 = mat0[,rownames(coldata1)]
#' mat2 = mat0[,rownames(coldata2)]
#' msc0 = MSC(mat1, mat2)

MSC <- function(X, Y){
  
  x = t(as.matrix(X))
  y = t(as.matrix(Y))
  ref = colMeans(x)
  z = cbind(1, ref)
  beta = t(solve(selfcrossprodCpp(z), t(multipleCpp(y, z))))
  res = (y - beta[,1]) / beta[,2]
  res = t(as.matrix(res))
  return(res)
  
}



####################################################################################
####################################################################################
correct.zero <- function(x){
  l1 = length(x[x > 0])
  l2 = length(x)
  ratio = 1 - l1/l2
  return(ratio)
}

SPLSDV <- function(Z, Lamda = Lamda){
  zrow = nrow(Z)
  zcol = ncol(Z)
  Znorm = median(abs(Z))
  Z = Z / Znorm
  M = selftcrossprodCpp(Z)
  
  dis = 10
  ctl = matrix(10, zrow, 1)
  ctl0 = ctl
  i = 1L
  while(dis > 1e-4 & i <= 300L){
    set.seed(123)
    mcsvd = svd(multipleCpp(M, ctl))
    sam = tcrossprodCpp(mcsvd$u, mcsvd$v)
    Ms = multipleCpp(M, sam)
    M0 = matrix(0, length(Ms), 1)
    index = abs(Ms) - Lamda*max(abs(Ms))
    M0[index >= 0] = index[index >= 0] * (sign(Ms))[index >= 0]
    ctl = M0
    dis = max(abs(ctl - ctl0))
    ctl0 = ctl
    i = i + 1L
  }
  return(ctl)
}

SPLS <- function(X, Y, K = K, Lamda = 0.5, scale = TRUE){
  x = as.matrix(X)
  y = as.matrix(Y)
  row0 = nrow(x)
  xcol = ncol(x)
  ycol = ncol(y)
  IP = c(1:xcol)
  weight = matrix(1, 1, row0)
  
  xsd = colSds(x)
  if(scale){
    x = t(t(x)/colSds(x))
  } else {
    x = x
  }
  ysd = rep(1, ycol)
  
  beta = matrix(0, xcol, ycol)
  x0  = x
  y0 = y
  i = 1L
  while(i <= K){
    Z = crossprodCpp(x0, y0)
    tag = SPLSDV(Z, Lamda = Lamda)
    sel = unique(IP[tag != 0 | beta[,1] != 0])
    xsel = x[, sel, drop = FALSE]
    coefxy = SIMPLS(xsel, y, npcs = min(i, length(sel)))
    beta[sel,] = matrix(coefxy, length(sel), ycol)
    y0 = y - multipleCpp(x, beta)
    i = i + 1L
  }
  return(beta)
}


