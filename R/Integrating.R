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
#' @param align The method for alignment of gene expression values: "Optimal" for 
#' alignment by experience, "Predict" for alignment by RPCI prediction, and "OLS" 
#' for alignment by the ordinary linear regression.
#' @param npc The number of the PCs returns from "scMultiIntegrate" function, 
#' they are usually used for the subsequent analyses, like cell embedding and 
#' cell clustering.
#' @param adjust Whether adjust the number of eigenvectors.
#' @param ncore The number of multiple cores for data integration.
#' @param seed The random seed to keep consistent result.
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
#' @importFrom Matrix rowMeans colMeans rowSums colSums spMatrix sparseMatrix
#' @importFrom irlba irlba
#' @importFrom stats embed model.matrix contr.sum contrasts<-
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom pbapply pblapply pboptions
#' @references Liu et al., Nature Biotech. (2021)
#' @name scMultiIntegrate
#' @export
#' @examples
#' obj1 = raw.mat[[3]]
#' obj2 = raw.mat[[4]]
#' obj0 = list(obj1, obj2)
#' var0 = intersect(obj1@vargene, obj2@vargene)
#' obj0 = scMultiIntegrate(obj0, eigens = 8, var.gene = var0, align = 'Predict', 
#'                         npc = 20, add.Id = c("Set1", "Set2"), ncore = 2)
#' obj0 = scUMAP(obj0, npc = 8, use = "PLS", dist = 0.001, neighbors = 15)
#' DimPlot(obj0, slot = "cell.umap", colFactor = "Set", size = 2)
#' DimPlot(obj0, slot = "cell.umap", colFactor = "Group", size = 2, label = TRUE)

scMultiIntegrate <- function(
  objects, 
  eigens = 10, 
  add.Id = NULL, 
  var.gene = NULL, 
  align = 'OLS', 
  npc = 50, 
  adjust = TRUE, 
  ncore = 1, 
  seed = 123
  ) {
  
  pboptions(type = "txt", style = 3)
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
    i = j = l = p = 1L
    logcount0 = Var0 = gene0 = genes = type0 = coldata0 = list()
    Group = design0 = coef0 = genel = ranki = vector()
    genem = matrix()
    
    while(i <= nset){
      
      object = objects[[i]]
      Var0[[i]] = object@vargene
      logcount0[[i]] = object@assay$logcount
      coldata0[[i]] = object@coldata
      gene0[[i]] = rownames(object@assay$logcount)
      type0[[i]] = colnames(object@coldata)
      i = i + 1L
      
    }
    
    names(logcount0) = names(Var0) = names(coldata0) = names(objects)
    type0 = Reduce(intersect, type0)
    rm(objects, object)
    genes = Reduce(intersect, gene0)
    gene0 = unique(unlist(gene0))
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
        var.gene = genes
      }
      
    } else {
      var.gene = var.gene
    }
    var.gene = intersect(var.gene, gene0)
    rm(Var0)
    
    # rank datasets
    while(l <= length(add.Id)){
      logcountl = logcount0[[l]]
      coldatal = coldata0[[l]][,type0]
      coldatal$Set = add.Id[l]
      colnames(logcountl) = rownames(coldatal) = coldatal[,'Barcode'] = paste0(add.Id[l], "_", colnames(logcountl))
      genel = gene0[!gene0 %in% rownames(logcountl)]
      genem = Matrix::spMatrix(length(genel), dim(logcountl)[2])
      rownames(genem) = genel
      logcountl = rbind(logcountl, genem)
      logcount0[[l]] = logcountl[gene0,]
      coldata0[[l]] = coldatal
      rm(logcountl, coldatal, genel, genem)
      l = l + 1L
    }
    
    
    ### core integration
    logcount1L = logcount0[[1L]]
    logcount1L = logcount1L[var.gene,]
    seed = as.numeric(seed)
    set.seed(seed)
    
    logcount1L = scale(logcount1L, center = TRUE, scale = TRUE)
    Var.pca1L = irlba(logcount1L, nv = eigen0)
    pb = txtProgressBar(
      min = 0, max = length(add.Id), 
      initial = length(add.Id), style = 3, char = "*",
      width = getOption("width") / 2
    )
    scale.beta = foreach(j = 1L:length(add.Id)) %dopar% {
      logcountl = logcount0[[j]][var.gene,]
      celll = colnames(logcountl)
      PCA1Lu = multiply_d_d(Var.pca1L$u, diag(Var.pca1L$d, nrow = eigen0))
      logcountl = scale(logcountl, center = TRUE, scale = TRUE)
      # logcountl = cent_sp_d(logcountl)
      PCAjv = crossprod_d_d(PCA1Lu, logcountl) / (Var.pca1L$d)^2
      rm(logcountl)
      beta0 = multiply_d_d(Var.pca1L$v, PCAjv)
      beta0 = as.matrix(beta0)
      rownames(beta0) = colnames(logcount1L)
      colnames(beta0) = celll
      rm(PCA1Lu, PCAjv)
      return(beta0)
    }
    close(pb)
    rm(Var.pca1L, logcount1L)
    
    
    ### logcount alignment
    if(align == 'Predict'){
      
      # RPCI prediction
      pb = txtProgressBar(
        min = 0, max = length(add.Id), 
        initial = length(add.Id), style = 3, char = "+",
        width = getOption("width") / 2
      )
      logcount0 = foreach(j = 1L:length(add.Id)) %dopar% {
        cell0 = colnames(logcount0[[j]])
        i = logcount0[[j]]@i
        p = logcount0[[j]]@p
        keep = Matrix::which(logcount0[[j]] > 0)
        x = multiply_sp_d_v(logcount0[[1L]], scale.beta[[j]])[keep]
        lmax = quantile(x, probs = (length(x) - 10) / length(x))
        x[x > lmax] = lmax
        x[x < 0] = 0
        logcount0[[j]] = sparseMatrix(i = i+1, p = p, x = x, dims = c(length(gene0), length(cell0)), dimnames = list(gene0, cell0))
        rm(keep, x)
        logcount0[[j]] = Matrix::drop0(logcount0[[j]])
        return(logcount0[[j]])
      }
      close(pb)
      
      names(logcount0) = add.Id
      
      scale.beta = do.call(cbind, scale.beta)
      cell_name0 = colnames(scale.beta)
      scale.beta = irlba(scale.beta, nv = npc)
      scale.beta = as.matrix(scale.beta$v)
      colnames(scale.beta) = paste0('PC', 1L:npc)
      rownames(scale.beta) = cell_name0
      
    } else if(align == 'Optimal'){
      
      # Scale by experience
      scale.beta = do.call(cbind, scale.beta)
      cell_name0 = colnames(scale.beta)
      scale.beta = irlba(scale.beta, nv = npc)
      scale.beta = as.matrix(scale.beta$v)
      colnames(scale.beta) = paste0('PC', 1L:npc)
      rownames(scale.beta) = cell_name0
      
      ratio1L = sparseMatrixStats::rowQuantiles(logcount0[[1L]], probs = 0.975)
      ratio1L = as.vector(ratio1L)
      rowmean1L = Matrix::rowMeans(logcount0[[1L]])
      
      logcount0[2L:length(add.Id)] = pblapply(2L:length(add.Id), FUN = function(x){range_id(ratio1L, logcount0[[x]])}, cl = ncore)
      names(logcount0) = add.Id
      rm(ratio1L, rowmean1L)
      
    } else if(align == 'OLS') {
      
      # OLS
      scale.beta = do.call(cbind, scale.beta)
      cell_name0 = colnames(scale.beta)
      scale.beta = irlba(scale.beta, nv = npc)
      scale.beta = as.matrix(scale.beta$v)
      colnames(scale.beta) = paste0('PC', 1L:npc)
      rownames(scale.beta) = cell_name0
      
      coldata.OLS = do.call(rbind, coldata0)
      coldata.OLS$rank = 1:nrow(coldata.OLS)
      Group = factor(coldata.OLS$Set, levels = add.Id)
      contrasts(Group) = contr.sum(levels(Group))
      design0 = model.matrix(~Group)
      # Group = design0[, -1, drop = FALSE]
      Group = list()
      for(i in add.Id){
        ranki = coldata.OLS$rank[coldata.OLS$Set == i]
        Group[[i]] = design0[ranki, -1, drop = FALSE]
      }
      
      coef0 = pblapply(1L:length(gene0), FUN = function(y){lm_coef(design0, unlist(lapply(logcount0, function(x){x[y,]})))}, cl = ncore)
      coef0 = do.call(cbind, coef0)
      coef0 = t(coef0)
      coef0 = coef0[, -1L, drop = FALSE]
      coef0[is.na(coef0)] = 0
      rm(coldata.OLS, design0)
      
      logcount0 = pblapply(1L:length(add.Id), FUN = function(x){logcount_lm_id(logcount0[[x]], Group[[x]], coef0)}, cl = ncore)
      names(logcount0) = add.Id
      rm(Group, coef0)
      
    } else {
      stop('Please input align, Optimal or Predict')
    }
    
    
    ### combine data
    coldata.integrate = do.call(rbind, coldata0)
    rownames(coldata.integrate) = coldata.integrate[,'Barcode']
    coldata.integrate$Set = as.factor(coldata.integrate$Set)
    rm(coldata0)
    rowdata.integrate = data.frame(Symbol = gene0, RNA = "Gene Expression", Shared_Exp = "False", stringsAsFactors = F, row.names = gene0)
    rowdata.integrate$Shared_Exp[rowdata.integrate$Symbol %in% genes] = "True"
    
    # Output
    object = SingleCellData(assay = list(logcount = logcount0), rowdata = rowdata.integrate, coldata = coldata.integrate)
    object@metadata[['normalise']] = 'Yes'
    object@metadata[['Integration']] = data.frame(eigens = eigen0, npc = npc, method = "RPCI", align = align, stringsAsFactors = FALSE)
    object@vargene = var.gene
    object@DimReduction[['cell.pls']] = as.matrix(scale.beta)
    return(object)
    
  }
  
}



####################################################################################
#' Integrating Multiple Large Datasets
####################################################################################
#' 
#' The "scPLS" function can be used for data integration of multiple 
#' datasets, it is basically based on our new algorithm: reference principal 
#' components integration (RPCI). RPCI decomposes all the target datasets based 
#' on the reference. The output of this function can be used for low dimension 
#' visualization.
#' 
#' @rdname PLS-Integrating
#' @param objects The list of multiple RISC objects: 
#' list{object1, object2, object3, ...}. The first set is the reference to generate 
#' gene-eigenvectors.
#' @param eigens The number of eigenvectors used for data integration.
#' @param add.Id Add a vector of Id to label different datasets, a character vector.
#' @param var.gene Define the variable genes manually. Here input a vector of gene 
#' names as variable genes
#' @param npc The number of the PCs returns from "scMultiIntegrate" function, 
#' they are usually used for the subsequent analyses, like cell embedding and 
#' cell clustering.
#' @param adjust Whether adjust the number of eigenvectors.
#' @param ncore The number of multiple cores for data integration.
#' @param seed The random seed to keep consistent result.
#' @references Liu et al., Nature Biotech. (2021)
#' @name scPLS
#' @export
#' @examples
#' obj1 = raw.mat[[3]]
#' obj2 = raw.mat[[4]]
#' obj0 = list(obj1, obj2)
#' var0 = intersect(obj1@vargene, obj2@vargene)
#' PLS0 = scPLS(obj0, var.gene = var0, npc = 20, add.Id = c("Set1", "Set2"), ncore = 1)

scPLS <- function(
    objects, 
    eigens = 10, 
    add.Id = NULL, 
    var.gene = NULL, 
    npc = 100, 
    adjust = TRUE, 
    ncore = 1, 
    seed = 123
) {
  
  pboptions(type = "txt", style = 3)
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
    logcount0 = Var0 = gene0 = genes = list()
    genel = vector()
    genem = matrix()
    
    while(i <= nset){
      
      object = objects[[i]]
      Var0[[i]] = object@vargene
      logcount0[[i]] = object@assay$logcount
      gene0[[i]] = rownames(object@assay$logcount)
      i = i + 1L
      
    }
    
    names(logcount0) = names(Var0) = names(objects)
    rm(objects, object)
    genes = Reduce(intersect, gene0)
    gene0 = unique(unlist(gene0))
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
        var.gene = genes
      }
      
    } else {
      var.gene = var.gene
    }
    var.gene = intersect(var.gene, genes)
    rm(Var0)
    
    # rank datasets
    while(l <= length(add.Id)){
      logcountl = logcount0[[l]]
      colnames(logcountl) = paste0(add.Id[l], "_", colnames(logcountl))
      genel = gene0[!gene0 %in% rownames(logcountl)]
      genem = Matrix::spMatrix(length(genel), dim(logcountl)[2])
      rownames(genem) = genel
      logcountl = rbind(logcountl, genem)
      logcount0[[l]] = logcountl[gene0,]
      rm(logcountl, genel, genem)
      l = l + 1L
    }
    
    
    ### core integration
    logcount1L = logcount0[[1L]]
    logcount1L = logcount1L[var.gene,]
    seed = as.numeric(seed)
    set.seed(seed)
    
    logcount1L = scale(logcount1L, center = TRUE, scale = TRUE)
    Var.pca1L = irlba(logcount1L, nv = eigen0)
    pb = txtProgressBar(
      min = 0, max = length(add.Id), 
      initial = length(add.Id), style = 3, char = "*",
      width = getOption("width") / 2
    )
    scale.beta = foreach(j = 1L:length(add.Id)) %dopar% {
      logcountl = logcount0[[j]][var.gene,]
      celll = colnames(logcountl)
      PCA1Lu = multiply_d_d(Var.pca1L$u, diag(Var.pca1L$d, nrow = eigen0))
      logcountl = scale(logcountl, center = TRUE, scale = TRUE)
      # logcountl = cent_sp_d(logcountl)
      PCAjv = crossprod_d_d(PCA1Lu, logcountl) / (Var.pca1L$d)^2
      rm(logcountl)
      beta0 = multiply_d_d(Var.pca1L$v, PCAjv)
      beta0 = as.matrix(beta0)
      rownames(beta0) = colnames(logcount1L)
      colnames(beta0) = celll
      rm(PCA1Lu, PCAjv)
      return(beta0)
    }
    close(pb)
    rm(Var.pca1L, logcount1L)
    
    scale.beta = do.call(cbind, scale.beta)
    cell_name0 = colnames(scale.beta)
    scale.beta = irlba(scale.beta, nv = npc)
    scale.beta = as.matrix(scale.beta$v)
    colnames(scale.beta) = paste0('PC', 1L:npc)
    rownames(scale.beta) = cell_name0
    return(as.matrix(scale.beta))
    
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
#' @param res The resolution of cluster searched for, works in "louvain" method.
#' @param method The method of cell clustering for individual datasets.
#' @param algorithm The algorithm for knn, the default is "kd_tree", all options: 
#' "kd_tree", "cover_tree", "CR", "brute".
#' @param ncore The number of multiple cores for testing.
#' @param minPC The minimal PCs for detecting cell clustering.
#' @param Std.cut The cutoff of standard deviation of the PCs.
#' @param bin The bin number for calculating cell clustering.
#' @importFrom stats ks.test reshape
#' @references Liu et al., Nature Biotech. (2021)
#' @name InPlot
#' @export

InPlot <- function(
  object = NULL, 
  var.gene = NULL, 
  Colors = NULL, 
  nPC = 20, 
  neighbor = 30, 
  res = 1.0,
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
    res = as.numeric(res)
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
    vari = scale(vari, center = TRUE, scale = TRUE)
    PCi = irlba(vari, nv = npc)
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
      E(ki)$weight = 1/(1 + as.vector(k0$nn.dist))
      ki = simplify(ki)
      ki = cluster_louvain(ki, resolution = res)
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
  cmean = sparseMatrixStats::rowMedians(matrix(unlist(clust0), nset, ))
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
#' @references De-Jong et al. (1993)
#' @name SIMPLS

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
  S = crossprod_d_d(X, Y)
  
  i = 1L
  while(i <= npcs){
    
    if(ycol < xcol){
      set.seed(seed)
      q.a = eigen(crossprod(S), symmetric = TRUE)$vectors[,1]
    } else {
      q.a = c(crossprod_d_d(S, eigen(tcrossprod(S), symmetric = TRUE)$vectors[,1]))
      q.a = q.a / sqrt(c(crossprod(q.a)))
    }
    
    r.a = multiply_d_d(S, q.a)
    t.a = multiply_d_d(X, r.a)
    t.a = t.a - mean(t.a)
    tnorm = sqrt(c(crossprod(t.a)))
    t.a = t.a / tnorm
    r.a = r.a / tnorm
    p.a = crossprod_d_d(X, t.a)
    q.a = crossprod_d_d(Y, t.a)
    v.a = p.a
    v.a = v.a - multiply_d_d(V, crossprod_d_d(V, p.a))
    v.a = v.a /sqrt(c(crossprod(v.a)))
    S = S - multiply_d_d(v.a, crossprod_d_d(v.a, S))
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

MSC <- function(X, Y){
  
  x = t(as.matrix(X))
  y = t(as.matrix(Y))
  ref = colMeans(x)
  z = cbind(1, ref)
  beta = t(solve(crossprod(z), t(multiply_d_d(y, z))))
  res = (y - beta[,1]) / beta[,2]
  res = t(as.matrix(res))
  return(res)
  
}



####################################################################################
####################################################################################
logcount_coef_id <- function(gene0, logcount1L, logcountl, beta0){
  cell0 = colnames(logcountl)
  keep = (logcountl == 0)
  logcountl = multiply_sp_d_sp(logcount1L, beta0)
  logcountl[logcountl <= 0] = 0
  logcountl[keep] = 0
  rownames(logcountl) = gene0
  colnames(logcountl) = cell0
  rm(keep)
  logcountl = as(logcountl, 'CsparseMatrix')
  return(logcountl)
}

range_id <- function(ratio1L, logcountl){
  ratiol = sparseMatrixStats::rowQuantiles(logcountl, probs = 0.975)
  ratiol = as.vector(ratiol)
  ratiom = ratio1L / ratiol
  ratiom[is.na(ratiom)] = ratiom[is.infinite(ratiom)] = ratiom[ratiom == 0] = 1
  rowmean0 = Matrix::rowMeans(logcountl)
  keep = (rowmean0 > 0)
  ratio.correct = glm(log2(ratiom[keep] + 1) ~ log2(rowmean0[keep] + 1), family = 'quasipoisson')$fitted.values
  max.correct = 2^max(ratio.correct) - 1
  min.correct = 2^min(ratio.correct) - 1
  ratiom[ratiom > max.correct] = max.correct
  ratiom[ratiom < min.correct] = min.correct
  logcountl = as(logcountl * ratiom, 'CsparseMatrix')
  rm(ratio1L, ratiol, ratiom, keep)
  return(logcountl)
}

lm_id <- function(id0, count0, design0){
  id0 = (id0[1]+1):id0[2]
  coef0 = sapply(id0, FUN = function(y){lm_coef(design0, unlist(lapply(count0, function(x){x[y,]})))})
  coef0 = t(as.matrix(coef0))
  return(coef0)
}

logcount_lm_id <- function(logcountl, Group, coef0){
  keep = Matrix::which(logcountl > 0)
  batch0 = multiply_d_d(coef0, t(Group))
  logcountl@x = logcountl@x - as.vector(batch0)[keep]
  logcountl@x[logcountl@x < 0] = 0
  rm(keep, batch0, Group, coef0)
  logcountl = Matrix::drop0(logcountl)
  return(logcountl)
}

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
  M = tcrossprod(Z)
  
  dis = 10
  ctl = matrix(10, zrow, 1)
  ctl0 = ctl
  i = 1L
  while(dis > 1e-4 & i <= 300L){
    set.seed(123)
    mcsvd = svd(multiply_d_d(M, ctl))
    sam = tcrossprod_d_d(mcsvd$u, mcsvd$v)
    Ms = multiply_d_d(M, sam)
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
  
  xsd = sparseMatrixStats::colSds(x)
  if(scale){
    x = t(t(x)/sparseMatrixStats::colSds(x))
  } else {
    x = x
  }
  ysd = rep(1, ycol)
  
  beta = matrix(0, xcol, ycol)
  x0  = x
  y0 = y
  i = 1L
  while(i <= K){
    Z = crossprod_d_d(x0, y0)
    tag = SPLSDV(Z, Lamda = Lamda)
    sel = unique(IP[tag != 0 | beta[,1] != 0])
    xsel = x[, sel, drop = FALSE]
    coefxy = SIMPLS(xsel, y, npcs = min(i, length(sel)))
    beta[sel,] = matrix(coefxy, length(sel), ycol)
    y0 = y - multiply_d_d(x, beta)
    i = i + 1L
  }
  return(beta)
}


