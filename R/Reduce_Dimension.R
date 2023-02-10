####################################################################################
#' Dimension Reduction.
####################################################################################
#' 
#' Based on highly variably expressed genes of the datasets, RISC calculates the 
#' principal components (PCs) of the cells using prcomp functions. The major PCs, 
#' which explain most gene expression variance, are used for dimension reduciton. 
#' 
#' @rdname PCA
#' @param object RISC object: a framework dataset.
#' @param npc The number of PCs will be generated based on highly variable genes 
#' (usually < 1,500), npc equal to the first 20 PCs as the default.
#' @return RISC single cell dataset, the DimReduction slot.
#' @references Jolliffe et al. (2016)
#' @references Alter et al., PNAS (2000)
#' @references Gonzalez et al., JSS (2008)
#' @references Mevik et al., JSS (2007)
#' @name scPCA
#' @export
#' @examples 
#' # RISC object
#' obj0 = raw.mat[[3]]
#' obj0 = scPCA(obj0, npc = 10)

scPCA <- function(object, npc = 20){
  
  if(!length(object@vargene) > 0){
    stop('Please disperse object first')
  } else {
    
    count = object@assay$logcount
    var = count[object@vargene,]
    varpc = irlba(var, nv = npc, center = TRUE)
    cell.pca = as.matrix(varpc$v)
    gene.pca = as.matrix(varpc$u)
    var.pca = varpc$d^2 / sum(varpc$d^2)
    rownames(cell.pca) = colnames(var)
    rownames(gene.pca) = rownames(var)
    colnames(cell.pca) = colnames(gene.pca) = names(var.pca) = paste0('PC', 1L:npc)
    object@DimReduction[['cell.pca']] = cell.pca
    object@DimReduction[['var.pca']] = var.pca
    object@DimReduction[['gene.pca']] = gene.pca
    return(object)
    
  }
}



####################################################################################
#' Dimension Reduction.
####################################################################################
#' 
#' The UMAP is calculated based on the eigenvectors of single cell dataset, and the 
#' user can select the eigenvectors manually. Of note, the selected eigenvectors 
#' directly affect UMAP values. 
#' For the integrated data (the result of "scMultiIntegrate" funciton), RISC utilizes
#' the PCR output "PLS" to calculate the UMAP, therefore, the user has to input "PLS"
#' in "use = ", instead of the default parameter "PCA".
#' 
#' @rdname UMAP
#' @param object RISC object: a framework dataset.
#' @param npc The number of the PCs (or the PLS) using for UMAP, the default is 20, 
#' but need to be modified by the users. The PCA for individual dataset, while PLS 
#' for the integrated data.
#' @param embedding The number of components UMAP output.
#' @param use What components used for UMAP: PCA or PLS.
#' @param neighbors The n_neighbors parameter of UMAP.
#' @param dist The min_dist parameter of UMAP.
#' @param seed The random seed to keep tSNE result consistent.
#' @return RISC single cell dataset, the DimReduction slot.
#' @importFrom umap umap
#' @references Becht et al., Nature Biotech. (2018)
#' @name scUMAP
#' @export
#' @examples 
#' # RISC object
#' obj0 = raw.mat[[3]]
#' obj0 = scPCA(obj0, npc = 10)
#' obj0 = scUMAP(obj0, npc = 3)
#' DimPlot(obj0, slot = "cell.umap", colFactor = 'Group', size = 2)

scUMAP <- function(
  object, npc = 20, 
  embedding = 2, 
  use = 'PCA', 
  neighbors = 15, 
  dist = 0.1, 
  seed = 123
  ) {
  
  if(use == 'PCA'){
    
    if(length(object@DimReduction$cell.pca) == 0){
      stop('Please disperse object first')
    } else {
      pca0 = FALSE
      pca_center0 = FALSE
      pca_scale0 = FALSE
      cell.pc = object@DimReduction$cell.pca[,1:npc]
    }
    
  } else if(use == 'PLS'){
    
    if(length(object@DimReduction$cell.pls) == 0){
      stop('Please integrate objects first')
    } else {
      pca0 = TRUE
      pca_center0 = TRUE
      pca_scale0 = TRUE
      cell.pc0 = object@DimReduction$cell.pls
      
      if(npc <= ncol(cell.pc0)){
        cell.pc = cell.pc0[,1:npc]
      } else {
        scale.beta = object@metadata$Beta
        scale.beta0 = irlba(scale.beta, nv = npc)
        cell.pc = rbind(scale.beta0$u, scale.beta0$v)
        rownames(cell.pc) = c(rownames(scale.beta), colnames(scale.beta))
        colnames(cell.pc) = paste0('PC', 1L:npc)
        # cell.pc = scale.beta1[order(rownames(scale.beta1), decreasing = F),]
      }
      
    }
    
  } else {
    stop('Input use, PCA or PCR')
  }
  
  set.seed(as.numeric(seed))
  embedding = as.integer(embedding)
  neighbor0 = as.integer(neighbors)
  dist0 = as.numeric(dist)
  umap0 = umap(as.matrix(cell.pc), method = 'naive', n_components = embedding, min_dist = dist0, n_neighbors = neighbor0)
  cell.umap = as.matrix(umap0$layout)
  rownames(cell.umap) = rownames(cell.pc)
  colnames(cell.umap) = paste0('UMAP', 1L:embedding)
  object@DimReduction[['cell.umap']] = cell.umap
  return(object)
  
}



####################################################################################
#' Dimension Reduction.
####################################################################################
#' 
#' The t-SNE is calculated based on the eigenvectors of single cell dataset, and 
#' the user can select the eigenvectors manually. Of note, the selected eigenvectors 
#' directly affect t-SNE values. 
#' For the integrated data (the result of "scMultiIntegrate" funciton), RISC utilizes
#' the PCR output "PLS" to calculate the t-SNE, therefore, the user has to input 
#' "PLS" in "use = ", instead of the defaut parameter "PCA".
#' 
#' @rdname tSNE
#' @param object RISC object: a framework dataset.
#' @param npc The number of PCs (or PLS) using for t-SNE, the default is 20, 
#' but need to be modified by the users. The PCA for individual dataset, while 
#' PLS for the integrated data.
#' @param embedding The number of components t-SNE output.
#' @param use What components used for t-SNE: PCA or PLS.
#' @param perplexity Perplexity parameter: if the cell numbers are small, 
#' decrease this parameter, otherwise tSNE cannot be calculated.
#' @param seed The random seed to keep tSNE result consistent.
#' @return RISC single cell dataset, the DimReduction slot.
#' @importFrom Rtsne Rtsne
#' @references Laurens van der Maaten, JMLR (2014)
#' @name scTSNE
#' @export
#' @examples 
#' # RISC object
#' obj0 = raw.mat[[3]]
#' obj0 = scPCA(obj0, npc = 10)
#' obj0 = scTSNE(obj0, npc = 4, perplexity = 10)
#' DimPlot(obj0, slot = "cell.tsne", colFactor = 'Group', size = 2)

scTSNE <- function(object, npc = 20, embedding = 2, use = 'PCA', perplexity = 30, seed = 123) {
  
  npc = as.integer(npc)
  perplexity = as.integer(perplexity)
  
  if(use == 'PCA'){
    
    if(length(object@DimReduction$cell.pca) == 0){
      stop('Please disperse object first')
    } else {
      pca0 = FALSE
      pca_center0 = FALSE
      pca_scale0 = FALSE
      cell.pc = object@DimReduction$cell.pca[,1:npc]
    }
    
  } else if(use == 'PLS'){
    
    if(length(object@DimReduction$cell.pls) == 0){
      stop('Please integrate objects first')
    } else {
      pca0 = TRUE
      pca_center0 = TRUE
      pca_scale0 = TRUE
      cell.pc0 = object@DimReduction$cell.pls
      
      if(npc <= ncol(cell.pc0)){
        cell.pc = cell.pc0[,1:npc]
      } else {
        scale.beta = object@metadata$Beta
        scale.beta0 = irlba(scale.beta, nv = npc)
        cell.pc = rbind(scale.beta0$u, scale.beta0$v)
        rownames(cell.pc) = c(rownames(scale.beta), colnames(scale.beta))
        colnames(cell.pc) = paste0('PC', 1L:npc)
        # cell.pc = scale.beta1[order(rownames(scale.beta1), decreasing = F),]
      }
      
    }
    
  } else {
    stop('Input use, PCA or PCR')
  }
  
  set.seed(as.numeric(seed))
  embedding = as.integer(embedding)
  tsne0 = Rtsne(as.matrix(cell.pc), dims = embedding, pca = pca0, pca_center = pca_center0, pca_scale = pca_scale0, perplexity = perplexity)
  cell.tsne = as.matrix(tsne0$Y)
  rownames(cell.tsne) = rownames(cell.pc)
  colnames(cell.tsne) = paste0('tSNE', 1L:embedding)
  object@DimReduction[['cell.tsne']] = cell.tsne
  return(object)
  
}


