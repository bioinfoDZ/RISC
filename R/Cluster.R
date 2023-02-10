####################################################################################
#' Clustering cells
####################################################################################
#' 
#' In RISC, two different methods are provided to cluster cells, all of them are 
#' widely used in single cells. The first method is "louvain" based on cell 
#' eigenvectors, and the other is "density" which calculates cell clusters using 
#' low dimensional space. 
#' 
#' @rdname Cluster
#' @param object RISC object: a framework dataset.
#' @param method The methods for clustering cells, density and louvain. 
#' The "density" is based on the slot "cell.umap" or other low dimensional 
#' space; while "louvain" based on "cell.pca" (individual data) or "cell.pls" 
#' (for integration data).
#' @param slot The dimension_reduction slot for cell clustering. The default is 
#' "cell.umap" under RISC object "DimReduction" item for UMAP method, but the 
#' customer can add new dimension_reduction method under DimReduction and use it.
#' @param neighbor The neighbor cells for "igraph" method.
#' @param algorithm The algorithm for knn, the default is "kd_tree", all options: 
#' "kd_tree", "cover_tree", "CR", "brute".
#' @param npc The number of PCA or PLS used for cell clustering.
#' @param k The number of cluster searched for, works in "density" method.
#' @param dc The distance used to generate random center points which affect 
#' clusters. If have no idea about this, do not input anything. Keep it as the 
#' default value for most users. Work for "density" method.
#' @param redo Whether re-cluster the cells.
#' @param random.seed The random seed, the default is 123.
#' @return RISC single cell dataset, the cluster slot. 
#' @name scCluster
#' @importFrom densityClust densityClust findClusters
#' @importFrom FNN get.knn
#' @importFrom igraph simplify graph_from_data_frame cluster_louvain
#' @references Blondel et al., JSTAT (2008)
#' @references Rodriguez et al., Sicence (2014)
#' @export
#' @examples 
#' # RISC object
#' obj0 = raw.mat[[3]]
#' obj0 = scPCA(obj0, npc = 10)
#' obj0 = scUMAP(obj0, npc = 3)
#' obj0 = scCluster(obj0, slot = "cell.umap", k = 3, method = 'density')
#' DimPlot(obj0, slot = "cell.umap", colFactor = 'Cluster', size = 2)

scCluster <- function(
  object, 
  slot = "cell.umap", 
  neighbor = 10, 
  algorithm = "kd_tree", 
  method = 'louvain', 
  npc = 20, 
  k = 10, 
  dc = NULL, 
  redo = TRUE, 
  random.seed = 123
  ) {
  
  set.seed(random.seed)
  k = as.integer(k)
  neighbor = as.integer(neighbor)
  algorithm = as.character(algorithm)
  npc = as.integer(npc)
  random.seed = as.integer(random.seed)
  
  if(is.null(dc)){
    dc = object@metadata$dcluster$dc
    if(isTRUE(redo)){
      dc = NULL
    } else {
      dc = dc
    }
  } else {
    dc = dc
  }
  
  if(!is.null(object@vargene) & !is.null(object@DimReduction)){
    
    slot0 = as.character(slot)
    dimReduce0 = object@DimReduction[[slot0]]
    if(is.null(dimReduce0)){
      stop("Do not include this dimention_reduction slot, try another one")
    } else if(ncol(dimReduce0) >= npc){
      count = as.matrix(dimReduce0[,1:npc])
    } else {
      count = as.matrix(dimReduce0)
    }
    
    if(method == 'density'){
      
      dist0 = dist(count)
      if(is.null(dc)){
        dataClust = densityClust(dist0, gaussian = TRUE)
      } else {
        dataClust = densityClust(dist0, dc = dc, gaussian = TRUE)
      }
      
      delta.rho = data.frame(rho = dataClust$rho, delta = dataClust$delta, stringsAsFactors = FALSE)
      delta.rho = delta.rho[order(delta.rho$delta, decreasing = TRUE),]
      delta.cut = delta.rho$delta[k + 1L]
      clust0 = findClusters(dataClust, 0, delta.cut)
      # object@metadata$dcluster = clust0
      object@cluster = object@coldata$Cluster = as.factor(clust0$clusters)
      names(object@cluster) = names(object@coldata$Cluster) = rownames(count)
      object@metadata[['clustering']] = data.frame(
        Method = 'densityClust', Distance = as.numeric(clust0$dc), 
        rho = as.numeric(clust0$threshold[1]), delta = as.numeric(clust0$threshold[2]), 
        stringsAsFactors = F
      )
      
    } else if(method == 'louvain') {
      
      clust0 = get.knn(count, k = neighbor, algorithm = algorithm)
      clust1 = data.frame(NodStar = rep(1L:nrow(count), neighbor), NodEnd = as.vector(clust0$nn.index), stringsAsFactors = FALSE)
      clust1 = graph_from_data_frame(clust1, directed = FALSE)
      clust1 = simplify(clust1)
      clust1 = cluster_louvain(clust1, weights = 1/(1 + as.vector(clust0$nn.dist)))
      object@cluster = object@coldata$Cluster = as.factor(clust1$membership)
      names(object@cluster) = names(object@coldata$Cluster) = rownames(count)
      object@metadata[['clustering']] = data.frame(Method = 'louvain', PCs = npc, Neighbors = neighbor, stringsAsFactors = F)
      
    } else {stop('A new method later')}
    
  } else {stop('Please calculate dispersion and dimention reduction first')}
  
  return(object)
  
}


