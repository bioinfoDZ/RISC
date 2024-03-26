####################################################################################
#' Utilities Subset data
####################################################################################
#' 
#' The "Subset" function can abstract a data subset from the full dataset, this 
#' function not only collect the subset of coldata and rowdata, but also abstract 
#' raw counts/UMIs. Meanwhile, after "Subset" function, RISC object need to be 
#' normalized and scaled one more time.
#' 
#' @rdname Subset
#' @param object RISC object: a framework dataset.
#' @param cells The cells are directly used for collecting a data subset.
#' @param genes The genes are directly used for collecting a data subset.
#' @name SubSet
#' @export
#' @examples 
#' # RISC object
#' obj0 = raw.mat[[5]]
#' obj0
#' cell1 = rownames(obj0@coldata)[1:15]
#' obj1 = SubSet(obj0, cells = cell1)
#' obj1

SubSet <- function(object, cells = NULL, genes = NULL){
  
  coldata0 = object@coldata
  rowdata0 = object@rowdata
  raw.assay = object@assay
  DimReduction0 = object@DimReduction
  
  if(is.null(object)){
    stop('Please input a RISC object')
  } else if(!is.null(cells) & is.null(genes)){
    
    coldata0 = coldata0[rownames(coldata0) %in% cells,]
    
    if('Integration' %in% names(object@metadata)){
      
      name0 = 'logcount'
      raw.assay = raw.assay$logcount
      raw.assay = lapply(raw.assay, FUN = function(y){y[, colnames(y) %in% rownames(coldata0), drop = F]})
      raw.assay[sapply(raw.assay, function(x){dim(x)[2] == 0})] = NULL
      rowsum0 = lapply(raw.assay, FUN = function(y){Matrix::rowSums(y > 0)})
      rowsum0 = do.call(cbind, rowsum0)
      keep = Matrix::rowSums(rowsum0) > 0
      gene0 = rownames(object@rowdata)[keep]
      raw.assay = lapply(raw.assay, FUN = function(y){y[gene0, , drop = FALSE]})
      raw.assay[sapply(raw.assay, function(x){dim(x)[1] == 0})] = NULL
      
    } else {
      
      name0 = names(raw.assay)
      raw.assay = raw.assay[[name0]]
      raw.assay = raw.assay[, rownames(coldata0), drop = FALSE]
      keep = Matrix::rowSums(raw.assay > 0) > 0
      gene0 = rownames(object@rowdata)[keep]
      raw.assay = raw.assay[gene0,]
      
      if(name0 == 'count'){
        coldata0$UMI = Matrix::colSums(raw.assay)
        coldata0$nGene = Matrix::colSums(raw.assay > 0)
      } else {
        coldata0 = coldata0
      }
      
    }
    
    rowdata0 = rowdata0[gene0,]
    
    if(length(DimReduction0) > 0){
      
      for(key0 in names(DimReduction0)){
        if(!key0 %in% c("var.pca", "gene.pca")){
          DimReduction0[[key0]] = DimReduction0[[key0]][rownames(DimReduction0[[key0]]) %in% cells,]
        }
        else{
          DimReduction0[[key0]] = DimReduction0[[key0]]
        }
      }
      
    } else {
      DimReduction0 = DimReduction0
    }
    
  } else if(!is.null(cells) & !is.null(genes)){
    
    coldata0 = coldata0[rownames(coldata0) %in% cells,]
    rowdata0 = rowdata0[rownames(rowdata0) %in% genes,]
    gene0 = rownames(rowdata0)
    
    if('Integration' %in% names(object@metadata)){
      
      name0 = 'logcount'
      raw.assay = raw.assay$logcount
      raw.assay = lapply(raw.assay, FUN = function(y){y[gene0, colnames(y) %in% rownames(coldata0), drop = FALSE]})
      raw.assay[sapply(raw.assay, function(x){dim(x)[1] == 0})] = NULL
      raw.assay[sapply(raw.assay, function(x){dim(x)[2] == 0})] = NULL
      
    } else {
      
      name0 = names(raw.assay)
      raw.assay = raw.assay[[name0]]
      raw.assay = raw.assay[gene0, rownames(coldata0), drop = FALSE]
      
      if(name0 == 'count'){
        coldata0$UMI = Matrix::colSums(raw.assay)
        coldata0$nGene = Matrix::colSums(raw.assay > 0)
      } else {
        coldata0 = coldata0
      }
      
    }
    
    if(length(DimReduction0) > 0){
      
      for(key0 in names(DimReduction0)){
        if(!key0 %in% c("var.pca", "gene.pca")){
          DimReduction0[[key0]] = DimReduction0[[key0]][rownames(DimReduction0[[key0]]) %in% cells,]
        }
        else{
          DimReduction0[[key0]] = DimReduction0[[key0]]
        }
      }
      
    } else {
      DimReduction0 = DimReduction0
    }
    
  } else if(is.null(cells) & !is.null(genes)){
  	
    rowdata0 = rowdata0[rownames(rowdata0) %in% genes,]
    gene0 = rownames(rowdata0)
    
    if('Integration' %in% names(object@metadata)){
      
      name0 = 'logcount'
      raw.assay = raw.assay$logcount
      raw.assay = lapply(raw.assay, FUN = function(y){y[gene0, , drop = FALSE]})
      raw.assay[sapply(raw.assay, function(x){dim(x)[1] == 0})] = NULL
      
    } else {
      
      name0 = names(raw.assay)
      raw.assay = raw.assay[[name0]]
      raw.assay = raw.assay[gene0, , drop = FALSE]
      
      if(name0 == 'count'){
        coldata0$UMI = Matrix::colSums(raw.assay)
        coldata0$nGene = Matrix::colSums(raw.assay > 0)
      } else {
        coldata0 = coldata0
      }
      
    }
    
  } else {stop('No parameters')}
  
  object@coldata = data.frame(coldata0)
  object@rowdata = data.frame(rowdata0)
  object@DimReduction = DimReduction0
  object@assay = list(raw.assay)
  names(object@assay) = name0
  object@cluster = factor()
  
  return(object)
  
}



####################################################################################
#' Utilities Add Factors
####################################################################################
#' 
#' The "AddFactor" function can add factors to the full dataset, this function 
#' can add one or more factors into coldata. Here the row.names/names of factor 
#' matrix/vector should be equal to the row.names of coldata of RISC object.
#' 
#' @rdname AddFactor
#' @param object RISC object: a framework dataset.
#' @param colData Input the names that will be added into coldata of RISC object, 
#' it should be characters, as the col.names of coldata.
#' @param rowData Input the names that will be added into rowdata of RISC object,
#' it should be characters, as the col.names of rowdata.
#' @param value The factor vector or data.frame that will be added into coldata 
#' or rowdata, the vector/data.frame should have equal names/row.names to the 
#' row.names coldata or rowdata of RISC object. The input: vector or data.frame.
#' @name AddFactor

AddFactor <- function(object, colData = NULL, rowData = NULL, value = NULL){
  
  coldata0 = as.data.frame(object@coldata)
  col.name0 = colnames(coldata0)
  rowdata0 = as.data.frame(object@rowdata)
  row.name0 = colnames(rowdata0)
  
  if(is.null(object)) {
    stop('Please input a RISC object')
  } else if(!is.null(colData)) {
    
    colData = as.character(colData)
    if(inherits(value) == "data.frame") {
      coldata1 = data.frame(coldata0, colData = value)
      colnames(coldata1) = c(col.name0, colData)
      object@coldata = coldata1
    } else {
      coldata1 = data.frame(coldata0, colData = value)
      colnames(coldata1) = c(col.name0, colData)
      object@coldata = coldata1
    }
    
  } else if(!is.null(rowData)){
    
    rowData = as.character(rowData)
    if(inherits(value) == "data.frame") {
      rowdata1 = data.frame(rowdata0, colData = value)
      colnames(rowdata1) = c(row.name0, rowData)
      object@rowdata = rowdata1
    } else {
      rowdata1 = data.frame(rowdata0, colData = value)
      colnames(rowdata1) = c(row.name0, rowData)
      object@rowdata = rowdata1
    }
    
  } else {
    stop('Value should be vector or data.frame, of which the name/row.names should be equal to row.names of coldata or rowdata')
  }
  
  return(object)
  
}


