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
#' obj0 = raw.mat[[3]]
#' obj0
#' cell1 = rownames(obj0@coldata)[1:25]
#' obj1 = SubSet(obj0, cells = cell1)
#' obj1

SubSet <- function(object, cells = NULL, genes = NULL){
  
  coldata0 = object@coldata
  rowdata0 = object@rowdata
  raw.count = object@assay$count
  DimReduction0 = object@DimReduction
  
  if(is.null(object)){
    stop('Please input a RISC object')
  } else if(!is.null(cells) & is.null(genes)){
    
    coldata0 = coldata0[rownames(coldata0) %in% cells,]
    raw.count = raw.count[,rownames(coldata0)]
    # min.cell = floor(nrow(coldata0)/500)
    keep = Matrix::rowSums(raw.count > 0) > 0
    raw.count = raw.count[keep,]
    rowdata0 = rowdata0[keep,]
    coldata0$scUMI = Matrix::colSums(raw.count)
    coldata0$ngene = Matrix::colSums(raw.count > 0)
    
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
    raw.count = raw.count[rownames(raw.count) %in% genes, rownames(coldata0)]
    coldata0$scUMI = Matrix::colSums(raw.count)
    coldata0$ngene = Matrix::colSums(raw.count > 0)
    
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
    raw.count = raw.count[rownames(raw.count) %in% genes,]
    coldata0$scUMI = Matrix::colSums(raw.count)
    coldata0$ngene = Matrix::colSums(raw.count > 0)

  } else {stop('No parameters')}
  
  object@coldata = data.frame(coldata0)
  object@rowdata = data.frame(rowdata0)
  object@DimReduction = DimReduction0
  assay = names(object@assay)
  
  for(i in assay){
    count = object@assay[[i]]
    count0 = count[rownames(count) %in% rownames(raw.count), colnames(count) %in% colnames(raw.count)]
    object@assay[[i]] = count0
  }
  
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
#' @export
#' @examples 
#' # RISC object
#' obj0 = raw.mat[[3]]
#' obj0
#' add0 = rep(0, 45)
#' obj0 = AddFactor(obj0, colData = 'add_factor', value = add0)
#' obj0

AddFactor <- function(object, colData = NULL, rowData = NULL, value = NULL){
  
  coldata0 = as.data.frame(object@coldata)
  col.name0 = colnames(coldata0)
  rowdata0 = as.data.frame(object@rowdata)
  row.name0 = colnames(rowdata0)
  
  if(is.null(object)) {
    stop('Please input a RISC object')
  } else if(!is.null(colData)) {
    
    colData = as.character(colData)
    if(class(value) == "data.frame") {
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
    if(class(value) == "data.frame") {
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


