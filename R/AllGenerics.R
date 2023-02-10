####################################################################################
#' RISC data
####################################################################################
#' 
#' The RISC object contains all the basic information used in single cell RNA-seq 
#' analysis, including raw counts/UMIs, normalized gene values, dimension reduction, 
#' cell clustering, and so on. The framework of RISC object is a S4 dataset, 
#' consisting of assay, coldata, rowdata, metadata, vargene, cluster, and 
#' DimReduction.
#'
#' @rdname setClass
#' @return RISC object: a S4 framework dataset
#' @slot assay The list of gene counts/UMIs: raw and normalized counts
#' @slot coldata The data.frame with cell information, such as cell types, stages, 
#' and other factors.
#' @slot rowdata The data.frame with gene information, such as coding or non-coding 
#' genes.
#' @slot metadata The data.frame with meta value.
#' @slot vargene The highly variable gene. 
#' These genes are utilized in dimension reduction.
#' @slot cluster The cell clustering information: include three algorithms for 
#' clustering.
#' @slot DimReduction The values of dimension reduction.
#' @docType class
#' @exportClass RISCdata

setClass(
  'RISCdata', slots = list(
    assay = 'list',
    coldata = 'data.frame',
    rowdata = 'data.frame',
    metadata = 'list',
    cluster = 'factor',
    DimReduction = 'list',
    vargene = 'vector'
  )
)

#' RISC data
#' 
#' This will show the full information of RISC object, including the number of cells, 
#' the number of genes, any biological or statistical information of cells or genes.
#' 
#' @rdname setMethod
#' @name RISCdata
#' @aliases RISC object
#' @param object RISC object: a S4 framework dataset
#' @docType methods

.RISC_show <- function(object){
  cat(
    "SingleCell-Dataset", '\n',
    "RISC v1.6", '\n',
    c('assay:', names(object@assay)), '\n',
    c(paste0('colData: ', '(', nrow(object@coldata), ')'), colnames(object@coldata)), '\n',
    c(paste0('rowData: ', '(', nrow(object@rowdata), ')'), colnames(object@rowdata)), '\n',
    'DimReduction', '\n',
    'Cell-Clustering', '\n'
  )
}

setMethod(
  f = 'show',
  signature = 'RISCdata',
  definition = .RISC_show
)




