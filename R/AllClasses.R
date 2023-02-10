####################################################################################
#' Import single cell data
####################################################################################
#' 
#' The single cell RNA-seq (scRNA-seq) data can be imported in three different ways.
#' Primarily, we could import from 10X Genomics output directly by using 
#' "read10Xgenomics". The user only need to provide the folder path. Secondly, 
#' we could read data from HT-seq output by "readHTSeqdata", the user have to input 
#' the folder path. Lastly, we could input matrix, cell and genes mannually 
#' using "readscdata".
#' 
#' @useDynLib RISC
#' @importFrom Rcpp evalCpp
#' @importFrom methods as new
#' @importFrom utils head read.table
#' @importFrom Matrix readMM colSums rowSums
#' @rdname SingleCellData
#' @return SingleCellData
#' @param assay The list of gene counts.
#' @param coldata The data.frame with cell information.
#' @param rowdata The data.frame with gene information.
#' @name SingleCellData

SingleCellData <- function(assay, coldata, rowdata){
  object <- new(Class = 'RISCdata', assay = assay, coldata = coldata, rowdata = rowdata)
}



####################################################################################
#' Example data
####################################################################################
#' 
#' @docType data
#' @usage data(raw.mat)
#' @format A list including a simulated cell-gene matrix, columns for cells and 
#' rows for genes, a cell group and a batch information. 
"raw.mat"



####################################################################################
#' Import data from matrix, cell and genes directly.
####################################################################################
#' 
#' Import data set from matrix, cell and genes directly, the customer needs three 
#' files: a matrix file including gene expression values: raw counts/UMIs (rows for 
#' genes while columns for cells), a cell file (whose row.name are equal to the 
#' col.name of the matrix), and a gene file whose row.name are the same as the 
#' row.name of the matrix. If row.names of the gene matrix are Ensembl ID, the 
#' customer need to transfer them to gene symbols manually.
#' 
#' @rdname Import-Matrix
#' @param count Matrix with raw counts/UMIs.
#' @param cell Data.frame with cell Barcode, whose row.name are equal to the 
#' col.name of the matrix.
#' @param gene Data.frame with gene symbol, whose row.name are the same as the 
#' row.name of the matrix.
#' @param is.filter Remove not expressed genes.
#' @return RISC single cell dataset, including count, coldata, and rowdata.
#' @name readsc
#' @export
#' @examples 
#' mat0 = as.matrix(raw.mat[[1]])
#' coldata0 = as.data.frame(raw.mat[[2]])
#' coldata.obj = coldata0[coldata0$Batch0 == 'Batch3',]
#' matrix.obj = mat0[,rownames(coldata.obj)]
#' obj0 = readsc(count = matrix.obj, cell = coldata.obj, 
#'        gene = data.frame(Symbol = rownames(matrix.obj), 
#'        row.names = rownames(matrix.obj)), is.filter = FALSE)

readsc <- function(
  count, 
  cell, 
  gene, 
  is.filter = TRUE
) {
  
  if(exists("count") & exists("cell") & exists("gene")){
    
    if(all(colnames(count) == rownames(cell)) & all(rownames(count) == rownames(gene))){
      
      run.matrix = as(count, 'CsparseMatrix')
      # run.matrix = as.matrix(count)
      
      mito.gene = grep(pattern = '^mt-', x = rownames(run.matrix), ignore.case = TRUE, value = TRUE)
      run.cell = data.frame(scBarcode = rownames(cell), UMI = Matrix::colSums(run.matrix), nGene = Matrix::colSums(run.matrix > 0), row.names = rownames(cell), stringsAsFactors = FALSE)
      run.cell$mito = Matrix::colSums(run.matrix[rownames(run.matrix) %in% mito.gene,]) / run.cell$UMI
      run.cell = cbind.data.frame(run.cell, cell)
      
      run.gene = data.frame(Symbol = rownames(run.matrix), RNA = "Gene Expression", row.names = rownames(run.matrix), stringsAsFactors = FALSE)
      run.gene$nCell = Matrix::rowSums(run.matrix > 0)
      if(is.filter){
        run.gene = run.gene[run.gene$nCell > 0,]
      } else {
        run.gene = run.gene
      }
      run.matrix = run.matrix[rownames(run.matrix) %in% run.gene$Symbol,]
      
      SingleCellData(assay = list(count = as(run.matrix, 'CsparseMatrix')), rowdata = data.frame(run.gene, stringsAsFactors = FALSE), coldata = data.frame(run.cell, stringsAsFactors = FALSE))
      
    } else {
      stop('Matrix colnames or rownames are not equal to cell name or gene name')
    }
    
  } else {
    stop('No matrix, cell or gene is found here')
  }
  
}



####################################################################################
#' Import data from 10X Genomics output (tsv-mtx).
####################################################################################
#' 
#' Import data directly from 10X Genomics output, usually using filtered gene 
#' matrices which contains three files: matrix.mtx, barcode.tsv and gene.tsv. 
#' The user only need to input the directory into "data.path". If not the original 
#' 10X Genomics output, the user have to make sure the barcode.tsv and gene.tsv 
#' without col.names, the barcode.tsv at least contains one column for cell 
#' barcode, and the gene.tsv has two columns for gene Ensembl ID and Symbol.
#' 
#' @rdname Import-10X-mtx
#' @param data.path Directory containing the filtered 10X Genomics output, 
#' including three files: matrix.mtx, barcode.tsv (without colnames) and gene.tsv 
#' (without colnames).
#' @param sep The sep can be changed by the users
#' @param is.filter Remove not expressed genes.
#' @return RISC single cell dataset, including count, coldata, and rowdata.
#' @name read10X_mtx
#' @export

read10X_mtx <- function(
  data.path, 
  sep = '\t', 
  is.filter = TRUE
  ) {
  
  if(!exists("data.path")){
    stop('Please input data.path')
  } else {
    data.path = as.character(data.path)
  }
  
  sep0 = sep
  files = list.files(path = data.path, full.names = TRUE)
  file.matrix = grep('matrix.mtx', files, ignore.case = TRUE, value = TRUE)
  file.gene = grep(pattern = 'features.tsv', files, ignore.case = TRUE, value = TRUE)
  file.cell = grep('barcodes.tsv', files, ignore.case = TRUE, value = TRUE)
  
  if(length(file.matrix) == 1 & length(file.gene) == 1 & length(file.cell) == 1){
    
    run.matrix = readMM(file = file.matrix)
    run.matrix = as(run.matrix, 'CsparseMatrix')
    
    run.gene = read.table(file = file.gene, header = FALSE, sep = sep0, stringsAsFactors = FALSE)
    run.gene = data.frame(run.gene, stringsAsFactors = FALSE)
    if(ncol(run.gene) > 2){
      colnames(run.gene) = c('Ensembl', 'Symbol', 'RNA')
    } else {
      colnames(run.gene) = c('Ensembl', 'Symbol')
      run.gene$RNA = "Gene Expression"
    }
    run.gene$nCell = Matrix::rowSums(run.matrix > 0)
    run.gene$Symbol = make.unique(run.gene$Symbol)
    rownames(run.matrix) = rownames(run.gene) = run.gene$Symbol
    if(is.filter){
      run.gene = run.gene[run.gene$nCell > 0,]
    } else {
      run.gene = run.gene
    }
    run.matrix = run.matrix[rownames(run.matrix) %in% run.gene$Symbol,]
    
    mito.gene = grep(pattern = '^mt-', x = rownames(run.matrix), ignore.case = TRUE, value = TRUE)
    run.cell = read.table(file = file.cell, header = FALSE, sep = sep0, stringsAsFactors = FALSE)
    # run.cell = sapply(run.cell$V1, function(x){strsplit(x, '-', fixed = T)[[1]][[1]]})
    run.cell = data.frame(scBarcode = as.character(run.cell$V1), UMI = Matrix::colSums(run.matrix), nGene = Matrix::colSums(run.matrix > 0), stringsAsFactors = FALSE)
    run.cell$mito = Matrix::colSums(run.matrix[rownames(run.matrix) %in% mito.gene,]) / run.cell$UMI
    colnames(run.matrix) = rownames(run.cell) = run.cell$scBarcode
    
  } else {
    stop('The direcotry is invalid, please input the dir including files: "barcodes.tsv", "features.tsv", "matrix.mtx"')
  }
  
  SingleCellData(assay = list(count = as(run.matrix, 'CsparseMatrix')), rowdata = data.frame(run.gene, stringsAsFactors = FALSE), coldata = data.frame(run.cell, stringsAsFactors = FALSE))
  
}



####################################################################################
#' Import data from 10X Genomics output (h5).
####################################################################################
#' 
#' Import data directly from 10X Genomics output, usually using filtered gene 
#' matrices which contains h5 file. 
#' The user only need to input the directory into "data.path". If not the original 
#' 10X Genomics output, the user can use 'readsc' function.
#' 
#' @rdname Import-10X-h5
#' @param file.path The path of the filtered 10X Genomics output (h5 file).
#' @param is.filter Remove not expressed genes.
#' @importFrom hdf5r H5File
#' @return RISC single cell dataset, including count, coldata, and rowdata.
#' @name read10X_h5
#' @export

read10X_h5 <- function(
  file.path, 
  is.filter = TRUE
) {
  
  if(!exists("file.path")){
    stop('Please input file.path')
  } else {
    file.path = as.character(file.path)
  }
  
  files = H5File$new(filename = file.path, mode = "r")
  key = files$names
  
  if(!is.null(files[[key]])){
    
    if(key == 'matrix'){
      
      run.cell = data.frame(scBarcode = files[['matrix/barcodes']][], row.names = files[['matrix/barcodes']][])
      run.gene = data.frame(Symbol = make.unique(files[['matrix/features/name']][]), Ensembl = files[['matrix/features/id']][], row.names = make.unique(files[['matrix/features/name']][]))
      run.matrix = Matrix::sparseMatrix(dims = files[['matrix/shape']][], x = files[['matrix/data']][], i = files[['matrix/indices']][], p = files[['matrix/indptr']][], index1 = FALSE)
      files$close_all()
      
      run.matrix = as(run.matrix, 'CsparseMatrix')
      run.gene$RNA = "Gene Expression"
      run.gene$nCell = Matrix::rowSums(run.matrix > 0)
      rownames(run.matrix) = rownames(run.gene)
      
      if(is.filter){
        run.gene = run.gene[run.gene$nCell > 0,]
      } else {
        run.gene = run.gene
      }
      run.matrix = run.matrix[rownames(run.matrix) %in% run.gene$Symbol,]
      
      mito.gene = grep(pattern = '^mt-', x = rownames(run.matrix), ignore.case = TRUE, value = TRUE)
      run.cell$UMI = Matrix::colSums(run.matrix)
      run.cell$nGene = Matrix::colSums(run.matrix > 0)
      run.cell$mito = Matrix::colSums(run.matrix[rownames(run.matrix) %in% mito.gene,]) / run.cell$UMI
      colnames(run.matrix) = rownames(run.cell)
      
    } else {
      
      run.cell = data.frame(scBarcode = files[[paste0(key, '/barcodes')]][], row.names = files[[paste0(key, '/barcodes')]][])
      run.gene = data.frame(Symbol = make.unique(files[[paste0(key, '/gene_names')]][]), Ensembl = files[[paste0(key, '/genes')]][], row.names = make.unique(files[[paste0(key, '/gene_names')]][]))
      run.matrix = Matrix::sparseMatrix(dims = files[[paste0(key, '/shape')]][], x = files[[paste0(key, '/data')]][], i = files[[paste0(key, '/indices')]][], p = files[[paste0(key, '/indptr')]][], index1 = FALSE)
      files$close_all()
      
      run.matrix = as(run.matrix, 'CsparseMatrix')
      run.gene$RNA = "Gene Expression"
      run.gene$nCell = Matrix::rowSums(run.matrix > 0)
      rownames(run.matrix) = rownames(run.gene)
      
      if(is.filter){
        run.gene = run.gene[run.gene$nCell > 0,]
      } else {
        run.gene = run.gene
      }
      run.matrix = run.matrix[rownames(run.matrix) %in% run.gene$Symbol,]
      
      mito.gene = grep(pattern = '^mt-', x = rownames(run.matrix), ignore.case = TRUE, value = TRUE)
      run.cell$UMI = Matrix::colSums(run.matrix)
      run.cell$nGene = Matrix::colSums(run.matrix > 0)
      run.cell$mito = Matrix::colSums(run.matrix[rownames(run.matrix) %in% mito.gene,]) / run.cell$UMI
      colnames(run.matrix) = rownames(run.cell)
      
    }
    
  } else {
    stop('The h5 file is invalid, please input correct file path.')
    files$close_all()
  }
  
  SingleCellData(assay = list(count = as(run.matrix, 'CsparseMatrix')), rowdata = data.frame(run.gene, stringsAsFactors = FALSE), coldata = data.frame(run.cell, stringsAsFactors = FALSE))
  
}


