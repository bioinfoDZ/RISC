####################################################################################
#' The Pre-processing Data.
####################################################################################
#' 
#' After input data, RISC preliminarily filter datasets by using three criteria: 
#' first, discard cells with too low/high raw counts/UMIs by distribution analysis. 
#' Second, remove cells with too low expressed genes. Lastly, filter out the genes 
#' only expressed in few cells.
#' 
#' @rdname Filter
#' @param object RISC object: a framework dataset.
#' @param min.UMI The min UMI for valid cells is usually based on the distribution 
#' analysis, with more than 2.5 percentage of UMI distribution, discarding the 
#' cells with too few UMIs. This parameter can be adjusted manually.
#' @param max.UMI The max UMI for valid cells is usually based on the distribution 
#' analysis, with less than 97.5 percentage of UMI distribution, discarding the 
#' cells with too many UMIs.This parameter can be adjusted manually.
#' @param min.gene The min number of expressed genes for valid cells. The default 
#' is based on the distribution analysis, with more than 0.5 percentage of gene 
#' distribution. This parameter can be adjusted manually.
#' @param min.cell The min number of cells for valid expressed genes. The default 
#' is based on the distribution analysis, with more than 0.5 percentage of cell 
#' distribution. This parameter can be adjusted manually.
#' @param mitochon The cutoff of the mitochondrial UMI proportion for valid cells.
#' @param is.filter Whether filter the data.
#' @return RISC single cell dataset, the coldata and rowdata slots.
#' @importFrom stats as.formula dist quantile qnorm
#' @references Liu et al., Nature Biotech. (2021)
#' @name scFilter
#' @export
#' @examples 
#' # RISC object
#' obj0 = raw.mat[[3]]
#' obj0 = scFilter(obj0, min.UMI = 0, max.UMI = Inf, min.gene = 10, min.cell = 3)

scFilter <- function(
  object, 
  min.UMI = NULL, 
  max.UMI = NULL, 
  min.gene = NULL, 
  min.cell = NULL, 
  mitochon = 1, 
  is.filter = TRUE
  ) {
  
  if(is.null(object)) {
    stop('The object is invalid, please input a scdataset.')
  } else {
    
    # size factor
    count = as.matrix(object@assay$count)
    coldata0 = coldata = object@coldata
    coldata$scUMI = Matrix::colSums(count)
    
    if(!is.null(min.UMI) & !is.null(max.UMI) & !is.null(min.cell) & !is.null(min.gene)) {
      
      # filter cells
      if(is.filter){
        coldata = coldata[coldata$scUMI >= min.UMI & coldata$scUMI <= max.UMI,]
        coldata = coldata[coldata$ngene >= min.gene,]
        coldata = coldata[coldata$mito <= mitochon,]
      } else {
        coldata = coldata0
      }
      
      # filter genes
      rowdata = object@rowdata
      count = count[,colnames(count) %in% coldata$scBarcode]
      
      if(is.filter){
        keep = Matrix::rowSums(count > 0) > min.cell
        count = count[keep,]
        rowdata = rowdata[rowdata$Symbol %in% rownames(count),]
        rowdata$nCell = Matrix::rowSums(count > 0)
      } else {
        count = count
        rowdata = rowdata[rowdata$Symbol %in% rownames(count),]
        rowdata$nCell = Matrix::rowSums(count > 0)
      }
      
    } else if(!is.null(object@metadata$filter$min.UMI) & !is.null(object@metadata$filter$max.UMI) & !is.null(object@metadata$filter$min.gene) & !is.null(object@metadata$filter$min.cell)) {
      
      min.UMI = object@metadata$filter$min.UMI
      max.UMI = object@metadata$filter$max.UMI
      min.gene = object@metadata$filter$min.gene
      min.cell = object@metadata$filter$min.cell
      
      # filter cells
      if(is.filter){
        coldata = coldata[coldata$scUMI >= min.UMI & coldata$scUMI <= max.UMI,]
        coldata = coldata[coldata$ngene >= min.gene,]
        coldata = coldata[coldata$mito <= mitochon,]
      } else {
        coldata = coldata0
      }
      
      # filter genes
      rowdata = object@rowdata
      count = count[,colnames(count) %in% coldata$scBarcode]
      
      if(is.filter){
        keep = Matrix::rowSums(count > 0) > min.cell
        count = count[keep,]
        rowdata = rowdata[rowdata$Symbol %in% rownames(count),]
        rowdata$nCell = Matrix::rowSums(count > 0)
      } else {
        count = count
        rowdata = rowdata[rowdata$Symbol %in% rownames(count),]
        rowdata$nCell = Matrix::rowSums(count > 0)
      }
      
    } else if(is.null(min.UMI) & is.null(max.UMI) & is.null(min.cell) & is.null(min.gene)) {
      
      # filter cells
      if(is.filter){
        min.UMI = quantile(coldata$scUMI, probs = seq(0, 1, 0.025))[[2]]
        max.UMI = quantile(coldata$scUMI, probs = seq(0, 1, 0.025))[[40]]
        coldata = coldata[coldata$scUMI >= min.UMI & coldata$scUMI <= max.UMI,]
        min.gene = min(quantile(coldata$ngene, probs = seq(0, 1, 0.005))[[2]], 200)
        coldata = coldata[coldata$ngene >= max(min.gene, 200),]
        coldata = coldata[coldata$mito <= mitochon,]
      } else {
        coldata = coldata0
      }
      
      # filter genes
      rowdata = object@rowdata
      count = count[,colnames(count) %in% coldata$scBarcode]
      
      if(is.filter){
        min.cell = min(floor(ncol(count)/200), 10)
        keep = Matrix::rowSums(count > 0) > min.cell
        count = count[keep,]
        rowdata = rowdata[rowdata$Symbol %in% rownames(count),]
        rowdata$nCell = Matrix::rowSums(count > 0)
      } else {
        count = count
        rowdata = rowdata[rowdata$Symbol %in% rownames(count),]
        rowdata$nCell = Matrix::rowSums(count > 0)
      }
      
    } else {
      stop('Please input all min.UMI, max.UMI, min.cell, mingene, or let the software detect by itself')
    }
  
  count = as(count, 'dgCMatrix')
  object@coldata = coldata
  object@rowdata = rowdata
  object@assay$count = count
  
  if(is.filter){
    object@metadata[['filter']] = data.frame(min.UMI = 0, max.UMI = Inf, min.gene = 0, min.cell = 0, stringsAsFactors = FALSE)
  } else {
    object@metadata[['filter']] = data.frame(min.UMI = min.UMI, max.UMI = max.UMI, min.gene = min.gene, min.cell = min.cell, stringsAsFactors = FALSE)
  }
  
  return(object)
  
  }
  
}



####################################################################################
#' The Processing Data.
####################################################################################
#' 
#' After data filtration, RISC will normalized the raw counts/UMIs by using size 
#' factors which are calculated by the raw counts/UMIs of each cell and will 
#' remove sequencing depth batch. The output will be transformed in log1p. The gene 
#' expression values of RISC object for the subsequent analyses is based on the 
#' normalized data. Here two kinds of normalization can be employed: one is based on 
#' the least absolute deviations, while the other is from the least square. 
#' 
#' @rdname Normalize
#' @param object RISC object: a framework dataset.
#' @param method A method for scdataset normalization, two options: "cosine" and 
#' "robust".
#' @param libsize The standard sum of the UMI in each cell.
#' @param remove.mito Remove mitochondrial genes from library size.
#' @param norm.dis Normalize the distribution of count data.
#' @param ncore The multiple cores for parallel calculating.
#' @return RISC single cell dataset, the assay and rowdata slots.
#' @name scNormalize
#' @importFrom Matrix mean
#' @references Boscovich, R.J. (1757)
#' @references Thompson, W.J., Computers in Physics (1992)
#' @export
#' @examples 
#' # RISC object
#' obj0 = raw.mat[[3]]
#' obj0 = scFilter(obj0, min.UMI = 0, max.UMI = Inf, min.gene = 10, min.cell = 3)
#' obj0 = scNormalize(obj0, method = 'robust')

scNormalize <- function(
  object, 
  method = 'robust', 
  libsize = 1e+06, 
  remove.mito = FALSE, 
  norm.dis = TRUE, 
  ncore = 1
  ){
  
  if(length(object@metadata) == 0 | nrow(object@metadata$filter) == 0){
    stop('Please filter the scdataset first')
  } else {
    
    ncore = as.integer(ncore)
    libsize = as.numeric(libsize)
    count = object@assay$count
    mito.gene = grep(pattern = '^mt-', x = rownames(count), ignore.case = TRUE, value = TRUE)
    coldata = object@coldata
    scUMI0 = NULL
    
    if(norm.dis){
      count = winsorizing2(count)
    } else {
      count = count
    }
    
    if(remove.mito){
      coldata$scUMI = Matrix::colSums(count[!rownames(count) %in% mito.gene,])
    } else {
      coldata$scUMI = Matrix::colSums(count)
    }
    
    if(method == 'cosine') {
      
      scUMI0 = sqrt(Matrix::colSums(count ^ 2))
      coldata$sizefactor = scUMI0 / (libsize * mean(scUMI0 / coldata$scUMI))
      logcount = log1p(count / coldata$sizefactor[col(count)])
      object@metadata[['normalise']] = 'cosine'
      
    } else if(method == "robust") {
      
      coldata$sizefactor = coldata$scUMI / libsize
      logcount = log1p(count / coldata$sizefactor[col(count)])
      object@metadata[['normalise']] = 'robust'
      
    } else {
      stop('Input method: lib.size or all counts')
    }
    
    logcount = as.matrix(logcount)
    if(length(logcount) < 4e+08) {
      logcount = t(apply(logcount, 1, winsorizing1))
      rownames(logcount) = rownames(count)
    } else {
      logcount = pbmcmapply(function(x){winsorizing1(logcount[x,])}, 1L:nrow(logcount), mc.cores = ncore)
      logcount = t(logcount)
      rownames(logcount) = rownames(count)
    }
    
    # keep = Matrix::rowSums(logcount > 0) > 0
    # logcount = logcount[keep,]
    colnames(logcount) = colnames(count)
    object@assay$logcount = as(logcount, 'dgCMatrix')
    object@coldata = coldata
    
  }
  
  return(object)
  
}



####################################################################################
#' Processing Data.
####################################################################################
#' 
#' After data normalization, RISC will perform root-mean-square scaling to the 
#' dataset and generated scaled counts which balance expression levels in each cell 
#' with empirical mean equal to 0. Therefore, only biological signal will be 
#' reserved in scaled counts. RISC utilizes scaled counts for dimension reduction 
#' and data integration.
#' 
#' @rdname Scale
#' @param object RISC object: a framework dataset.
#' @param method A model used for scale scdataset, the default is root-mean-square 
#' scaling.
#' @param center Whether to center the matrix.
#' @param scale Whether using standard deviation to scale the matrix.
#' @return The scaled matrix.
#' @name scScale
#' @importFrom Matrix colMeans
#' @importFrom matrixStats colSds
#' @references Juszczak et al., CiteSeer (2002)
#' @references Jiawei et al. (2011)
#' @export
#' @examples 
#' # RISC object
#' obj0 = raw.mat[[3]]
#' obj0 = scFilter(obj0, min.UMI = 0, max.UMI = Inf, min.gene = 10, min.cell = 3)
#' obj0 = scNormalize(obj0)
#' scale.mat = scScale(obj0)

scScale <- function(object, method = 'scale', center = TRUE, scale = FALSE){
  
  if(length(object@metadata$normalise) == 0) {stop('Please normalize data first')}
  
  else if(method == 'scale'){
    
    count0 = as.matrix(object@assay$logcount)
    if(scale){
      count = colScale(count0, center = center, scale = TRUE)
    } else {
      count = scale(count0, center = TRUE, scale = FALSE)
    }
    
  } else {stop('Waiting for a new method')}
  
  return(as.matrix(count))
  
}



####################################################################################
#' Processing Data.
####################################################################################
#' 
#' After data scaling, RISC will identify highly variably expressed genes, based on 
#' Quasi-Poinsson model, where the coefficient of variation is calculated for each 
#' gene (C > 0.5 as a cutoff for highly variable genes). Then, to controlling for 
#' the relationship between S (standard deviation) and mean (average value), 
#' Quasi-Poisson regression is used to further filter the genes with over-
#' dispersion C caused by small mean. Lastly, RISC estimates the corresponding 
#' ratio with r between the observed C and the predicted C of each gene, with a 
#' threshold r > 1.
#' 
#' @rdname Disperse
#' @param object RISC object: a framework dataset.
#' @param mean.cut A cutoff of the average value of each gene, the genes expression 
#' outside the mean.cut range will be discarded from highly variable genes. The 
#' input is a range vector, like c(0.1, 5).
#' @param ncut The number of fragments using to fit dispersion. The Quasi-Poinsson 
#' model divides genes into a number of bins based on the gene expression. This 
#' parameter indicated how many bins are used for Quasi-Poinsson model.
#' @param top.var The maximum number of highly variable genes, the default is NULL, 
#' including all the highly variable genes.
#' @param method What method is used to define dispersion: quasipoisson or residual.
#' The default method is "Quasi-Poisson".
#' @param ratio To change mean.cut by ratio, and make RISC robustly to select 
#' variable genes. The user can use the default value.
#' @param FDR The cutoff for variable genes.
#' @return RISC single cell dataset, the metadata slot.
#' @name scDisperse
#' @importFrom Matrix rowMeans
#' @importFrom matrixStats rowSds
#' @importFrom stats gaussian glm median p.adjust pchisq quasipoisson poisson sd
#' @references Liu et al., Nature Biotech. (2021)
#' @export
#' @examples 
#' # RISC object
#' obj0 = raw.mat[[3]]
#' obj0 = scFilter(obj0, min.UMI = 0, max.UMI = Inf, min.gene = 10, min.cell = 3)
#' obj0 = scNormalize(obj0)
#' obj0 = scDisperse(obj0)

scDisperse <- function(
  object, 
  method = "quasipoisson", 
  mean.cut = NULL, 
  ncut = 20, 
  top.var = NULL, 
  FDR = 0.1, 
  ratio = 1.25
){
  
  if(length(object@metadata) == 0 | length(object@metadata$normalise) == 0){
    stop('Please filter and normalize the scdataset first')
  } else {
    
    method = as.character(method)
    ncut = as.integer(ncut)
    ratio = as.numeric(ratio)
    FDR0 = as.numeric(FDR)
    
    count = as.matrix(object@assay$logcount)
    keep = rowSums(count > 0) > min(floor(ncol(count)/400), 10)
    count = count[keep,]
    
    sc.mean = abs(rowMeans(count))
    sc.sd = abs(rowSds(count))
    sc.disp = sc.sd/sc.mean
    sc.mean[is.na(sc.mean)] = sc.sd[is.na(sc.sd)] = sc.disp[is.na(sc.disp)] = 0
    
    if(is.null(mean.cut)){
      mean.min = quantile(sc.mean, probs = 0.05)[[1]]
      mean.max = quantile(sc.mean, probs = 0.95)[[1]]
    } else {
      mean.cut = as.numeric(mean.cut)
      mean.min = mean.cut[1]
      mean.max = mean.cut[2]
    }
    
    if(method == "quasipoisson"){
      
      bin = cut(sc.mean, c(-1, quantile(sc.mean, probs = seq(0, 1, 1/ncut))[2:(ncut+1)]), labels = FALSE)
      fit = rep(0, length(sc.mean))
      names(bin) = names(fit) = names(sc.sd) = names(sc.disp) = names(sc.mean)
      
      i = 1L
      while(i <= length(table(bin))){
        bin0 = names(bin)[bin == i]
        fit[bin0] = glm(sc.sd[bin0] ~ sc.mean[bin0], family = quasipoisson)$fitted.values
        i = i + 1L
      }
      
      disp = data.frame(Mean = sc.mean, SD = sc.sd, Dispersion = sc.disp, Disperse.fit = sc.disp/(fit/sc.mean), stringsAsFactors = FALSE)
      disp0 = disp[disp$Dispersion > 0.5 & disp$Disperse.fit > 1 & disp$Mean >= mean.min & disp$Mean <= mean.max,]
      var.gene = rownames(disp0)
      
      if(is.null(top.var)){
        var.gene = var.gene
      } else {
        
        while(length(var.gene) > top.var){
          mean.min = mean.min*ratio
          mean.max = mean.max*ratio
          FDR.cut = FDR0
          disp1 = disp[disp$Dispersion > 0.5 & disp$Disperse.fit > 1 & disp$Mean >= mean.min & disp$Mean <= mean.max,]
          pval = pchisq(disp1$Disperse.fit*(length(disp1$Disperse.fit)-1), (length(disp1$Disperse.fit)-1), lower.tail = FALSE)
          FDR = p.adjust(pval, method = 'BH')
          names(FDR) = names(pval) = rownames(disp1)
          FDR = FDR[order(FDR, decreasing = FALSE)]
          var.gene = names(FDR)[FDR < FDR.cut]
        }
        
      }
      
    } else {
      stop('Input method quasipoisson or residual')
    }
    
    object@metadata[['dispersion']] = disp
    object@metadata[['dispersion.var']] = disp0[order(disp0$Disperse.fit, decreasing = TRUE),]
    object@vargene = var.gene
    return(object)
    
  }
  
}



####################################################################################
####################################################################################
winsorizing1 <- function(x){
  x = as.matrix(x)
  lmax = quantile(x, probs = (1 - 5/length(x)))
  x[x > lmax] = lmax
  return(x)
}

winsorizing2 <- function(x){
  x = as.matrix(x)
  lmax = qnorm(seq(0, 1, length.out = length(x) / 2), mean = mean(x), sd = sd(x), lower.tail = FALSE)[2]
  x[x > lmax] = ceiling(lmax)
  return(x)
}

colScale = function(x, center = TRUE, scale = TRUE){
  colmean = colMeans(x, na.rm = TRUE)
  if(scale){
    colsd = colSds(x, center = colmean)
  } else {
    colsd = rep(1, length(colmean))
  }
  x = t((t(x) - colmean) / colsd)
  return(x)
}


