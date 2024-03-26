####################################################################################
#' UMAP Plots
####################################################################################
#' 
#' The UMAP plots are widespread in scRNA-seq data analysis. Here, the "UMAPlot" 
#' function not only can make plots for factor labels of individual cells but also 
#' can show gene expression values of each cell.
#' 
#' @rdname UMAPlot
#' @param object RISC object: a framework dataset.
#' @param colFactor Use the factor (column name) in the coldata to make a UMAP plot, 
#' but each time only one column name can be inputted.
#' @param genes Use the gene expression values (gene symbol) to make UMAP plots, 
#' each time more than one genes can be inputted.
#' @param legend Whether a legend shown at UMAP plot.
#' @param Colors The users can use their own colors (color vector). The default of 
#' the "UMAPlot" funciton will assign colors automatically.
#' @param size Choose the size of dots at UMAP plots, the default size is 0.5.
#' @param Alpha Whether show transparency of individual points, the default is 0.8.
#' @param plot.ncol If the users input more than one genes, the arrangement of 
#' multiple UMAP plots depends on this parameter.
#' @param raw.count If use normalized or raw counts.
#' @param exp.range The gene expression cutoff for plot, e.g. "c(0, 1.5)" for 
#' expression level between 0 and 1.5.
#' @param exp.col The gradient color for gene expression.
#' @import ggplot2
#' @importFrom gridExtra grid.arrange
#' @import RColorBrewer
#' @importFrom grDevices col2rgb colorRampPalette
#' @references Wickham, H. (2016)
#' @references Auguie, B. (2015)
#' @name UMAPlot

UMAPlot <- function(
  object, 
  colFactor = NULL, 
  genes = NULL, 
  legend = TRUE, 
  Colors = NULL, 
  size = 0.5, 
  Alpha = 0.8, 
  plot.ncol = NULL,
  raw.count = FALSE, 
  exp.range = NULL, 
  exp.col = "firebrick2"
  ) {
  
  X = Y = as.character()
  exp.col = as.character(exp.col)
  size0 = as.numeric(size)
  
  if(is.null(object@DimReduction$cell.umap)){
    stop("Please run scUMAP first")
  } else {
    
    m0 = data.frame(X = object@DimReduction$cell.umap[,1], Y = object@DimReduction$cell.umap[,2])
    draw = NULL
    
    if(is.null(colFactor) & is.null(genes)){
      stop('Please input either colFactor or gene symbol')
    } else if(!is.null(colFactor) & is.null(genes)){
      
      colFactor = as.character(colFactor)
      colData = object@coldata
      m0$factor = as.factor(colData[,colFactor])
      
      if(is.null(Colors)){
        Cols0 = Color0(length(levels(m0$factor)))
      } else {
        Cols0 = Colors
      }
      
      if(length(levels(m0$factor)) <= 35 | !is.null(Colors)){
        plot = ggplot(m0, aes(X, Y)) + geom_point(aes(color = as.factor(factor)), alpha = Alpha, size = size0) + scale_color_manual(values = Cols0) + labs(x = 'UMAP-1', y = 'UMAP-2', color = colFactor) + theme_bw(base_size = 12, base_line_size = 0) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + guides(colour = guide_legend(override.aes = list(size = 4), ncol = ifelse(length(levels(m0$factor)) <= 15, 1, 2)))
      } else {
        plot = ggplot(m0, aes(X, Y)) + geom_point(aes(color = as.factor(factor)), alpha = Alpha, size = size0) + labs(x = 'UMAP-1', y = 'UMAP-2', color = colFactor) + theme_bw(base_size = 12, base_line_size = 0) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + guides(colour = guide_legend(override.aes = list(size = 4), ncol = ifelse(length(levels(m0$factor)) <= 15, 1, 2)))
      }
      return(plot)
      
    } else if (is.null(colFactor) & !is.null(genes)){
      
      if(raw.count == FALSE) {
        count = as.matrix(object@assay$logcount)
      } else {
        count = as.matrix(object@assay$count)
      }
      genes = intersect(genes, rownames(count))
      
      if(length(genes) == 1){
        
        count = count[genes,]
        
        if(is.null(exp.range)){
          m0$draw = as.numeric(count)
        } else {
          m0$draw = as.numeric(count)
          min0 = exp.range[1]
          max0 = exp.range[2]
          m0$draw[m0$draw < min0] = min0
          m0$draw[m0$draw > max0] = max0
        }
        
        if(isTRUE(legend)){
          plot = ggplot(m0, aes(X, Y)) + geom_point(aes(color = draw), alpha = Alpha, size = size0) + scale_color_gradient(low = 'gray80', high = exp.col) + labs(x = 'UMAP-1', y = 'UMAP-2', color = 'Expre.', title = genes) + theme_bw(base_size = 12, base_line_size = 0) + theme(plot.title = element_text(hjust = 0.95, size = 16))
        } else {
          plot = ggplot(m0, aes(X, Y)) + geom_point(aes(color = draw), alpha = Alpha, size = size0) + scale_color_gradient(low = 'gray80', high = exp.col, guide = 'none') + labs(x = 'UMAP-1', y = 'UMAP-2', color = 'Expre.', title = genes) + theme_bw(base_size = 12, base_line_size = 0) + theme(plot.title = element_text(hjust = 0.95, size = 16))
        }
        
        return(plot)
        
      } else{
        
        count = count[genes,]
        
        if(raw.count == FALSE) {
          count = t(apply(count, 1, winsorizing))
        } else {
          count = count
        }
        
        if(is.null(exp.range)){
          count = count
        } else {
          min0 = exp.range[1]
          max0 = exp.range[2]
          count[count < min0] = min0
          count[count > max0] = max0
        }
        
        if(isTRUE(legend)){
          plot0 = lapply(genes, function(z){ggplot(m0, aes(X, Y)) + geom_point(aes(color = count[z,]), alpha = Alpha, size = size0) + scale_color_gradient(low = 'gray80', high = exp.col) + labs(x = 'UMAP-1', y = 'UMAP-2', color = 'Expre.', title = z) + theme_bw(base_size = 12, base_line_size = 0) + theme(plot.title = element_text(hjust = 0.95, size = 16))})
        } else {
          plot0 = lapply(genes, function(z){ggplot(m0, aes(X, Y)) + geom_point(aes(color = count[z,]), alpha = Alpha, size = size0) + scale_color_gradient(low = 'gray80', high = exp.col, guide = 'none') + labs(x = 'UMAP-1', y = 'UMAP-2', color = 'Expre.', title = z) + theme_bw(base_size = 12, base_line_size = 0) + theme(plot.title = element_text(hjust = 0.95, size = 16))})
        }
        
        names(plot0) = genes
        
        if(is.null(plot.ncol)){
          return(do.call(grid.arrange, c(plot0, ncol = ceiling(sqrt(length(genes))))))
        } else {
          return(do.call(grid.arrange, c(plot0, ncol = plot.ncol)))
        }
        
      }
      
    } else {
      stop('Please input colFactor or gene symbol')
    }
    
  }
  
}



####################################################################################
#' Dimension Reduction Plots
####################################################################################
#' 
#' The Dimension Reduction plots are widespread in scRNA-seq data analysis. Here, 
#' the "DimPlot" function not only can make plots for factor labels of individual 
#' cells but also can show gene expression values of each cell.
#' 
#' @rdname DimPlot
#' @param object RISC object: a framework dataset.
#' @param slot The dimension_reduction slot for drawing the plots. The default is 
#' "cell.umap" under RISC object "DimReduction" item for UMAP plot, but the customer 
#' can add new dimension_reduction method under DimReduction and use it.
#' @param colFactor Use the factor (column name) in the coldata to make a 
#' dimension_reduction plot, but each time only one column name can be inputted.
#' @param genes Use the gene expression values (gene symbol) to make dimension 
#' reduction plot, each time more than one genes can be inputted.
#' @param legend Whether a legend shown at dimension_reduction plot.
#' @param Colors The users can use their own colors (color vector). The default of 
#' the "tSNEPlot" funciton will assign colors automatically.
#' @param size Choose the size of dots at dimension_reduction plot, the default 
#' size is 0.5.
#' @param Alpha Whether show transparency of individual points, the default is 0.8.
#' @param plot.ncol If the users input more than one genes, the arrangement of 
#' multiple dimension_reduction plot depends on this parameter.
#' @param exp.range The gene expression cutoff for plot, e.g. "c(0, 1.5)" for 
#' expression level between 0 and 1.5.
#' @param exp.col The gradient color for gene expression.
#' @param label Whether label the clusters or cell populations in the plot.
#' @param adjust.label The adjustment of the label position.
#' @param label.font The font size for the label.
#' @import ggplot2
#' @importFrom gridExtra grid.arrange
#' @import RColorBrewer
#' @importFrom grDevices col2rgb colorRampPalette
#' @references Wickham, H. (2016)
#' @references Auguie, B. (2015)
#' @name tSNEPlot
#' @export
#' @examples 
#' # RISC object
#' obj0 = raw.mat[[3]]
#' obj0 = scPCA(obj0, npc = 10)
#' obj0 = scUMAP(obj0, npc = 3)
#' DimPlot(obj0, slot = "cell.umap", colFactor = 'Group', size = 2, label = TRUE)
#' DimPlot(obj0, genes = c('Gene718', 'Gene325', 'Gene604'), size = 2)

DimPlot <- function(
  object, 
  slot = "cell.umap", 
  colFactor = NULL, 
  genes = NULL, 
  legend = TRUE, 
  Colors = NULL, 
  size = 0.5, 
  Alpha = 0.8, 
  plot.ncol = NULL,
  exp.range = NULL, 
  exp.col = "firebrick2", 
  label = FALSE, 
  adjust.label = 0.25, 
  label.font = 5
  ) {
  
  X = Y = label0 = as.character()
  exp.col = as.character(exp.col)
  slot = as.character(slot)
  adjust0 = as.numeric(adjust.label)
  font0 = as.integer(label.font)
  dimReduce0 = object@DimReduction[[slot]]
  size0 = as.numeric(size)
  
  if(is.null(dimReduce0)){
    stop("Do not include this dimention_reduction slot, try another one")
  } else {
    
    m0 = data.frame(X = dimReduce0[,1], Y = dimReduce0[,2])
    draw = NULL
    
    if(is.null(colFactor) & is.null(genes)){
      stop('Please input either colFactor or gene symbol')
    } else if(!is.null(colFactor) & is.null(genes)){
      
      colFactor = as.character(colFactor)
      colData = object@coldata
      m0$factor = as.factor(colData[,colFactor])
      
      if(is.null(Colors)){
        Cols0 = Color0(length(levels(m0$factor)))
      } else {
        Cols0 = Colors
      }
      
      if(length(size0) == length(levels(m0$factor))){
        size0 = size0
      } else {
        size0 = rep(size0[1], length(levels(m0$factor)))
      }
      
      if(length(levels(m0$factor)) <= 35 | !is.null(Colors)){
        if(isTRUE(label)){
          m1 = data.frame(X = rep(0, length(levels(m0$factor))), Y = rep(0, length(levels(m0$factor))), label0 = levels(m0$factor))
          m1$X = aggregate(m0$X, list(m0$factor), mean)[,2]
          m1$Y = aggregate(m0$Y, list(m0$factor), mean)[,2]
          plot = ggplot(m0, aes(X, Y)) + geom_point(aes(color = as.factor(factor), size = as.factor(factor)), alpha = Alpha) + geom_text(data = m1, aes(X, Y, label = label0), nudge_x = adjust0, nudge_y = adjust0, size = font0) + scale_color_manual(values = Cols0, na.value = NA) + scale_size_manual(values = size0) + labs(x = 'Dimension-1', y = 'Dimension-2', color = colFactor, size = colFactor) + theme_bw(base_size = 12, base_line_size = 0) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + guides(colour = guide_legend(override.aes = list(size = 4), ncol = ifelse(length(levels(m0$factor)) <= 15, 1, 2)))
        } else {
          plot = ggplot(m0, aes(X, Y)) + geom_point(aes(color = as.factor(factor), size = as.factor(factor)), alpha = Alpha) + scale_color_manual(values = Cols0, na.value = NA) + scale_size_manual(values = size0) + labs(x = 'Dimension-1', y = 'Dimension-2', color = colFactor, size = colFactor) + theme_bw(base_size = 12, base_line_size = 0) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + guides(colour = guide_legend(override.aes = list(size = 4), ncol = ifelse(length(levels(m0$factor)) <= 15, 1, 2)))
        }
      } else {
        if(isTRUE(label)){
          m1 = data.frame(X = rep(0, length(levels(m0$factor))), Y = rep(0, length(levels(m0$factor))), label0 = levels(m0$factor))
          m1$X = aggregate(m0$X, list(m0$factor), mean)[,2]
          m1$Y = aggregate(m0$Y, list(m0$factor), mean)[,2]
          plot = ggplot(m0, aes(X, Y)) + geom_point(aes(color = as.factor(factor), size = as.factor(factor)), alpha = Alpha) + geom_text(data = m1, aes(X, Y, label = label0), nudge_x = adjust0, nudge_y = adjust0, size = font0) + scale_size_manual(values = size0) + labs(x = 'Dimension-1', y = 'Dimension-2', color = colFactor, size = colFactor) + theme_bw(base_size = 12, base_line_size = 0) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + guides(colour = guide_legend(override.aes = list(size = 4), ncol = ifelse(length(levels(m0$factor)) <= 15, 1, 2)))
        } else {
          plot = ggplot(m0, aes(X, Y)) + geom_point(aes(color = as.factor(factor), size = as.factor(factor)), alpha = Alpha) + scale_size_manual(values = size0) + labs(x = 'Dimension-1', y = 'Dimension-2', color = colFactor, size = colFactor) + theme_bw(base_size = 12, base_line_size = 0) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + guides(colour = guide_legend(override.aes = list(size = 4), ncol = ifelse(length(levels(m0$factor)) <= 15, 1, 2)))
        }
      }
      return(plot)
      
    } else if (is.null(colFactor) & !is.null(genes)){
      
      genes = intersect(genes, rownames(object@rowdata))
      
      if('Integration' %in% names(object@metadata)){
        
        count = object@assay$logcount
        count = lapply(count, FUN = function(y){y[genes,]})
        if(length(genes) > 1){
          count = do.call(cbind, count)
        } else if(length(genes) == 1){
          count = Matrix::Matrix(unlist(count), nrow = 1)
        } else {
          stop("Input valid gene symbol")
        }

      } else {
        count = object@assay$logcount[genes, , drop = F]
      }
      
      count = winsorizing(count)
      
      if(length(genes) == 1){
        
        if(is.null(exp.range)){
          m0$draw = as.numeric(count)
        } else {
          m0$draw = as.numeric(count)
          min0 = exp.range[1]
          max0 = exp.range[2]
          m0$draw[m0$draw < min0] = min0
          m0$draw[m0$draw > max0] = max0
        }
        
        if(isTRUE(legend)){
          plot = ggplot(m0, aes(X, Y)) + geom_point(aes(color = draw), alpha = Alpha, size = size0) + scale_color_gradient(low = 'gray80', high = exp.col) + labs(x = 'Dimension-1', y = 'Dimension-2', color = 'Expre.', title = genes) + theme_bw(base_size = 12, base_line_size = 0) + theme(plot.title = element_text(hjust = 0.95, size = 16))
        } else {
          plot = ggplot(m0, aes(X, Y)) + geom_point(aes(color = draw), alpha = Alpha, size = size0) + scale_color_gradient(low = 'gray80', high = exp.col, guide = 'none') + labs(x = 'Dimension-1', y = 'Dimension-2', color = 'Expre.', title = genes) + theme_bw(base_size = 12, base_line_size = 0) + theme(plot.title = element_text(hjust = 0.95, size = 16))
        }
        
        return(plot)
        
      } else {
        
        if(is.null(exp.range)){
          count = count
        } else {
          min0 = exp.range[1]
          max0 = exp.range[2]
          count@x[count@x < min0] = min0
          count@x[count@x > max0] = max0
        }
        
        if(isTRUE(legend)){
          plot0 = lapply(genes, function(z){ggplot(m0, aes(X, Y)) + geom_point(aes(color = count[z,]), alpha = Alpha, size = size0) + scale_color_gradient(low = 'gray80', high = exp.col) + labs(x = 'Dimension-1', y = 'Dimension-2', color = 'Expre.', title = z) + theme_bw(base_size = 12, base_line_size = 0) + theme(plot.title = element_text(hjust = 0.95, size = 16))})
        } else {
          plot0 = lapply(genes, function(z){ggplot(m0, aes(X, Y)) + geom_point(aes(color = count[z,]), alpha = Alpha, size = size0) + scale_color_gradient(low = 'gray80', high = exp.col, guide = 'none') + labs(x = 'Dimension-1', y = 'Dimension-2', color = 'Expre.', title = z) + theme_bw(base_size = 12, base_line_size = 0) + theme(plot.title = element_text(hjust = 0.95, size = 16))})
        }
        
        names(plot0) = genes
        
        if(is.null(plot.ncol)){
          return(do.call(grid.arrange, c(plot0, ncol = ceiling(sqrt(length(genes))))))
        } else {
          return(do.call(grid.arrange, c(plot0, ncol = plot.ncol)))
        }
        
      }
      
    } else {
      stop('Please input colFactor or gene symbol')
    }
    
  }
  
}



####################################################################################
#' Processing Plot
####################################################################################
#' 
#' The "FilterPlot" function makes the plots to show the UMIs and expressed genes 
#' of individual cells. These plots are usually used to estimate the data before 
#' and after pre-processing, so the users can visually select the optimal 
#' parameters to filter data.
#' 
#' @rdname FilterPlot
#' @param object RISC object: a framework dataset.
#' @param colFactor Use the factor (column name) in the coldata to make a processing 
#' plot, but each time only one column name can be inputted.
#' @references Wickham, H. (2016)
#' @references Auguie, B. (2015)
#' @references Liu et al., Nature Biotech. (2021)
#' @name FilterPlot
#' @export
#' @examples 
#' # RISC object
#' obj0 = raw.mat[[3]]
#' FilterPlot(obj0, colFactor = 'Group')

FilterPlot <- function(
  object, 
  colFactor = NULL
  ) {
  
  UMI =  nGene = as.character()
  coldata = as.data.frame(object@coldata)
  
  if(is.null(colFactor)){
    
    plot1 = ggplot(coldata, aes(UMI,  nGene)) + geom_point(color = '#80B1D3', size = 0.5, alpha = 0.8) + labs(x = 'UMI', y = ' nGene', color = '', title = ' nGene ~ UMI') + theme_bw(base_size = 12, base_line_size = 0.1) + theme(plot.title = element_text(hjust = 1, size = 16))
    return(plot1)
    
  } else {
    
    colFactor = as.character(colFactor)
    plot1 = ggplot(coldata, aes(UMI,  nGene)) + geom_point(color = '#FB8072', size = 0.5, alpha = 0.8) + labs(x = 'UMI', y = ' nGene', color = '', title = ' nGene ~ UMI') + theme(plot.title = element_text(hjust = 1, size = 16))
    plot2 = ggplot(coldata, aes(coldata[,colFactor],  nGene)) + geom_violin(aes(fill = coldata[,colFactor]), scale = 'width') + guides(fill = 'none') + labs(x = '', y = '', fill = '', title = ' nGene') + theme(plot.title = element_text(hjust = 1, size = 16))
    plot3 = ggplot(coldata, aes(coldata[,colFactor], UMI)) + geom_violin(aes(fill = coldata[,colFactor]), scale = 'width') + guides(fill = 'none') + labs(x = '', y = '', fill = '', title = 'UMI') + theme(plot.title = element_text(hjust = 1, size = 16))
    return(grid.arrange(plot1, plot2, plot3, ncol = 3))
    
  }
  
}



####################################################################################
#' Processing Plot
####################################################################################
#' 
#' The "PCPlot" function makes the plot to show How the PCs explain the variance. 
#' This plot helps the users to select the optimal PCs to perform dimension 
#' reduction and data integration. 
#' 
#' @rdname PCPlot
#' @param object RISC object: a framework dataset.
#' @references Wickham, H. (2016)
#' @references Auguie, B. (2015)
#' @references Liu et al., Nature Biotech. (2021)
#' @name PCPlot

PCPlot <- function(object) {
  
  PC = Variance = as.character()
  PC0 = object@DimReduction$var.pca
  PC0c = cumsum(PC0)
  
  PCs = data.frame(PC = 1:length(PC0), Variance = PC0)
  PCcs = data.frame(PC = 1:length(PC0c), Variance = PC0c)
  plot1 = ggplot(PCs, aes(PC, Variance)) + geom_point(color = 'firebrick1', size = 3, alpha = 0.8) + geom_line(color = 'navy', size = 1) + labs(x = 'PCs', y = 'Variance', color = '', title = 'The PCs vs variance') + theme_bw(base_size = 12, base_line_size = 0.1) + theme(plot.title = element_text(hjust = 1, size = 16), axis.text.x = element_text(angle = 90, vjust = 0.5))
  plot2 = ggplot(PCs, aes(PC, Variance)) + geom_point(color = 'firebrick1', size = 3, alpha = 0.8) + labs(x = 'PCs', y = 'Variance', color = '') + ylim(0, 0.2) + theme_bw(base_size = 12, base_line_size = 0.1) + theme(plot.title = element_text(hjust = 1, size = 16), axis.text.x = element_text(angle = 90, vjust = 0.5))
  plot3 = ggplot(PCs, aes(PC, Variance)) + geom_point(color = 'firebrick1', size = 3, alpha = 0.8) + labs(x = 'PCs', y = 'Variance', color = '') + ylim(0, 0.1) + theme_bw(base_size = 12, base_line_size = 0.1) + theme(plot.title = element_text(hjust = 1, size = 16), axis.text.x = element_text(angle = 90, vjust = 0.5))
  plot4 = ggplot(PCs, aes(PC, Variance)) + geom_point(color = 'firebrick1', size = 3, alpha = 0.8) + labs(x = 'PCs', y = 'Variance', color = '') + ylim(0, 0.05) + theme_bw(base_size = 12, base_line_size = 0.1) + theme(plot.title = element_text(hjust = 1, size = 16), axis.text.x = element_text(angle = 90, vjust = 0.5))
  return(grid.arrange(plot1, plot2, plot3, plot4, ncol = 2))
  
}



####################################################################################
#' Heatmap
####################################################################################
#' 
#' The "Heat" map makes heatmap to show gene expression patterns of single cells. 
#' The default groups cells into clusters, so the column of heatmap represents 
#' genes while the row of heatmap for the clusters of all the cells.
#' 
#' @rdname Heatmap
#' @param object RISC object: a framework dataset.
#' @param colFactor Use the factor (column name) in the coldata to make heatmap, but 
#' be factors.
#' @param genes Use the gene expression values (gene symbol) to make heatmap, need to 
#' be inputted by the users.
#' @param cells Use the subset cells of the whole coldata (cells) to make heatmap, 
#' the default is NULL and including all the cells.
#' @param gene.lab Whether label gene names for the heatmap.
#' @param gene.cluster The cluster numbers for gene clustering in the heatmap. 
#' The default is 0, without clustering genes.
#' @param sample_bin The cell aggregating in samples, the default is FALSE.
#' @param ann_col The annotation colors for colFactors, the input is a list.
#' @param lim The gene expression range shown at heat-maps.
#' @param smooth If use smooth to adjust heatmap, the default is "smooth" and another
#' choice is "loess".
#' @param span The loess span.
#' @param degree The loess degree.
#' @param palette The color palette used for heatmap. The default is 
#' brewer.pal(n = 7, name = "RdYlBu").
#' @param num The cells for individual bin spans.
#' @param con.bin Whether use consistent bin span.
#' @param cell.lab.size The font size for column.
#' @param gene.lab.size The font size for row.
#' @param value_only Only return values.
#' @importFrom pheatmap pheatmap
#' @importFrom stats aggregate cutree smooth loess.smooth
#' @references Kolde, R. (2015)
#' @name Heat
#' @export
#' @examples 
#' # RISC object
#' obj0 = raw.mat[[3]]
#' gene0 = c('Gene718', 'Gene120', 'Gene313', 'Gene157', 'Gene30', 
#'           'Gene325', 'Gene415', 'Gene566', 'Gene990', 'Gene13', 
#'           'Gene604', 'Gene934', 'Gene231', 'Gene782', 'Gene10')
#' Heat(obj0, colFactor = 'Group', genes = gene0, gene.lab = TRUE, gene.cluster = 3, 
#' sample_bin = TRUE, lim = 2, gene.lab.size = 8)

Heat <- function(
  object, 
  colFactor = NULL, 
  genes = NULL, 
  cells = NULL, 
  gene.lab = FALSE, 
  gene.cluster = 0, 
  sample_bin = FALSE, 
  ann_col = NULL, 
  lim = NULL, 
  smooth = "smooth", 
  span = 0.75,
  degree = 1,
  palette = NULL, 
  num = 50, 
  con.bin = TRUE, 
  cell.lab.size = 10, 
  gene.lab.size = 5,
  value_only = FALSE,
  ...
  ) {
  
  if(length(colFactor) == 0 | length(genes) == 0){
    stop('Please input both colFactor and genes')
  } else if(!is.null(cells)) {
    colFactor = as.character(colFactor)
    genes = as.character(genes)
    coldata = as.data.frame(object@coldata)
    coldata = coldata[rownames(coldata) %in% cells,]
    colFactor0 = colFactor[colFactor %in% colnames(coldata)]
  } else {
    colFactor = as.character(colFactor)
    genes = as.character(genes)
    coldata = as.data.frame(object@coldata)
    colFactor0 = colFactor[colFactor %in% colnames(coldata)]
  }
  
  if(is.null(palette)){
    col0 = rev(c("#D73027", "#FC8D59", "#FEE090", "#FFFFBF", "#E0F3F8", "#91BFDB", "#4575B4"))
  } else {
    col0 = palette
  }
  
  if(is.character(smooth)){
    smooth = smooth
  } else {
    smooth = "none"
  }
  
  mat0 = object@assay$logcount
  genes = intersect(genes, rownames(object@rowdata))
  cell0 = rownames(coldata)
  
  if('Integration' %in% names(object@metadata)){
    mat0 = lapply(mat0, FUN = function(y){y[genes, colnames(y) %in% cell0]})
    mat0 = do.call(cbind, mat0)
  } else {
    mat0 = mat0[genes, cell0]
  }
  
  # mat0 = winsorizing(mat0)
  mat0 = as.matrix(mat0)
  
  if(length(colFactor0) == 0){
    stop('Please input valid colFactor')
  } else if(length(colFactor0) == 1){
    Group0 = data.frame(G1 = coldata[,colFactor0], G2 = coldata[,colFactor0])
    colnames(Group0) = c(colFactor0, 'Seq')
    colFactor = c(colFactor0, 'Seq')
  } else {
    Group0 = coldata[,colFactor0]
    for(i in 1L:length(colFactor0)){
      Group0[,i] = as.character(Group0[,i])
      Group0[,i] = as.factor(Group0[,i])
    }
    colFactor = colFactor0
  }
  
  All0 = data.frame(Group0, t(mat0))
  Order = do.call(order, c(data.frame(All0[,colFactor], All0[,!colnames(All0) %in% colFactor])))
  All0 = All0[Order,]
  Group0 = All0[,colFactor]
  
  ######
  Bin2 = 0
  i = 1L
  num = as.numeric(num)
  account = 0
  
  while(i <= length(levels(Group0[,1]))){
    
    Group1 = Group0[Group0[,1] == levels(Group0[,1])[i],]
    if(sample_bin == TRUE){
      num = table(Group1[,-1])
      Bin0 = rep((account + 1L):(account + length(num)), num)
    } else if(sample_bin == FALSE & con.bin == TRUE) {
      Bin0 = c(rep((account + 1L):(account + (num - 1L)), each = floor(nrow(Group1)/num)), rep((account + num), (nrow(Group1) - floor(nrow(Group1)/num) * (num - 1L))))
    } else {
      Bin0 = c(rep((account + 1L):(account + floor(nrow(Group1)/num) - 1L), each = num), rep((account + floor(nrow(Group1)/num)), (nrow(Group1) - (floor(nrow(Group1)/num) - 1L) * num)))
    }
    Bin2 = c(Bin2, Bin0)
    account = (account + max(Bin2))
    i = i + 1L
    
  }
  #####
  
  Bin = Bin2[-1]
  Group0$Bin = Bin
  mat.bin0 = data.frame(Bin = Bin, All0[,!colnames(All0) %in% colFactor])
  mat.bin0$Bin = as.factor(mat.bin0$Bin)
  mat.bin1 = aggregate(mat.bin0[,2L:ncol(mat.bin0)], list(mat.bin0$Bin), mean)
  mat2 = as.matrix(mat.bin1[,-1])
  Group1 = Group0[!duplicated(Group0$Bin),]
  Group0 = Group1[,-ncol(Group1)]
  rownames(mat2) = rownames(Group0)
  # cmin0 = as.matrix(table(Group0))
  # keep = colSums(mat2 > 0) > min(cmin0[cmin0 > 0])
  keep = colSums(mat2 > 0) > 0
  mat.sel = mat2[,keep]
  
  if(smooth == "loess" & sample_bin == FALSE){
    group = 1:nrow(mat.sel)
    mat.sel0 = apply(mat.sel, 2, FUN = function(x){loess.smooth(group, x, span = span, degree = degree, evaluation = nrow(mat.sel))$y})
    rownames(mat.sel0) = rownames(mat.sel)
    colnames(mat.sel0) = colnames(mat.sel)
    keep = colSums(mat.sel0 > 0) > 0
    mat.sel = mat.sel0[,keep]
  } else if(smooth == "smooth" & sample_bin == FALSE) {
    mat.sel0 = apply(mat.sel, 2, FUN = function(x){smooth(x, twiceit = TRUE, kind = '3RS3R')})
    rownames(mat.sel0) = rownames(mat.sel)
    colnames(mat.sel0) = colnames(mat.sel)
    keep = colSums(mat.sel0 > 0) > 0
    mat.sel = mat.sel0[,keep]
  } else {
    mat.sel = mat.sel
  }
  
  mat3 = scale(mat.sel)
  mat3[is.na(mat3)] = 0
  
  if(length(colFactor0) > 1){
    # Group1 = Group0
    Group1 = Group0[,seq(dim(Group0)[2], 1)]
    ann_col0 = ann_col
  } else {
    Group1 = data.frame(Group = Group0[,1], row.names = rownames(Group0))
    if(is.null(ann_col)){
      ann_col0 = NULL
    } else if(nrow(summary(ann_col)) == 1) {
      names(ann_col) = 'Group'
      ann_col0 = ann_col
    } else {
      ann_col0 = NULL
    }
  }
  
  if(is.null(lim)){
    lim0 = 3.1
  } else {
    lim0 = as.numeric(lim)
  }
  
  mat3[mat3 > lim0] = lim0
  mat3[mat3 < -lim0] = -lim0
  bks = seq(-lim0, lim0, by = 0.1)
  
  if(is.null(palette)){
    cols = colorRampPalette(rev(c("#D73027", "#F46D43", "#FDAE61", "#FEE090", "#FFFFBF", "#E0F3F8", "#ABD9E9", "#74ADD1", "#4575B4")))(length(bks) - 1)
  } else {
    cols = colorRampPalette(palette)(length(bks) - 1)
  }
  
  gene.cluster = as.integer(gene.cluster)
  
  if(value_only){
    return(mat3)
  } else {
    
    if(gene.cluster > 1){
      
      out = pheatmap(t(mat3), scale = 'none', useRaster = TRUE, cluster_cols = FALSE, cluster_rows = TRUE, treeheight_row = 0, treeheight_col = 0, show_rownames = gene.lab, show_colnames = FALSE, breaks = bks, color = cols, annotation_col = Group1, silent = TRUE)
      gene.cluster = cutree(out$tree_row, k = gene.cluster)
      gene.cluster0 = data.frame(Gene_K = as.factor(gene.cluster))
      rownames(gene.cluster0) = names(gene.cluster)
      
      out0 = pheatmap(t(mat3), scale = 'none', useRaster = TRUE, cluster_cols = FALSE, cluster_rows = TRUE, treeheight_row = 0, treeheight_col = 0, show_rownames = gene.lab, show_colnames = FALSE, breaks = bks, color = cols, annotation_col = Group1, annotation_row = gene.cluster0, annotation_colors = ann_col0, border_color = NA, fontsize_row = gene.lab.size, fontsize_col = cell.lab.size,  ... = ...)
      
    } else {
      
      out0 = pheatmap(t(mat3), scale = 'none', useRaster = TRUE, cluster_cols = FALSE, cluster_rows = FALSE, treeheight_row = 0, treeheight_col = 0, show_rownames = gene.lab, show_colnames = FALSE, breaks = bks, color = cols, annotation_col = Group1, annotation_colors = ann_col0, border_color = NA, fontsize_row = gene.lab.size, fontsize_col = cell.lab.size, ... = ...)
      
    }
    
    return(out0)
    
  }
  
}



####################################################################################
#' Violin Plot
####################################################################################
#' 
#' The "ViolinPlot" map makes plots to show gene expression patterns of the clusters 
#' or other factors. The default groups cells into the clusters, but the users can 
#' input the factor (column name) of coldata.  
#' 
#' @rdname Violin-Plot
#' @param object RISC object: a framework dataset.
#' @param colFactor Use the factor (column name) in the coldata to make heatmap, but 
#' each time only one column name can be inputted.
#' @param genes The gene expression pattern: gene symbol
#' @param legend Whether a legend shown at heatmap.
#' @param trim Whether trim the violin plot
#' @param Colors The users can use their own colors (color vector). The default of 
#' the "UMAPlot" funciton will assign colors automatically.
#' @param Alpha Whether show transparency of individual points, the default is 0.8.
#' @param dots Adding jitter dots to the violin plot. The default is TRUE.
#' @param wid The scale format, options: "area", "width", "count". 
#' The default: "area"
#' @references Wickham, H. (2016)
#' @references Auguie, B. (2015)
#' @name ViolinPlot
#' @export
#' @examples 
#' # RISC object
#' obj0 = raw.mat[[3]]
#' ViolinPlot(obj0, colFactor = 'Group', genes = 'Gene718')

ViolinPlot <- function(
  object, 
  colFactor = NULL, 
  genes = NULL, 
  legend = TRUE, 
  trim = TRUE, 
  Colors = NULL, 
  Alpha = 0.8, 
  dots = TRUE, 
  wid = "area"
  ) {
  
  Factor = NULL
  if(is.null(genes)){
    stop("Please input one gene symbol")
  } else {
    
    coldata = as.data.frame(object@coldata)
    alpha0 = as.numeric(Alpha)
    genes = intersect(genes, rownames(object@rowdata))
    
    if('Integration' %in% names(object@metadata)){
      
      count = object@assay$logcount
      count = lapply(count, FUN = function(y){y[genes,]})
      if(length(genes) > 1){
        count = do.call(cbind, count)
      } else if(length(genes) == 1){
        count = Matrix::Matrix(unlist(count), nrow = 1)
      } else {
        stop("Input valid gene symbol")
      }
      
    } else {
      count = object@assay$logcount[genes, , drop = F]
    }
    
    count = winsorizing(count)
    
    if(!is.null(colFactor)){
      
      if(length(genes) > 1) {
        m0 = data.frame(Factor = coldata[,colFactor], t(as.matrix(count)))
      } else {
        m0 = data.frame(Factor = coldata[,colFactor], Gene = count[1,])
      }
      colnames(m0) = c('Factor', genes)
      
      if(is.null(Colors)){
        Cols0 = Color0(length(levels(m0$Factor)))
      } else {
        Cols0 = Colors
      }

      if(isTRUE(legend)){
        if(isTRUE(dots)){
          plot0 = lapply(c(2L:ncol(m0)), function(z){
            mi = m0[,c(1, z)]
            # mi = mi[mi[,2] > 0,]
            ggplot(mi, aes(Factor, mi[,2])) + 
              geom_violin(aes(fill = Factor), scale = wid, trim = trim, show.legend = TRUE, alpha = alpha0, position = position_dodge()) + 
              geom_jitter(aes(fill = Factor), position = position_jitterdodge(jitter.width = 0.9), size = 0.5, show.legend = FALSE) + 
              scale_fill_manual(values = Cols0) + 
              theme_classic(base_size = 12, base_line_size = 1) + 
              labs(x = '', y = 'Expression', fill = '', title = colnames(mi)[2]) + 
              theme(plot.title = element_text(hjust = 0.5, size = 16), axis.text.x = element_blank())
          })
        } else {
          plot0 = lapply(c(2L:ncol(m0)), function(z){
            mi = m0[,c(1, z)]
            # mi = mi[mi[,2] > 0,]
            ggplot(mi, aes(Factor, mi[,2])) + 
              geom_violin(aes(fill = Factor), scale = wid, trim = trim, show.legend = TRUE, alpha = alpha0, position = position_dodge()) + 
              scale_fill_manual(values = Cols0) + 
              theme_classic(base_size = 12, base_line_size = 1) + 
              labs(x = '', y = 'Expression', fill = '', title = colnames(mi)[2]) + 
              theme(plot.title = element_text(hjust = 0.5, size = 16), axis.text.x = element_blank())
          })
        }
      } else {
        if(isTRUE(dots)){
          plot0 = lapply(c(2L:ncol(m0)), function(z){
            mi = m0[,c(1, z)]
            # mi = mi[mi[,2] > 0,]
            ggplot(mi, aes(Factor, mi[,2])) + 
              geom_violin(aes(fill = Factor), scale = wid, trim = trim, show.legend = FALSE, alpha = alpha0, position = position_dodge()) + 
              geom_jitter(aes(fill = Factor), position = position_jitterdodge(jitter.width = 0.9), size = 0.5, show.legend = FALSE) + 
              scale_fill_manual(values = Cols0) + 
              theme_classic(base_size = 12, base_line_size = 1) + 
              labs(x = '', y = 'Expression', fill = '', title = colnames(mi)[2]) + 
              theme(plot.title = element_text(hjust = 0.5, size = 16), axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5))
          })
        } else {
          plot0 = lapply(c(2L:ncol(m0)), function(z){
            mi = m0[,c(1, z)]
            # mi = mi[mi[,2] > 0,]
            ggplot(mi, aes(Factor, mi[,2])) + 
              geom_violin(aes(fill = Factor), scale = wid, trim = trim, show.legend = FALSE, alpha = alpha0, position = position_dodge()) + 
              scale_fill_manual(values = Cols0) + 
              theme_classic(base_size = 12, base_line_size = 1) + 
              labs(x = '', y = 'Expression', fill = '', title = colnames(mi)[2]) + 
              theme(plot.title = element_text(hjust = 0.5, size = 16), axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5))
          })
        }
      }
      
      return(do.call(grid.arrange, c(plot0, ncol = ceiling(sqrt(length(genes))))))
      
    } else if(is.null(colFactor) & !is.null(object@coldata$Cluster)){
      
      if(length(genes) > 1) {
        m0 = data.frame(Factor = coldata[,'Cluster'], t(as.matrix(count)))
      } else {
        m0 = data.frame(Factor = coldata[,'Cluster'], Gene = count)
      }
      colnames(m0) = c('Factor', genes)
      
      if(is.null(Colors)){
        Cols0 = Color0(length(levels(m0$Factor)))
      } else {
        Cols0 = Colors
      }
      
      if(isTRUE(legend)){
        if(isTRUE(dots)){
          plot0 = lapply(c(2L:ncol(m0)), function(z){
            mi = m0[,c(1, z)]
            # mi = mi[mi[,2] > 0,]
            ggplot(mi, aes(Factor, mi[,2])) + 
              geom_violin(aes(fill = Factor), scale = wid, trim = trim, show.legend = TRUE, alpha = alpha0, position = position_dodge()) + 
              geom_jitter(aes(fill = Factor), position = position_jitterdodge(jitter.width = 0.9), size = 0.5, show.legend = FALSE) + 
              scale_fill_manual(values = Cols0) + 
              theme_classic(base_size = 12, base_line_size = 1) + 
              labs(x = '', y = 'Expression', fill = '', title = colnames(mi)[2]) + 
              theme(plot.title = element_text(hjust = 0.5, size = 16), axis.text.x = element_blank())
          })
        } else {
          plot0 = lapply(c(2L:ncol(m0)), function(z){
            mi = m0[,c(1, z)]
            # mi = mi[mi[,2] > 0,]
            ggplot(mi, aes(Factor, mi[,2])) + 
              geom_violin(aes(fill = Factor), scale = wid, trim = trim, show.legend = TRUE, alpha = alpha0, position = position_dodge()) + 
              scale_fill_manual(values = Cols0) + 
              theme_classic(base_size = 12, base_line_size = 1) + 
              labs(x = '', y = 'Expression', fill = '', title = colnames(mi)[2]) + 
              theme(plot.title = element_text(hjust = 0.5, size = 16), axis.text.x = element_blank())
          })
        }
      } else {
        if(isTRUE(dots)){
          plot0 = lapply(c(2L:ncol(m0)), function(z){
            mi = m0[,c(1, z)]
            # mi = mi[mi[,2] > 0,]
            ggplot(mi, aes(Factor, mi[,2])) + 
              geom_violin(aes(fill = Factor), scale = wid, trim = trim, show.legend = FALSE, alpha = alpha0, position = position_dodge()) + 
              geom_jitter(aes(fill = Factor), position = position_jitterdodge(jitter.width = 0.9), size = 0.5, show.legend = FALSE) + 
              scale_fill_manual(values = Cols0) + 
              theme_classic(base_size = 12, base_line_size = 1) + 
              labs(x = '', y = 'Expression', fill = '', title = colnames(mi)[2]) + 
              theme(plot.title = element_text(hjust = 0.5, size = 16), axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5))
          })
        } else {
          plot0 = lapply(c(2L:ncol(m0)), function(z){
            mi = m0[,c(1, z)]
            # mi = mi[mi[,2] > 0,]
            ggplot(mi, aes(Factor, mi[,2])) + 
              geom_violin(aes(fill = Factor), scale = wid, trim = trim, show.legend = FALSE, alpha = alpha0, position = position_dodge()) + 
              scale_fill_manual(values = Cols0) + 
              theme_classic(base_size = 12, base_line_size = 1) + 
              labs(x = '', y = 'Expression', fill = '', title = colnames(mi)[2]) + 
              theme(plot.title = element_text(hjust = 0.5, size = 16), axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5))
          })
        }
      }
      
      return(do.call(grid.arrange, c(plot0, ncol = ceiling(sqrt(length(genes))))))
      
    } else {
      stop("Please run scCluster first")
    }
    
  }
  
}



####################################################################################
####################################################################################
Color0 <- function(n) {
  qual = brewer.pal.info[rownames(brewer.pal.info) %in% c('Set1', 'Set3', 'Dark2', 'Paired'),]
  Col0 = unique(unlist(mapply(brewer.pal, qual$maxcolors, rownames(qual))))
  Col0 = Col0[order(Col0, decreasing = FALSE)]
  Col0 = Col0[-c(37:40)]
  Col1 = col2rgb(Col0)
  dist.color = as.matrix(dist(t(Col1)))
  diag(dist.color) = 1e10
  while(length(Col0) > n) {
    minCol = apply(dist.color, 1, FUN = min)
    ids = which(minCol == min(minCol))[1]
    dist.color = dist.color[-ids, -ids]
    Col0 = Col0[-ids]
  }
  return(Col0)
}

Color1 <- function(Cols0, g){
  l0 = g[g > 0]
  l1 = (l0 - min(l0))/(max(l0) - min(l0))
  Cols0[g > 0] = colorRampPalette(c("gray80", "#FEE0D2", "#FCBBA1", "#FC9272", "#FB6A4A", "#EF3B2C", "#CB181D", "#A50F15"))(length(l1))
  return(Cols0)
}


