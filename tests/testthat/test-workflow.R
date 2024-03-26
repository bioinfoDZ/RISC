context("test-workflow")

# Load raw data
# load("../testdata/testdata.rda")

mat0 = as.matrix(raw.mat[[1]])
coldata0 = as.data.frame(raw.mat[[2]])
coldata1 <- coldata0[coldata0$Batch0 == 'Batch1',]
coldata2 <- coldata0[coldata0$Batch0 == 'Batch4',]
mat1 <- mat0[,rownames(coldata1)]
mat2 <- mat0[,rownames(coldata2)]


#####################################################################################
context("Creat RISC object")

sce1 <- readsc(count = mat1, cell = coldata1, gene = data.frame(Symbol = rownames(mat1), row.names = rownames(mat1)), is.filter = FALSE)
sce2 <- readsc(count = mat2, cell = coldata2, gene = data.frame(Symbol = rownames(mat2), row.names = rownames(mat2)), is.filter = FALSE)

test_that("Whether objects are scdataset objects", {
  expect_is(sce1, 'RISCdata')
  expect_is(sce2, 'RISCdata')
})

test_that("col.names of matrix in objects equal to row.names of coldata in objects", {
  expect_equal(colnames(sce1@assay$count), rownames(sce1@coldata))
})

test_that("row.names of matrix in objects equal to row.names of rowdata in objects", {
  expect_equal(rownames(sce2@assay$count), rownames(sce2@rowdata))
})


#####################################################################################
context("Preprocess data")

sce1 = scFilter(sce1, min.UMI = 0, max.UMI = Inf, min.gene = 0, min.cell = 0, is.filter = FALSE)
sce2 = scFilter(sce2, min.UMI = 0, max.UMI = Inf, min.gene = 0, min.cell = 0, is.filter = FALSE)

test_that("Select the correct cells", {
  expect_equal(min(sce1@coldata$nGene), 464)
  expect_equal(max(sce2@coldata$nGene), 522)
})

test_that("Select the correct genes", {
  expect_equal(sce1@rowdata$Symbol[1], "Gene1")
  expect_equal(max(sce2@rowdata$nCell), 25)
  expect_false(is.null(sce1@metadata$filter))
})

sce1 = scNormalize(sce1)
sce2 = scNormalize(sce2)

test_that("Normalize raw counts/UMIs", {
  expect_equal(min(sce1@assay$logcount), 0)
  expect_gt(max(sce2@assay$logcount), 6.4)
  expect_equal(length(sce1@metadata$normalise), 1)
})


#####################################################################################
context("Identify highly variable genes")

sce1 = scDisperse(sce1)
sce2 = scDisperse(sce2)

test_that("highly variable genes", {
  expect_false(is.null(sce1@vargene))
  expect_false(is.null(sce2@metadata$dispersion.var))
})


#####################################################################################
context("Dimension reduction")

sce1 = scPCA(sce1, npc = 5)
sce2 = scPCA(sce2, npc = 5)

test_that("PCA", {
  expect_false(is.null(sce1@DimReduction$cell.pca))
  expect_gt(max(sce2@DimReduction$var.pca), 0.30)
})

sce1 = scTSNE(sce1, npc = 5, perplexity = 5)
sce2 = scTSNE(sce2, npc = 5, perplexity = 5)

test_that("tSNE", {
  expect_false(is.null(sce1@DimReduction$cell.tsne))
})

sce1 = scUMAP(sce1, npc = 5)
sce2 = scUMAP(sce2, npc = 5)

test_that("UMAP", {
  expect_false(is.null(sce1@DimReduction$cell.umap))
})


#####################################################################################
context("Data integration")

var0 = intersect(rownames(sce1@assay$logcount), rownames(sce2@assay$logcount))
idat = list(sce1, sce2)
idat = scMultiIntegrate(idat, eigens = 4, npc = 5, var.gene = var0, add.Id = c('Set1', 'Set2'))

test_that("Integrate datasets", {
  expect_false(is.null(idat@DimReduction$cell.pls))
  expect_equal(length(levels(idat@coldata$Set)), 2)
})


#####################################################################################
context("Cell clustering")

idat = scUMAP(idat, npc = 5, use = 'PLS')
idat = scCluster(idat, slot = 'cell.umap', k = 4, method = 'density', dc = 0.3253)

test_that("Integrate datasets", {
  expect_false(is.null(idat@DimReduction$cell.umap))
  expect_equal(length(levels(idat@coldata$Cluster)), 4)
})


