library(testthat)
library(SingleCellExperiment)
library(chai)

data

test_that("get_clust_assignments runs all clustering algorithms and adds results to SCE object", {
  # Make mock sce data 
  data("baron_mouse_1")

  baron_mouse_1 <- baron_mouse_1[, c(1:100)]

  # Create SingleCellExperiment object
  sce <- SingleCellExperiment(assays = list(counts = as.matrix(baron_mouse_1)))
# Add logcounts 
  sce <- scuttle::logNormCounts(sce)

  # Testing get_clust_assignments with mock data 
  sce <- get_clust_assignments(sce, n_cores = 1, svm_max = 50, max_k = 10)
  
  # Check if clustering results are added to the SCE object's colData
  expect_true("seurat_louvain_assign" %in% colnames(colData(sce)))
  expect_true("raceid_assign" %in% colnames(colData(sce)))
  expect_true("spectrum_assign" %in% colnames(colData(sce)))
  expect_true("scSHC_assign" %in% colnames(colData(sce)))
  expect_true("sc3_assign" %in% colnames(colData(sce)))
  expect_true("choir_assign" %in% colnames(colData(sce)))
  expect_true("seurat_slc_assign" %in% colnames(colData(sce)))
  
})

test_that("CHAI_AvgSim assigns clusters correctly", {
  # Mock  SingleCellExperiment object
  sce <- SingleCellExperiment::SingleCellExperiment(list(counts = matrix(runif(100), ncol = 10)))
  best_k <- 3

  
  result_sce <- CHAI_AvgSim(sce, best_k, eval = FALSE)

  # Check if the result has the expected properties
  expect_true(all(SingleCellExperiment::colData(result_sce)$CHAI_AvgSim_assign %in% 1:best_k))
  expect_equal(length(unique(SingleCellExperiment::colData(result_sce)$CHAI_AvgSim_assign)), best_k)
})


test_that("CHAI_SNF assigns clusters correctly", {
  # Mock SingleCellExperiment object
  sce <- SingleCellExperiment::SingleCellExperiment(list(counts = matrix(runif(100), ncol = 10)))
  best_k <- 3

  result_sce <- CHAI_SNF(sce, best_k, eval = FALSE)

  # Check if the result has the expected properties
  expect_true(all(SingleCellExperiment::colData(result_sce)$CHAI_SNF_assign %in% 1:best_k))
  expect_equal(length(unique(SingleCellExperiment::colData(result_sce)$CHAI_SNF_assign)), best_k)
})
