library(testthat)
library(SingleCellExperiment)
library(chai)

test_that("get_clust_assignments runs all clustering algorithms and adds results to SCE object", {
  # Make mock sce data 
  n_genes = 100
  n_cells = 50
  # Create a matrix of counts
  counts <- matrix(rpois(n_genes * n_cells, lambda = 1), ncol = n_cells)
  # Create a SingleCellExperiment object with these counts
  sce <- SingleCellExperiment(list(counts = counts))
  
  rowData(sce)$gene_id <- paste0("Gene", seq_len(n_genes))
  colData(sce)$cell_id <- paste0("Cell", seq_len(n_cells))
  
  # Testing get_clust_assignments with mock data 
  sce <- get_clust_assignments(sce, n_cores = 1, svm_max = 200, max_k = 15)
  
  # Check if clustering results are added to the SCE object's colData
  expect_true("seurat_louvain_assign" %in% colnames(colData(sce)))
  expect_true("raceid_assign" %in% colnames(colData(sce)))
  expect_true("spectrum_assign" %in% colnames(colData(sce)))
  expect_true("scSHC_assign" %in% colnames(colData(sce)))
  expect_true("sc3_assign" %in% colnames(colData(sce)))
  expect_true("choir_assign" %in% colnames(colData(sce)))
})
