#' Run Spectrum
#'
#' Function to run Spectrum clustering algorithm
#' @param sce SingleCellExperiment object of the dataset
#' @keywords Spectrum
#' @import Spectrum
#' @export
#' @import SingleCellExperiment
#' @examples
#' spectrum_assign_func(sce)
spectrum_assign_func <- function(sce) {

    sim_df <- as.data.frame(as.matrix(counts(sce)))
    colnames(sim_df) <- 1:length(colnames(sim_df))
    rownames(sim_df) <- 1:length(rownames(sim_df))

    test3 <- Spectrum::Spectrum(sim_df,showres=FALSE,runrange=TRUE,krangemax=10)
    spectrum_assignments <- cbind(1:length(colnames(sim_df)),test3[[2]]$assignments)
    spectrum_assignments <- as.data.frame(spectrum_assignments)
    colnames(spectrum_assignments) <- c("cells", "clust_assign")
    SingleCellExperiment::colData(sce)$spectrum_assign <- spectrum_assignments$clust_assign
    return(sce)
}

#' Run Seurat
#'
#' Function to run Seurat clustering algorithm. Runs Seurat SLC and Seurat Louvain clustering 
#' @param sce SingleCellExperiment object of the dataset
#' @keywords Seurat
#' @import Seurat
#' @export
#' @import SingleCellExperiment
#' @examples
#' seurat_assign_func(sce)
seurat_assign_func <- function(sce) {
    
    sim_groups_convert <- Seurat::CreateSeuratObject(counts = counts(sce))
    #sim_groups_convert <- Seurat::as.Seurat(sce, counts = counts(sce))
    sim_groups_convert <- Seurat::NormalizeData(sim_groups_convert, normalization.method = "LogNormalize", scale.factor = 10000)
    sim_groups_convert <- Seurat::FindVariableFeatures(sim_groups_convert, selection.method = "vst", nfeatures = 2000)
    all.genes <- rownames(sim_groups_convert)
    sim_groups_convert <- Seurat::ScaleData(sim_groups_convert, features = all.genes)
    sim_groups_convert <- Seurat::RunPCA(sim_groups_convert, features = VariableFeatures(object = sim_groups_convert))
    sim_groups_convert <- Seurat::FindNeighbors(sim_groups_convert, dims = 1:10)
    sim_groups_convert <- Seurat::FindClusters(sim_groups_convert, resolution = 0.5)


    seurat_louvain_assign <- as.data.frame(Seurat::Idents(sim_groups_convert))
    seurat_louvain_assign$cells <- rownames(seurat_louvain_assign)
    colnames(seurat_louvain_assign) <- c("clust_assign", "cells")

    seurat_louvain_assign <- seurat_louvain_assign[, c(2,1)]

    sim_slc <- Seurat::FindNeighbors(sim_groups_convert, dims = 1:10)
    sim_slc <- Seurat::FindClusters(sim_groups_convert, resolution = 0.5, algorithm = 3)

    seurat_slc_assign <- as.data.frame(Seurat::Idents(sim_slc))
    seurat_slc_assign$cells <- rownames(seurat_slc_assign)
    colnames(seurat_slc_assign) <- c("clust_assign", "cells")
    seurat_slc_assign <- seurat_slc_assign[, c(2,1)]

    SingleCellExperiment::colData(sce)$seurat_louvain_assign <- seurat_louvain_assign$clust_assign
    SingleCellExperiment::colData(sce)$seurat_slc_assign <- seurat_slc_assign$clust_assign
    return(sce)
}

#' Run scSHC
#'
#' Function to run scSHC clustering algorithm.
#' @param sce SingleCellExperiment object of the dataset
#' @param n_cores CPU cores for parallel library allocated for running functions. 
#' @keywords scSHC
#' @export
#' @import SingleCellExperiment
#' @examples
#' scSHC_assign_func(sce)
scSHC_assign_func <- function(sce, n_cores=1) {
    
    clusters <- scSHC::scSHC(SingleCellExperiment::counts(sce), cores=n_cores)
    scSHC_clust_assign <- as.data.frame(clusters[[1]])
    scSHC_clust_assign$cells <- rownames(scSHC_clust_assign)
    colnames(scSHC_clust_assign)[1] <- "clust_assign"

    scSHC_clust_assign <- scSHC_clust_assign[, c(2, 1)]

    SingleCellExperiment::colData(sce)$scSHC_assign <- scSHC_clust_assign$clust_assign
    return(sce)
}

#' Run RaceID
#'
#' Function to run RaceID clustering algorithm.
#' @param sce SingleCellExperiment object of the dataset
#' @keywords RaceID
#' @export
#' @import SingleCellExperiment
#' @examples
#' raceid_assign_func(sce)
raceid_assign_func <- function(sce) {

    sc <- RaceID::SCseq(counts(sce))
    sc <- RaceID::filterdata(sc,mintotal=10)
    fdata <- RaceID::getfdata(sc)
    sc <- RaceID::compdist(sc,metric="pearson")
    sc <- RaceID::clustexp(sc)
    raceid_ground_truth <- as.data.frame(sc@cluster$kpart)
    print("RACE ID RESULTS HEAD")
    raceid_ground_truth$cells <- rownames(raceid_ground_truth)
    print(head(raceid_ground_truth))
    #colnames(raceid_ground_truth) <- c("cells", "clust_assign")
    raceid_ground_truth <- raceid_ground_truth[, c(2, 1)]
    colnames(raceid_ground_truth) <- c("cells", "clust_assign")
    SingleCellExperiment::colData(sce)$raceid_assign <- raceid_ground_truth$clust_assign
    return(sce)
}

#' Run SC3
#'
#' Function to run SC3 clustering algorithm.
#' @param sce SingleCellExperiment object of the dataset
#' @param n_cores Computer cores for parallel library allocated for running functions. 
#' @param svm_max Max number of cells to use for SVM. Default is 1000 
#' @param max_k max_k for sc3 to run to determine optimal clustering.
#' @keywords SC3
#' @export
#' @import SingleCellExperiment
#' @examples
#' sce3_assign_func(sce)
sc3_assign_func <- function(sce, n_cores=1, svm_max = 1000, max_k = 15) {
    sce_sc3 <- sce
    SingleCellExperiment::rowData(sce_sc3)$feature_symbol <- rownames(sce_sc3)
    sce_sc3 <- sce_sc3[!duplicated(SingleCellExperiment::rowData(sce_sc3)$feature_symbol), ]
    sce_sc3 <- SC3::sc3(sce_sc3, ks = 2:as.integer(max_k), n_cores = n_cores, svm_max=svm_max)
    #print("Ran sc3_prepare")
    sce_sc3 <- SC3::sc3_run_svm(sce_sc3, ks = 2:as.integer(max_k))
    col_data <- SingleCellExperiment::colData(sce_sc3)
    last_col <- col_data[, tail(seq_len(ncol(col_data)), 1)]
    
    clust_assign_df <- cbind(rownames(col_data), last_col)
    clust_assign_df <- as.data.frame(clust_assign_df)
    colnames(clust_assign_df) <- c("cells", "clust_assign")
    SingleCellExperiment::colData(sce)$sc3_assign <- clust_assign_df$clust_assign 
    return(sce)
}

#' Run Choir
#'
#' Function to run choir clustering algorithm 
#' @param sce SingleCellExperiment object of the dataset
#' @keywords choir
#' @export
#' @import SingleCellExperiment
#' @examples
#' choir_assign_func(sce)
choir_assign_func <- function(sce, n_cores=1) {
  
  seurat_object <- CHOIR::CHOIR(sce, n_cores = n_cores)
  choir_assignments <- cbind(colnames(sce), seurat_object$CHOIR_clusters_0.05)
  choir_assignments <- as.data.frame(choir_assignments)
  colnames(choir_assignments) <- c("cells", "clust_assign")
  SingleCellExperiment::colData(sce)$choir_assign <- choir_assignments$clust_assign
  return(sce)
}

#' Run All Clustering Assignment Algorithms currently in the pipeline
#'
#' Function to run all clustering algorithms (Seurat, scSHC, Spectrum, RaceID, SC3).
#' @param sce SingleCellExperiment object of the dataset
#' @param n_cores Computer cores for parallel library allocated for running functions. 
#' @param svm_max define the maximum number of cells below which SVM is not run.
#' @param max_k maxk range for sc3
#' @keywords All Algorithm Wrapper Function
#' @export
#' @examples
#' get_clust_assignments(sce)
get_clust_assignments <- function(sce, n_cores=1, svm_max=500, max_k=15) {
    print("Running Seurat")
    sce <- seurat_assign_func(sce)
    print("Running RaceID")
    sce <- raceid_assign_func(sce)
    print("Running Spectrum")
    sce <- spectrum_assign_func(sce)
    print("Running SC3")
    sce <- sc3_assign_func(sce, n_cores=n_cores, svm_max=svm_max, max_k=max_k)
    print("Running scSHC")
    sce <- (scSHC_assign_func(sce, n_cores=n_cores))
    print("Running CHOIR")
    sce <- choir_assign_func(sce, n_cores=n_cores)
    print("All clustering assignments are completed.")
    return(sce)
}