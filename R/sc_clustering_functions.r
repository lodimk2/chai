#' Run Spectrum
#'
#' Function to run Spectrum clustering algorithm 
#' @param sce SingleCellExperiment object of the dataset
#' @param out_dir Directory where output should be written to. Defaults to the current working directory. Will automatically create a directory called 'alg_clust_assign' to house all algorithm assignments
#' @keywords Spectrum
#' @import Spectrum
#' @export 
#' @examples
#' spectrum_assign_func(sce, out_dir = getwd())
spectrum_assign_func <- function(sce, out_dir = getwd()) {

    if (!dir.exists(paste0(out_dir, "/alg_clust_assign"))){
    dir.create(paste0(out_dir, "/alg_clust_assign"))
    }else{
    print("alg_clust_assign dir already exists")
    }

    sim_df <- as.data.frame(as.matrix(counts(sce)))
    colnames(sim_df) <- 1:length(colnames(sim_df))
    rownames(sim_df) <- 1:length(rownames(sim_df))

    test3 <- Spectrum::Spectrum(sim_df,showres=FALSE,runrange=TRUE,krangemax=10)
    spectrum_assignments <- cbind(1:length(colnames(sim_df)),test3[[2]]$assignments)
    colnames(spectrum_assignments) <- c("cells", "clust_assign")


    write.csv(spectrum_assignments, paste0(out_dir, "/alg_clust_assign/spectrum_assign.csv"))

}

#' Run Seurat
#'
#' Function to run Seurat clustering algorithm. Runs Seurat SLC and Seurat Louvain clustering 
#' @param sce SingleCellExperiment object of the dataset
#' @param out_dir Directory where output should be written to. Defaults to the current working directory. Will automatically create a directory called 'alg_clust_assign' to house all algorithm assignments
#' @keywords Seurat
#' @import Seurat
#' @export
#' @examples
#' seurat_assign_func(sce, out_dir = getwd())
seurat_assign_func <- function(sce, out_dir = getwd()) {
    
    if (!dir.exists(paste0(out_dir, "/alg_clust_assign"))){
    dir.create(paste0(out_dir, "/alg_clust_assign"))
    }else{
    #print("alg_clust_assign dir already exists")
    }
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

    write.csv(seurat_louvain_assign, paste0(out_dir, "/alg_clust_assign/seurat_louvain_assign.csv"))
    write.csv(seurat_slc_assign, paste0(out_dir, "/alg_clust_assign/seurat_slc_assign.csv"))

}

#' Run scSHC
#'
#' Function to run scSHC clustering algorithm.
#' @param sce SingleCellExperiment object of the dataset
#' @param out_dir Directory where output should be written to. Defaults to the current working directory. Will automatically create a directory called 'alg_clust_assign' to house all algorithm assignments
#' @keywords scSHC
#' @export 
#' @examples
#' scSHC_assign_func(sce, out_dir = getwd())
scSHC_assign_func <- function(sce, out_dir = getwd()) {

    if (!dir.exists(paste0(out_dir, "/alg_clust_assign"))){
    dir.create(paste0(out_dir, "/alg_clust_assign"))
    }else{
    #print("alg_clust_assign dir already exists")
    }

    clusters <- scSHC::scSHC(counts(sce))
    #View(clusters[[1]])

    scSHC_clust_assign <- as.data.frame(clusters[[1]])
    scSHC_clust_assign$cells <- rownames(scSHC_clust_assign)
    colnames(scSHC_clust_assign)[1] <- "clust_assign"

    scSHC_clust_assign <- scSHC_clust_assign[, c(2, 1)]

    write.csv(scSHC_clust_assign,paste0(out_dir, "/alg_clust_assign/scSHC_clust_assign.csv"))

}

#' Run RaceID
#'
#' Function to run RaceID clustering algorithm.
#' @param sce SingleCellExperiment object of the dataset
#' @param out_dir Directory where output should be written to. Defaults to the current working directory. Will automatically create a directory called 'alg_clust_assign' to house all algorithm assignments
#' @keywords RaceID
#' @export
#' @examples
#' raceid_assign_func(sce, out_dir = getwd())
raceid_assign_func <- function(sce, out_dir = getwd()) {

    if (!dir.exists(paste0(out_dir, "/alg_clust_assign"))){
    dir.create(paste0(out_dir, "/alg_clust_assign"))
    } else{
    #print("alg_clust_assign dir already exists")
    }

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
    write.csv(raceid_ground_truth, paste0(out_dir, "/alg_clust_assign/raceid_assign.csv"))


}

#' Run SC3
#'
#' Function to run SC3 clustering algorithm.
#' @param sce SingleCellExperiment object of the dataset
#' @param out_dir Directory where output should be written to. Defaults to the current working directory. Will automatically create a directory called 'alg_clust_assign' to house all algorithm assignments
#' @keywords RaceID
#' @export
#' @examples
#' sce3_assign_func(sce, out_dir = getwd())
sc3_assign_func <- function(sce, out_dir = getwd()) {

    if (!dir.exists(paste0(out_dir, "/alg_clust_assign"))){
    dir.create(paste0(out_dir, "/alg_clust_assign"))
    }else{
    #print("alg_clust_assign dir already exists")
    }

    # STARTING NEW SC3 FUNCTION
    sce_sc3 <- SingleCellExperiment::SingleCellExperiment(
    assays = list(
        counts = as.matrix(counts(sce)),
        logcounts = log2(as.matrix(counts(sce)) + 1)
    ),
    colData = colData(sce),
    )
    
    rowData(sce_sc3)$feature_symbol <- rownames(sce_sc3)
    sce_sc3 <- sce_sc3[!duplicated(rowData(sce_sc3)$feature_symbol), ]
    sce_sc3 <- SC3::sc3_estimate_k(sce_sc3)
   
    sce_sc3 <- SC3::sc3(sce_sc3, ks = 2:as.integer(metadata(sce_sc3)$sc3$k_estimation), n_cores = 4)
    col_data <- colData(sce_sc3)
    last_col <- col_data[, tail(seq_len(ncol(col_data)), 1)]
    
    clust_assign_df <- cbind(rownames(col_data), last_col)
    colnames(clust_assign_df) <- c("cells", "clust_assign")
    write.csv(clust_assign_df,paste0(out_dir, "/alg_clust_assign/sc3_clust_assign.csv"))
    
}

#' Run Choir
#'
#' Function to run choir clustering algorithm 
#' @param sce SingleCellExperiment object of the dataset
#' @param out_dir Directory where output should be written to. Defaults to the current working directory. Will automatically create a directory called 'alg_clust_assign' to house all algorithm assignments
#' @keywords choir
#' @export 
#' @examples
#' choir_assign_func(sce, out_dir = getwd())
choir_assign_func <- function(sce, out_dir = getwd()) {
  
  if (!dir.exists(paste0(out_dir, "/alg_clust_assign"))){
    dir.create(paste0(out_dir, "/alg_clust_assign"))
  }else{
    #print("alg_clust_assign dir already exists")
  }
  

  seurat_object <- Seurat::CreateSeuratObject(counts(sce))
  seurat_object <- Seurat::NormalizeData(seurat_object)
  # TEST LINE TO SEE IF THIS WORKS
  seurat_object <- seurat_object
  seurat_object <- CHOIR::CHOIR(seurat_object, n_cores = 2)
  seurat_object <- CHOIR::buildTree(seurat_object, n_cores = 2)
  seurat_object <- CHOIR::pruneTree(seurat_object, n_cores = 2)

  choir_assignments <- cbind(colnames(sce), seurat_object$CHOIR_clusters_0.05)

  colnames(choir_assignments) <- c("cells", "clust_assign")
  write.csv(choir_assignments,paste0(out_dir, "/alg_clust_assign/choir_clust_assign.csv"))

}

#' Run All Clustering Assignment Algorithms currently in the pipeline
#'
#' Function to run all clustering algorithms (Seurat, scSHC, Spectrum, RaceID, SC3).
#' @param sce SingleCellExperiment object of the dataset
#' @param out_dir Directory where output should be written to. Defaults to the current working directory. Will automatically create a directory called 'alg_clust_assign' to house all algorithm assignments
#' @keywords All Algorithm Wrapper Function
#' @export
#' @examples
#' get_clust_assignments(sce, out_dir = getwd())
get_clust_assignments <- function(sce, out_dir = getwd(), windows = FALSE) {

    if (!dir.exists(paste0(out_dir, "/alg_clust_assign"))){
    dir.create(paste0(out_dir, "/alg_clust_assign"))
    } else{
    print("alg_clust_assign dir already exists")
    }

    print("Running Seurat")
    invisible(capture.output(seurat_assign_func(sce, out_dir)))
    print("Running Spectrum")
    invisible(capture.output(spectrum_assign_func(sce, out_dir)))
    print("Running RaceID")
    invisible(capture.output(raceid_assign_func(sce, out_dir)))
    print("Running SC3")
    invisible(capture.output(sc3_assign_func(sce, out_dir)))
    print("Running CHOIR")
    invisible(capture.output(choir_assign_func(sce, out_dir)))
    if (windows == TRUE) {
        "Skipping scSHC since Windows machine specified"
    }
    else {
    print("Running scSHC")
    invisible(capture.output(scSHC_assign_func(sce, out_dir)))
    }

    print("All clustering assignments are completed.")
    print(paste0("Output written to ", out_dir, "/alg_clust_assign"))
}