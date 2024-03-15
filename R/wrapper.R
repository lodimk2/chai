#'Creates a similarity matrix and does clustering for avg_sim
#' To be used internally for the create_matrix_list function 
#' 
#' @param sce SingleCellExperiment object containing all of the clustering assignments in colData
#' @param best_k Either is the max_k to evaluate k until, or the true_k in the dataset, depending on if EVAL is true
#' @param eval if eval is true, evaluate the best_k, else use best_k to be true_k
#' @return sce
#' @export
#' @import SingleCellExperiment
CHAI_AvgSim <- function(sce,best_k,eval = TRUE) {

  
  # Generate list of binary similarity matrices 
  similarity_matrix_list <- create_matrix_list(sce)
  AvgSim_matrix <- create_average_matrix(similarity_matrix_list)

  # Determines silhouette scores to find the best k
  if (eval == TRUE) {
  print("Evaluating best_k using Silhouette Score")
  AvgSim_best_k <- calc_silhouette_scores(AvgSim_matrix, best_k)
  }
  else {
    print("Using user inputted best_k as best_k")
    AvgSim_best_k <- best_k
  }
  print("Running spectral clustering")
  AvgSim_clusters <- spectral_clustering(AvgSim_matrix, AvgSim_best_k)
  
  colData(sce)$CHAI_AvgSim_assign <- AvgSim_clusters@.Data
  print("CHAI_AvgSim assignment added to colData(sce)")
  return(sce)
}

#'Creates a similarity matrix and does clustering for snf_matrix
#' To be used internally for the create_matrix_list function 
#' 
#' @param sce SingleCellExperiment object containing all of the clustering assignments in colData
#' @param best_k Either is the max_k to evaluate k until, or the true_k in the dataset, depending on if EVAL is true
#' @param eval if eval is true, evaluate the best_k, else use best_k to be true_k
#' @return sce 
#' @export
#' @import SingleCellExperiment
CHAI_SNF <- function(sce,best_k, eval = TRUE) {

  # Generate list of binary similarity matrices 
  similarity_matrix_list <- create_matrix_list(sce)
  print("Calculating SNF matrix. This step may take a while depending on dataset size")
  snf_matrix <- create_snf_matrix(similarity_matrix_list)
  
  if (eval == TRUE) {
  # Detremines silhouette scores to find the best k
  print("Evaluating best_k using Silhouette Score")
  snf_best_k <- calc_silhouette_scores(snf_matrix, best_k)
  
  }else{
    print("Using user inputted best_k as best_k")
    snf_best_k <- best_k
  }
  
  snf_clusters <- spectral_clustering(snf_matrix, snf_best_k)
  colData(sce)$chai_snf_assign <- snf_clusters@.Data
  print("CHAI_SNF assignment added to colData(sce)")
  return(sce)
} 

