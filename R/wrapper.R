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
  AvgSim_matrix <- create_average_matrix(sce)
  
  # Determines silhouette scores to find the best k
  
  if (eval == TRUE) {
  AvgSim_best_k <- calc_silhouette_scores(AvgSim_matrix, best_k)
  }
  else {
    AvgSim_best_k <- best_k
  }
  AvgSim_clusters <- spectral_clustering(AvgSim_matrix, AvgSim_best_k)
  
  SingleCellExperiment::colData(sce)$CHAI_AvgSim_assign <- AvgSim_clusters@.Data
  
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
  snf_matrix <- create_average_matrix(sce)
  
  if (eval == TRUE) {
  # Detremines silhouette scores to find the best k
  snf_best_k <- calc_silhouette_scores(snf_matrix, best_k)
  }else{
    snf_best_k <- best_k
  }
  
  snf_clusters <- spectral_clustering(snf_matrix, snf_best_k)
  SingleCellExperiment::colData(sce)$chai_snf_assign <- snf_clusters@.Data
  

  return(sce)
} 

