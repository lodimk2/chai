#' Add ARI Function
#'
#' Function to run ARI
#' @param clusters Dataframe column containing clustering assignments 
#' @param ground_truth Dataframe column containing ground_truth assignments
#' @keywords ARI
#' @import mclust
#' @export 


ari_function <- function(clusters, ground_truth) {
    return(mclust::adjustedRandIndex(ground_truth, clusters))
}


#' NMI Function
#'
#' Function to run NMI
#' @param clusters Dataframe column containing clustering assignments 
#' @param ground_truth Dataframe column containing ground_truth assignments
#' @keywords NMI
#' @import aricode
#' @export 
#' @examples
#' nmi_function(snf_clusters@.Data, ground_truth$clust_assign)

nmi_function <- function(clusters, ground_truth) {
    return(aricode::NMI(ground_truth, clusters))
}




#' Create Eval Table
#'
#' Function to run ARI
#' @param clusters Dataframe column containing clustering assignments 
#' @param ground_truth Dataframe column containing ground_truth assignments
#' @keywords NMI
#' @import aricode
#' @import mclust
#' @export 
#' @examples
#' evaluate_table(sce, ground_truth$clust_assign)
evaluation_table <- function(sce, ground_truth) {
  df <- as.data.frame(colData(sce))
  # Initialize lists to store results
  algorithms <- character()
  ari_scores <- numeric()
  nmi_scores <- numeric()
  
  # Iterate over columns
  for (col_name in names(df)) {
    # Check if column name contains "_assign"
    if (grepl("_assign", col_name)) {
      # Calculate ARI
      print(paste0("Calculating evaluation metrics for", col_name))
      ari_score <- ari_function(df[[col_name]], ground_truth)
      # Calculate NMI
      nmi_score <- nmi_function(df[[col_name]], ground_truth)
      
      # Store column name (algorithm)
      algorithms <- c(algorithms, col_name)
      # Store ARI score
      ari_scores <- c(ari_scores, ari_score)
      # Store NMI score
      nmi_scores <- c(nmi_scores, nmi_score)
    }
  }
  
  # Create dataframe with results
  results_df <- data.frame(Algorithms = algorithms,
                           ARI = ari_scores,
                           NMI = nmi_scores)
  results_df$Algorithms <- gsub("_assign", "", results_df$Algorithms)
  
  return(results_df)
}