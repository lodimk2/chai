#' Add ARI Function
#'
#' Function to run ARI
#' @param clusters Dataframe column containing clustering assignments 
#' @param ground_truth Dataframe column containing ground_truth assignments
#' @keywords ARI
#' @import mclust
#' @export 
#' @examples
#' eval_function(snf_clusters@.Data, ground_truth$clust_assign)

eval_function <- function(clusters, ground_truth) {
    return(mclust::adjustedRandIndex(ground_truth, clusters))
}