#' Determine best k for spectral clustering using silhouette score
#'
#' On either CSPA or SNF matrix, determine what the best silhouette score is and choose it as best K. 
#' @param matrix Either CSPA or SNF matrix 
#' @param max_k Run silhouette score evaluation up to the max_k. Should be upper bound of # of clusters user expects in the dataset
#' @keywords Silhouette Score Determination
#' @import SNFtool
#' @import cluster
#' @import kernlab
#' @export 
#' @examples
#' best_k <- calc_silhouette_score(matrix, 15)
calc_silhouette_scores <- function(matrix, max_k=15) {
    best_score <- -1000000
    best_k <- 0
    dissimilarity_matrix <- 1 - matrix
    for (i in 2:max_k) {
        clusters <- SNFtool::spectralClustering(kernlab::as.kernelMatrix(matrix), K = i)
        sil_score <- cluster::silhouette(clusters@.Data, dissimilarity_matrix)
        mean_sil_score <- mean(sil_score[, "sil_width"])
        #print(mean_sil_score)
        if (mean_sil_score > best_score) {
            best_score <- mean_sil_score
            best_k <- i
        }
    }
    print(paste0("Best K based on Silhouette score is", best_k))
    return(best_k)
}

#' Run Spectral Clustering
#'
#' On either CSPA or SNF matrix, run spectral clustering with user defined k. 
#' @param matrix Either CSPA or SNF matrix 
#' @param best_k Number of clusters that should be run for spectral clustering 
#' @keywords Run Spectral Clustering
#' @import SNFtool
#' @export 
#' @examples
#' best_k <- calc_silhouette_score(matrix, 15)

spectral_clustering <- function(matrix, best_k) {
    clusters <- SNFtool::spectralClustering(as.kernelMatrix(matrix), K = best_k)
}
