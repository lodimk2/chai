#' Creates a binary similarity matrix from spatial transcriptomics coordinates based on GraphST's method
#' 
#' @param spat_data  NAME: the cell  - X: x coordinates  - Y: y coordinates - Cluster: cell type annotation. Ensure that "X" coordinates and "Y" coordinates columns are named with "X" and "Y". 
#' @return adj
#' @import FNN
#' @export
ST_Sim <- function(spat_data) {
    position <- spat_data[-1, c("X", "Y")] #exclude top non-data row, and X, Y columns for position info
    position$X <- as.numeric(position$X)
    position$Y <- as.numeric(position$Y) #convert position data to be numeric

    n_spot <- nrow(position) #number of cells
    n_neighbors <- 3 #chosen value of k
    knn_res <- FNN::get.knn(position, n_neighbors)
    nbrs <- knn_res$nn.index #n x k, ith row contains indices of neighbors to ith cell sorted in order of closeness
    
    interaction <- matrix(0, n_spot, n_spot) #initialize interaction/neighbor matrix
    for (i in 1:n_spot){
        interaction[i, nbrs[i, ]] <- 1 #cell i and cell j neighbors => interaction[i][j] = 1
    } 

    # Transform adj to symmetrical adj
    adj <- interaction
    adj <- adj + t(adj)
    adj <- ifelse(adj > 1, 1, adj)

    return(adj)
}

#' Runs CHAI with ST integration with SNF - Two Level
#' 
#' @param st_matrix binary ST matrix 
#' @param similarity_matrix_list similarity matrix list created from clustering assignment list 
#' @return adj
#' @import SNFtool
#' @export
CHAI_ST <- function(st_matrix, sce) {
    
    similarity_matrix_list <- append(similarity_matrix_list, st_matrix)

    snf_matrix <- create_snf_matrix(similarity_matrix_list)

    # Determine best K through silhouette score. 
    snf_best_k <- calc_silhouette_scores(snf_matrix, 6)

    # To run on snf_matrix, replace AvgSim_matrix with SNF_matrix in this function, as well as either the true number of clusters in the dataset or evaluate for best K. 
    snf_clusters <- spectral_clustering(snf_matrix, snf_best_k)

    colData(sce)$CHAI_ST_assign <- snf_clusters@.Data

    return(sce)

}

