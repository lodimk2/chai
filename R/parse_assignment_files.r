#' Create List of Similarity Matrices based on Clustering Assignments
#'
#' Create Similarity Matrices for each clustering assignment and add it to a list. Parses the "alg_clust_assign" directory. List can be used for both CSPA and SNF
#' @param home_dir Directory where subdirectories are stored, such as ground truth, alg_clust_assign, data. 
#' @keywords Similarity Matrix List Creation
#' @export 
#' @examples
#' create_matrix_list(getwd())
create_matrix_list <- function(home_dir = getwd()) {
    sim_matrix_list <- list()
    setwd(paste0(home_dir, "/alg_clust_assign"))
    for (file in list.files(getwd())) {
            assignment_df <- read.csv(file)
            sim_mat <- create_matrix(assignment_df)
            sim_matrix_list[[file]] <- sim_mat
    }
    setwd(home_dir)
    print("Similarity Matrix List Created")
    return(sim_matrix_list)
}
#' Creates a Binary Similarity Matrix for Clustering Assignments. Internal Function. 
#'
#' To be used internally for the create_matrix_list function 
#' @param home_dir SingleCellExperiment object of the dataset
#' @keywords Binary Similarity Matrix Creation
#' @examples
#' create_matrix_list(getwd())
create_matrix <- function(assignment_df) {

    clusters <- unique(assignment_df$clust_assign)
    # Create an empty matrix
    similarity_matrix <- matrix(0, nrow = nrow(assignment_df), ncol = nrow(assignment_df))

    # Assign similarity scores
    for (i in 1:length(clusters)) {
    idx <- assignment_df$clust_assign == clusters[i]
    similarity_matrix[idx, idx] <- 1
    }

   diag(similarity_matrix) <- 1

   return(similarity_matrix)
}