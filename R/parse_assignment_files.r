#' Create List of Similarity Matrices based on Clustering Assignments
#'
#' Create Similarity Matrices for each clustering assignment and add it to a list. Parses the "alg_clust_assign" directory. List can be used for both CSPA and SNF
#' @param sce Directory where subdirectories are stored, such as ground truth, alg_clust_assign, data. 
#' @keywords Similarity Matrix List Creation
#' @export 
#' @examples
#' create_matrix_list(getwd())
create_matrix_list <- function(sce) {
    dataframe <- as.data.frame(colData(sce))
    sim_matrix_list <- list()
    
    for (col_name in names(dataframe)) {
        if (grepl("_assign", col_name)) {
            sim_mat <- create_matrix(dataframe[[col_name]])
            sim_matrix_list[[col_name]] <- sim_mat
        }
    }
    
    print("Similarity Matrix List Created")
    return(sim_matrix_list)
}
#' Creates a Binary Similarity Matrix for Clustering Assignments. Internal Function. 
#'
#' To be used internally for the create_matrix_list function 
#' @param assignment_col column containing clustering assignments
#' @keywords Binary Similarity Matrix Creation
#' @examples
#' create_matrix_list(getwd())
create_matrix <- function(assignment_col) {
    clusters <- unique(assignment_col)
    n <- length(assignment_col)
    # Create an empty matrix
    similarity_matrix <- matrix(0, nrow = n, ncol = n)

    # Assign similarity scores
    for (i in 1:length(clusters)) {
        idx <- assignment_col == clusters[i]
        similarity_matrix[idx, idx] <- 1
    }

    diag(similarity_matrix) <- 1

    return(similarity_matrix)
}


