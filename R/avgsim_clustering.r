#' Average Similarity Matrix
#'
#' Function to create average similarity matrix for CHAI_AvgSim
#' @param similarity_matrix_list list of binary similarity matrices
#' @keywords Average Similarity Matrix
#' @export 
#' @examples 
#' \dontrun{
#' create_average_matrix(similarity_matrix_list)
#' }
create_average_matrix <- function(similarity_matrix_list) {
    matrix_dim <- dim(similarity_matrix_list[[1]])
    #print(paste0("Creating empty matrix of dimensions", matrix_dim[[1]], " and ", matrix_dim[[2]]))
    average_matrix <- matrix(NA, nrow = matrix_dim[[1]], ncol = matrix_dim[[2]])

    average_matrix <- Reduce("+", similarity_matrix_list) / length(similarity_matrix_list)
    print("Average similarity matrix created")
    return(average_matrix)
}