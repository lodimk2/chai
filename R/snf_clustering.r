#' SNF Matrix
#'
#' Function to create SNF matrix
#' @param similarity_matrix_list list of binary similarity matrices
#' @keywords SNF Matrix
#' @import SNFtool
#' @export 
#' @examples \dontrun{
#' snf_matrix <- create_snf_matrix(similarity_matrix_list)
#' }
create_snf_matrix <- function(similarity_matrix_list) {
    SNF_mat <- SNF(similarity_matrix_list, 20, 20)
    return(SNF_mat)
}