#' Compute Gram matrix for a list of vectors using a specified kernel function.
#'
#' @param vec_list A list of vectors for which the Gram matrix needs to be computed.
#' @param kernel_type The type of kernel function to be used (default: "rbf").
#' @param kernel_params A list of parameters for the kernel function (default values provided).
#'
#' @return A Gram matrix computed using the specified kernel function.
#' @export
#' @examples
#' vectors_list <- list(c(1, 2, 3), c(4, 5, 6), c(7, 8, 9))
#' result_gram_matrix <- gram_matrix(vectors_list, kernel_type = "linear")
#' print(result_gram_matrix)


gram_matrix <- function(vec_list,
                        kernel_type = "rbf",
                        kernel_params = list(length_scale = 1, degree = 2,
                                             free_add = 0, free_mult = 1,
                                             nu_matern = 1, centering_param)) {

  # Get all combinations of vector pairs
  combinations <- combn(1:length(vec_list), 2, simplify = TRUE)

  # Apply the kernel function to all combinations and store in an upper triangular vector
  upper_tri_vec <- apply(combinations, 2, function(pair) {
    kernel(vec_list[[pair[1]]], vec_list[[pair[2]]], type = kernel_type,
           length_scale = kernel_params$length_scale, degree = kernel_params$degree,
           free_add = kernel_params$free_add, free_mult = kernel_params$free_mult,
           nu_matern = kernel_params$nu_matern, centering_param = kernel_params$centering_param)
  })

  # Convert the upper triangular vector to a symmetric matrix
  result_gram_matrix <- upper_tri_to_mat(upper_tri_vec)

  # Compute diagonal elements using the kernel function with each vector against itself
  diag_comb <- rbind(1:length(vec_list), 1:length(vec_list))
  diag_vec <- apply(diag_comb, 2, function(pair) {
    kernel(vec_list[[pair[1]]], vec_list[[pair[2]]], type = kernel_type,
           length_scale = kernel_params$length_scale, degree = kernel_params$degree,
           free_add = kernel_params$free_add, free_mult = kernel_params$free_mult,
           nu_matern = kernel_params$nu_matern, centering_param = kernel_params$centering_param)
  })

  # Add the diagonal elements to the result_gram_matrix
  result_gram_matrix <- result_gram_matrix + diag(diag_vec)

  return(result_gram_matrix)
}
