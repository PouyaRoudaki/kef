#' Centered Kernel Matrix for Two Vectors
#'
#' This function calculates the centered kernel matrix between two input vectors,
#' based on the Brownian kernel with a specified Hurst coefficient. The resulting matrix
#' contains the centered kernel values for each pair of elements from the input vectors.
#'
#' @param first_vec_kernel Numeric vector. The first input vector.
#' @param second_vec_kernel Numeric vector. The second input vector.
#' @param centering_grid Numeric vector. The centering grid used for kernel centering.
#' @param hurst_coef Numeric. The Hurst coefficient for the Brownian kernel.
#'
#' @return A matrix of size \code{length(first_vec_kernel)} by \code{length(second_vec_kernel)},
#' where each element represents the centered kernel between the corresponding elements
#' of the first and second vectors.
#' @export
#'
#' @examples
#' vec1 <- c(1, 2, 7)
#' vec2 <- c(1, 2)
#' cent_grid <- sort(runif(1000, min = -1, max = 1))
#' h_result_matrix <- centered_kernel_matrix(vec1, vec2, cent_grid, hurst_coef = 0.5)
centered_kernel_matrix <- function(first_vec_kernel, second_vec_kernel, centering_grid, hurst_coef) {
  # Length of the first input vector
  n0 <- length(first_vec_kernel)

  # Length of the second input vector
  n1 <- length(second_vec_kernel)

  # Length of the centering grid
  n2 <- length(centering_grid)

  # Compute the difference matrix for the first and second vectors raised to the power of 2 * hurst_coef
  term1_matrix <- outer(first_vec_kernel, second_vec_kernel, FUN = function(x, y) abs(x - y)^(2 * hurst_coef))

  # Compute the difference matrix between the first vector and the centering grid
  term2_matrix <- outer(first_vec_kernel, centering_grid, FUN = function(x, y) abs(x - y)^(2 * hurst_coef))

  # Compute the difference matrix between the second vector and the centering grid
  term3_matrix <- outer(second_vec_kernel, centering_grid, FUN = function(x, y) abs(x - y)^(2 * hurst_coef))

  # Compute the difference matrix within the centering grid itself
  term4_matrix <- outer(centering_grid, centering_grid, FUN = function(x, y) abs(x - y)^(2 * hurst_coef))

  # Calculate the centered kernel matrix
  result_matrix <- -1/2 * (term1_matrix - outer(rowMeans(term2_matrix), rep(1, n1)) -
                             outer(rep(1, n0), rowMeans(term3_matrix)) +
                             ones(n0, n1) * mean(term4_matrix))

  return(result_matrix)
}



# Example usage:
#vec1 <- c(1, 2, 7)
#vec2 <- c(1,2)
#cent_grid <- sort(runif(1000,min = -1,max =1))

# Call the function
#h_result_matrix <- centered_kernel_matrix(vec1, vec2,cent_grid)

# Print the result
#print(h_result_matrix)


#c(5,4) %*% h_result_matrix


#term1_matrix <-outer(sampled_x, x_grid, FUN = function(x, y) abs(x - y)^(1))
#term2_matrix <-outer(sampled_x, centering_grid, FUN = function(x, y) abs(x - y)^(1))
#term3_matrix <-outer(x_grid, centering_grid, FUN = function(x, y) abs(x - y)^(1))
#term4_matrix <-outer(centering_grid, centering_grid, FUN = function(x, y) abs(x - y)^(2 * 0.5))

#sum_terms <- matrix(0,nrow = 5, ncol = 10)

#sum_terms <- sum_terms + term1_matrix - outer(rowMeans(term2_matrix),rep(1, 10)) -
#  outer(rep(1, 5),rowMeans(term3_matrix)) + ones(5,10) * mean(term4_matrix)


#' Parallely Compute a Centered Kernel Matrix
#'
#' This function computes a centered kernel matrix using input vectors and a centering grid.
#' The computation involves multiple difference matrices raised to a power determined by the Hurst coefficient.
#'
#' @param first_vec_kernel A numeric vector representing the first set of input points for the kernel computation.
#' @param second_vec_kernel A numeric vector representing the second set of input points for the kernel computation.
#' @param centering_grid A numeric vector used as the centering grid for adjusting the kernel matrix.
#' @param hurst_coef A numeric value (typically between 0 and 1) representing the Hurst coefficient, which determines
#' the scaling of the absolute differences in the kernel computations.
#'
#' @return A numeric matrix of size \code{length(first_vec_kernel) x length(second_vec_kernel)}, representing the
#' centered kernel matrix after all computations.
#'
#' @export
#'
#' @examples
#' # Example usage of centered_kernel_matrix_parallel
#' first_vec_kernel <- seq(1, 10, length.out = 100)
#' second_vec_kernel <- seq(1, 5, length.out = 50)
#' centering_grid <- seq(0, 2, length.out = 20)
#' hurst_coef <- 0.5
#'
#' # Compute the centered kernel matrix in parallel
#' result <- centered_kernel_matrix_parallel(first_vec_kernel, second_vec_kernel, centering_grid, hurst_coef)
#' print(result)

library(parallel)
centered_kernel_matrix_parallel <- function(first_vec_kernel, second_vec_kernel, centering_grid, hurst_coef) {
  # Lengths of input vectors
  n0 <- length(first_vec_kernel)
  n1 <- length(second_vec_kernel)
  n2 <- length(centering_grid)

  # Create a cluster
  num_cores <- detectCores() - 1
  cl <- makeCluster(num_cores)

  # Export necessary variables and functions to the cluster
  clusterExport(cl, c("first_vec_kernel", "second_vec_kernel", "centering_grid", "hurst_coef"), envir = environment())

  # Define helper functions for each matrix computation
  compute_term1_chunk <- function(indices) {
    outer(first_vec_kernel[indices], second_vec_kernel, FUN = function(x, y) abs(x - y)^(2 * hurst_coef))
  }

  compute_term2_chunk <- function(indices) {
    outer(first_vec_kernel[indices], centering_grid, FUN = function(x, y) abs(x - y)^(2 * hurst_coef))
  }

  compute_term3_chunk <- function(indices) {
    outer(second_vec_kernel[indices], centering_grid, FUN = function(x, y) abs(x - y)^(2 * hurst_coef))
  }

  compute_term4_chunk <- function(indices) {
    outer(centering_grid[indices], centering_grid, FUN = function(x, y) abs(x - y)^(2 * hurst_coef))
  }

  # Compute term1_matrix in parallel
  chunk_indices1 <- split(seq_len(n0), cut(seq_len(n0), num_cores, labels = FALSE))
  term1_chunks <- parLapply(cl, chunk_indices1, compute_term1_chunk)
  term1_matrix <- do.call(rbind, term1_chunks)

  # Compute term2_matrix in parallel
  term2_chunks <- parLapply(cl, chunk_indices1, compute_term2_chunk)
  term2_matrix <- do.call(rbind, term2_chunks)

  # Compute term3_matrix in parallel
  chunk_indices2 <- split(seq_len(n1), cut(seq_len(n1), num_cores, labels = FALSE))
  term3_chunks <- parLapply(cl, chunk_indices2, function(indices) {
    outer(second_vec_kernel[indices], centering_grid, FUN = function(x, y) abs(x - y)^(2 * hurst_coef))
  })
  term3_matrix <- do.call(rbind, term3_chunks)

  # Compute term4_matrix in parallel
  chunk_indices3 <- split(seq_len(n2), cut(seq_len(n2), num_cores, labels = FALSE))
  term4_chunks <- parLapply(cl, chunk_indices3, compute_term4_chunk)
  term4_matrix <- do.call(rbind, term4_chunks)

  # Final calculations
  result_matrix <- -1 / 2 * (
    term1_matrix -
      outer(rowMeans(term2_matrix), rep(1, n1)) -
      outer(rep(1, n0), rowMeans(term3_matrix)) +
      matrix(mean(term4_matrix), nrow = n0, ncol = n1)
  )

  # Stop the cluster
  stopCluster(cl)

  return(result_matrix)
}
