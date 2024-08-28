library(pracma)

#' Compute Probabilities at Sampled Points
#'
#' This function computes the probabilities at sampled points based on a centered kernel matrix,
#' estimated lambda, and a weight vector.
#'
#' @param centered_kernel_mat_at_sampled A square matrix (n x n) where n is the number of sampled points.
#'        The matrix should be a centered kernel matrix evaluated at the sampled points.
#' @param lambda_hat A scalar value representing the estimated lambda.
#' @param weight_hat_vec A vector of weights (length n) corresponding to the sampled points.
#'
#' @return A vector of probabilities (length n) corresponding to the sampled points.
#' @export
#'
#' @examples
#' # Example usage (assuming centered_kernel_mat_at_sampled, lambda_hat, and weight_hat_vec are defined):
#' prob_vec <- prob_at_sampled_x(centered_kernel_mat_at_sampled, lambda_hat, weight_hat_vec)
prob_at_sampled_x <- function(centered_kernel_mat_at_sampled, lambda_hat, weight_hat_vec) {
  # Extract the diagonal elements of the centered kernel matrix
  diag_vals <- diag(centered_kernel_mat_at_sampled)

  # Compute the probability for each sampled point using the specified formula
  p_vec <- sapply(1:nrow(centered_kernel_mat_at_sampled), function(i) {
    exp(lambda_hat * (weight_hat_vec %*% centered_kernel_mat_at_sampled[,i] - 0.5 * diag_vals[i]))
  })

  return(p_vec)
}

#' Compute Probabilities at Grid Points
#'
#' This function computes the probabilities at grid points based on a centered kernel matrix,
#' estimated lambda, and a weight vector.
#'
#' @param centered_kernel_mat_at_grid A matrix (n x m) where n is the number of sampled points
#'        and m is the number of grid points. The matrix should be a centered kernel matrix
#'        evaluated at the grid points.
#' @param centerd_kernel_self_grid A vector (length m) representing the diagonal of the centered kernel
#'        matrix evaluated at the grid points.
#' @param lambda_hat A scalar value representing the estimated lambda.
#' @param weight_hat_vec A vector of weights (length n) corresponding to the sampled points.
#'
#' @return A vector of probabilities (length m) corresponding to the grid points.
#' @export
#'
#' @examples
#' # Example usage (assuming centered_kernel_mat_at_grid, centerd_kernel_self_grid, lambda_hat, and weight_hat_vec are defined):
#' prob_vec <- prob_at_grid(centered_kernel_mat_at_grid, centerd_kernel_self_grid, lambda_hat, weight_hat_vec)
prob_at_grid <- function(centered_kernel_mat_at_grid, centerd_kernel_self_grid, lambda_hat, weight_hat_vec) {
  # Compute the probability for each grid point using the specified formula
  p_vec <- sapply(1:ncol(centered_kernel_mat_at_grid), function(i) {
    exp(lambda_hat * (weight_hat_vec %*% centered_kernel_mat_at_grid[,i] - 0.5 * centerd_kernel_self_grid[i]))
  })

  return(p_vec)
}

#' Compute Normalized Probabilities for Sampled and Grid Points
#'
#' This function computes normalized probabilities for both sampled points and grid points.
#' The normalization is done by dividing the probabilities by the integral over the grid points.
#'
#' @param centered_kernel_mat_at_sampled A square matrix (n x n) where n is the number of sampled points.
#'        The matrix should be a centered kernel matrix evaluated at the sampled points.
#' @param centered_kernel_mat_at_grid A matrix (n x m) where n is the number of sampled points
#'        and m is the number of grid points. The matrix should be a centered kernel matrix
#'        evaluated at the grid points.
#' @param centerd_kernel_self_grid A vector (length m) representing the diagonal of the centered kernel
#'        matrix evaluated at the grid points.
#' @param x_grid A vector of grid points where the probabilities are evaluated.
#' @param lambda_hat A scalar value representing the estimated lambda.
#' @param weight_hat_vec A vector of weights (length n) corresponding to the sampled points.
#'
#' @return A list containing two elements:
#'         \item{sampled_x}{A vector of normalized probabilities at the sampled points.}
#'         \item{grid_x}{A vector of normalized probabilities at the grid points.}
#' @export
#'
#' @examples
#' # Example usage (assuming inputs are defined):
#' probs <- get_probs(centered_kernel_mat_at_sampled, centered_kernel_mat_at_grid,
#'                    centerd_kernel_self_grid, x_grid, lambda_hat, weight_hat_vec)
get_probs <- function(centered_kernel_mat_at_sampled, centered_kernel_mat_at_grid,
                      centerd_kernel_self_grid, x_grid, lambda_hat, weight_hat_vec){

  # Compute the probabilities at the sampled points
  prob_sampled_x <- prob_at_sampled_x(centered_kernel_mat_at_sampled, lambda_hat, weight_hat_vec)

  # Compute the probabilities at the grid points
  prob_grid <- prob_at_grid(centered_kernel_mat_at_grid, centerd_kernel_self_grid, lambda_hat, weight_hat_vec)

  # Normalize the probabilities by the integral over the grid
  normalizing_cte <- trapz(x_grid, prob_grid)  # trapz is from pracma package

  # Prepare the output as a list of normalized probabilities
  prob_list <- list()
  prob_list$sampled_x <- prob_sampled_x / normalizing_cte
  prob_list$grid_x <- prob_grid / normalizing_cte

  return(prob_list)
}
