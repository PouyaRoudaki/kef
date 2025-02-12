library(pracma)

#' Compute Densities at Sampled Points
#'
#' This function computes the densities at sampled points based on a centered kernel matrix,
#' estimated lambda, and a weight vector.
#'
#' @param centered_kernel_mat_at_sampled A square matrix (n x n) where n is the number of sampled points.
#'        The matrix should be a centered kernel matrix evaluated at the sampled points.
#' @param lambda_hat A scalar value representing the estimated lambda.
#' @param weight_hat_vec A vector of weights (length n) corresponding to the sampled points.
#'
#' @return A vector of densities (length n) corresponding to the sampled points.
#' @export
#'
#' @examples
#' # Example usage (assuming centered_kernel_mat_at_sampled, lambda_hat, and weight_hat_vec are defined):
#' den_vec <- density_at_sampled_x(centered_kernel_mat_at_sampled, lambda_hat, weight_hat_vec)
density_at_sampled_x <- function(centered_kernel_mat_at_sampled, lambda_hat, weight_hat_vec) {
  # Extract the diagonal elements of the centered kernel matrix
  diag_vals <- diag(centered_kernel_mat_at_sampled)

  # Compute the density for each sampled point using the specified formula
  d_vec <- sapply(1:nrow(centered_kernel_mat_at_sampled), function(i) {
    exp(lambda_hat * (weight_hat_vec %*% centered_kernel_mat_at_sampled[,i] - 0.5 * diag_vals[i]))
  })

  return(d_vec)
}

#' Compute Densities at Grid Points
#'
#' This function computes the densities at grid points based on a centered kernel matrix,
#' estimated lambda, and a weight vector.
#'
#' @param centered_kernel_mat_at_grid A matrix (n x m) where n is the number of sampled points
#'        and m is the number of grid points. The matrix should be a centered kernel matrix
#'        evaluated at the grid points.
#' @param centered_kernel_self_grid A vector (length m) representing the diagonal of the centered kernel
#'        matrix evaluated at the grid points.
#' @param lambda_hat A scalar value representing the estimated lambda.
#' @param weight_hat_vec A vector of weights (length n) corresponding to the sampled points.
#'
#' @return A vector of densities (length m) corresponding to the grid points.
#' @export
#'
#' @examples
#' # Example usage (assuming centered_kernel_mat_at_grid, centered_kernel_self_grid, lambda_hat, and weight_hat_vec are defined):
#' den_vec <- densities_at_grid(centered_kernel_mat_at_grid, centered_kernel_self_grid, lambda_hat, weight_hat_vec)
densities_at_grid <- function(centered_kernel_mat_at_grid, centered_kernel_self_grid, lambda_hat, weight_hat_vec) {
  # Compute the density for each grid point using the specified formula
  d_vec <- sapply(1:ncol(centered_kernel_mat_at_grid), function(i) {
    exp(lambda_hat * (weight_hat_vec %*% centered_kernel_mat_at_grid[,i] - 0.5 * centered_kernel_self_grid[i]))
  })

  return(d_vec)
}

#' Compute Densities for Sampled and Grid Points
#'
#' This function computes densities for both sampled points and grid points.
#' The normalization is done by dividing the densities by the integral over the grid points.
#'
#' @param centered_kernel_mat_at_sampled A square matrix (n x n) where n is the number of sampled points.
#'        The matrix should be a centered kernel matrix evaluated at the sampled points.
#' @param centered_kernel_mat_at_grid A matrix (n x m) where n is the number of sampled points
#'        and m is the number of grid points. The matrix should be a centered kernel matrix
#'        evaluated at the grid points.
#' @param centered_kernel_self_grid A vector (length m) representing the diagonal of the centered kernel
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
#'                    centered_kernel_self_grid, x_grid, lambda_hat, weight_hat_vec)
get_densities <- function(centered_kernel_mat_at_sampled, centered_kernel_mat_at_grid,
                      centered_kernel_self_grid, x_grid, lambda_hat, weight_hat_vec){

  # Compute the probabilities at the sampled points
  den_sampled_x <- density_at_sampled_x(centered_kernel_mat_at_sampled, lambda_hat, weight_hat_vec)

  # Compute the densities at the grid points
  den_grid <- density_at_grid(centered_kernel_mat_at_grid, centered_kernel_self_grid, lambda_hat, weight_hat_vec)

  # Normalize the density by the integral over the grid
  normalizing_cte <- trapz(x_grid, den_grid)  # trapz is from pracma package

  # Prepare the output as a list of normalized probabilities
  den_list <- list()
  den_list$sampled_x <- den_sampled_x / normalizing_cte
  den_list$grid_x <- den_grid / normalizing_cte

  return(den_list)
}
