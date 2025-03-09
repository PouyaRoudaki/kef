#' Compute the Mean Squared Error (MSE) in an RKHS Given the True Weights
#'
#' This function calculates the MSE using a Reproducing Kernel Hilbert Space (RKHS) formulation.
#'
#' @param w_hat_vec Numeric vector. Estimated weights.
#' @param w_vec Numeric vector. True weights.
#' @param kernel_matrix_at_samples Numeric matrix. The kernel matrix evaluated at the sample points.
#'
#' @return A numeric scalar representing the computed MSE.
#' @export
rkhs_se <- function(w_hat_vec, w_vec, kernel_matrix_at_samples){
  se <- t(w_hat_vec - w_vec) %*%
    kernel_matrix_at_samples %*%
    (w_hat_vec - w_vec)
  return(as.numeric(se))
}


#' Compute the Integrated Squared Error (ISE) for L2 Norm Given the True Densities
#'
#' This function calculates the Integrated Squared Error (ISE) between the true and estimated density functions using numerical integration.
#'
#' @param grid Numeric vector. The grid points where the densities are evaluated.
#' @param true_density_grid Numeric vector. The true density values evaluated at `grid` points.
#' @param estimated_density_grid Numeric vector. The estimated density values evaluated at `grid` points.
#'
#' @return A numeric scalar representing the computed ISE.
#' @export
l2_ise <- function(grid, true_density_grid, estimated_density_grid){
  l2_error <- (true_density_grid - estimated_density_grid)^2
  ise <- pracma::trapz(grid, l2_error)
  return(as.numeric(ise))
}



