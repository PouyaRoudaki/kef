#' Compute Density at Sampled Points (Fast C++ Version)
#' @param centered_kernel_mat_at_sampled Matrix of centered kernel values.
#' @param lambda_hat Scalar lambda parameter.
#' @param weight_hat_vec Vector of weights.
#' @return Density values at sampled points.
#' @export
density_at_sampled_x <- function(centered_kernel_mat_at_sampled, lambda_hat, weight_hat_vec) {
  .Call(`_kef_density_at_sampled_x`, centered_kernel_mat_at_sampled, lambda_hat, weight_hat_vec)
}

#' Compute Density at Grid Points (Fast C++ Version)
#' @param centered_kernel_mat_at_grid Matrix of kernel values at grid points.
#' @param centered_kernel_self_grid Self-kernel values at grid points.
#' @param lambda_hat Scalar lambda parameter.
#' @param weight_hat_vec Vector of weights.
#' @return Density values at grid points.
#' @export
density_at_grid <- function(centered_kernel_mat_at_grid, centered_kernel_self_grid, lambda_hat, weight_hat_vec) {
  .Call(`_kef_density_at_grid`, centered_kernel_mat_at_grid, centered_kernel_self_grid, lambda_hat, weight_hat_vec)
}

#' Compute Densities for Sampled Points Without a Grid (Fast C++ Version)
#' @param centered_kernel_mat_at_sampled Centered kernel matrix at sampled points.
#' @param min_x Minimum x value.
#' @param max_x Maximum x value.
#' @param sampled_x Sampled x values.
#' @param lambda_hat Scalar lambda parameter.
#' @param weight_hat_vec Vector of weights.
#' @return Normalized densities at sampled points.
#' @export
get_dens_wo_grid <- function(centered_kernel_mat_at_sampled, min_x, max_x, sampled_x, lambda_hat, weight_hat_vec) {
  .Call(`_kef_get_dens_wo_grid`, centered_kernel_mat_at_sampled, min_x, max_x, sampled_x, lambda_hat, weight_hat_vec)
}

#' Compute Densities for Sampled and Grid Points (Fast C++ Version)
#' @param centered_kernel_mat_at_sampled Centered kernel matrix at sampled points.
#' @param centered_kernel_mat_at_grid Centered kernel matrix at grid points.
#' @param centered_kernel_self_grid Self-kernel values at grid points.
#' @param sampled_x Sampled x values.
#' @param x_grid Grid x values.
#' @param lambda_hat Scalar lambda parameter.
#' @param weight_hat_vec Vector of weights.
#' @param type_of_p_is_prob Logical: if TRUE, return probabilities.
#' @param type_of_q_is_prob Logical: if TRUE, return probabilities.
#' @return List with densities at sampled and grid points.
#' @export
get_dens_or_prob <- function(centered_kernel_mat_at_sampled, centered_kernel_mat_at_grid, centered_kernel_self_grid, sampled_x, x_grid, lambda_hat, weight_hat_vec, type_of_p_is_prob = TRUE, type_of_q_is_prob = TRUE) {
  .Call(`_kef_get_dens_or_prob`, centered_kernel_mat_at_sampled, centered_kernel_mat_at_grid, centered_kernel_self_grid, sampled_x, x_grid, lambda_hat, weight_hat_vec, type_of_p_is_prob, type_of_q_is_prob)
}
