# R/centered_kernel.R
#' Compute the Centered Kernel Matrix
#'
#' This function computes the centered kernel matrix efficiently using C++.
#' @param first_vec_kernel Numeric vector.
#' @param second_vec_kernel Numeric vector.
#' @param centering_grid Numeric vector.
#' @param hurst_coef Numeric scalar.
#' @return A numeric matrix.
#' @export
centered_kernel_matrix <- function(first_vec_kernel, second_vec_kernel, centering_grid, hurst_coef) {
  .Call(`_kef_centered_kernel_matrix`, first_vec_kernel, second_vec_kernel, centering_grid, hurst_coef)
}
