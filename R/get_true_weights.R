#' Compute Approximation of True Weights for Density Estimation
#'
#' This function estimates the true weights (`w`) used in a density estimation model.
#' It computes the centered kernel matrices, evaluates the density function, and
#' optimizes the weights using least squares error minimization.
#'
#' @param sample A numeric vector of sampled points.
#' @param grid A numeric vector representing the grid over which the density is estimated.
#' @param true_density_at_grid A numeric vector of true density values at each grid point.
#'
#' @return A numeric vector of estimated weights (`w`) corresponding to the sampled points.
#' @importFrom pracma trapz Norm
#' @importFrom stats optim
#'
#' @examples
#' sample <- seq(-3, 3, length.out = 10)
#' grid <- seq(-3, 3, length.out = 100)
#' true_density_at_grid <- dnorm(grid)  # Example: Gaussian density
#' w_estimates <- get_true_weights(sample, grid, true_density_at_grid)
#' print(w_estimates)
#'
#' @export
get_true_weights <- function(sample, grid, true_density_at_grid) {

  # Compute the centered kernel matrix at sampled points
  centered_kernel_mat_at_sampled <- centered_kernel_matrix(
    first_vec_kernel = sample,
    second_vec_kernel = sample,
    centering_grid = grid,
    hurst_coef = 0.5
  )

  # Compute the centered kernel matrix at grid points
  centered_kernel_mat_at_grid <- centered_kernel_matrix(
    first_vec_kernel = sample,
    second_vec_kernel = grid,
    centering_grid = grid,
    hurst_coef = 0.5
  )

  # Compute the diagonal of the centered kernel matrix at grid points
  centered_kernel_self_grid <- diag(centered_kernel_matrix(
    first_vec_kernel = grid,
    second_vec_kernel = grid,
    centering_grid = grid,
    hurst_coef = 0.5
  ))

  # Given true f(x) for all grid points
  f_x <- true_density_at_grid

  # Function to compute predicted f(x) given weights w
  predicted_f <- function(w) {
    numerator <- get_dens_or_prob(
      centered_kernel_mat_at_sampled,
      centered_kernel_mat_at_grid,
      centered_kernel_self_grid,
      sample, grid,
      lambda_hat = 1,
      weight_hat_vec = w,
      type_of_p_is_prob = FALSE,
      type_of_q_is_prob = FALSE
    )$grid_x

    # Compute the denominator using numerical integration
    denom <- pracma::trapz(grid, numerator)

    return(numerator / denom)
  }

  # Define loss function for optimization
  loss_function <- function(w) {
    return((pracma::Norm(predicted_f(w) - f_x, 2))^2)  # Least squares error
  }

  # Choose method for solving the nonlinear system of equations
  initial_w <- rep(0, length(sample))  # Initial guess for weights


  result <- optim(initial_w, loss_function, method = "L-BFGS-B")
  approx_true_w <- result$par



  return(approx_true_w)
}
