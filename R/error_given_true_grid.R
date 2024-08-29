#' RKHS norm of error based on the true weights.
#'
#' @param weight_true A vector of true weights.
#' @param centered_kernel_mat_at_sampled A centered kernel matrix at sampled points.
#' @param centered_kernel_mat_at_grid A centered kernel matrix at grid points.
#' @param centered_kernel_self_grid A centered kernel matrix for the grid points themselves.
#' @param x_grid A grid of x values where the function is evaluated.
#' @param lambda_hat_grid A vector of lambda hat values.
#' @param tau_hat_grid A vector of tau hat values.
#'
#' @return A data frame with lambda_hat, tau_hat, and Norm_Diff_Err values.
#' @export
#'
#' @examples
#' # Example usage (assuming appropriate variables are defined):
#' # result <- weights_error_given_true_grid(weight_true, centered_kernel_mat_at_sampled, centered_kernel_mat_at_grid, centered_kernel_self_grid, x_grid, lambda_hat_grid, tau_hat_grid)
weights_error_given_true_grid <- function(weight_true,
                       centered_kernel_mat_at_sampled,
                       centered_kernel_mat_at_grid,
                       centered_kernel_self_grid,
                       x_grid,
                       lambda_hat_grid,
                       tau_hat_grid) {

  # Create a grid of lambda_hat and tau_hat values
  grid <- expand.grid(lambda_hat = lambda_hat_grid, tau_hat = tau_hat_grid)

  # Compute the norm difference error for each combination in the grid using mapply
  err <- mapply(function(lambda_hat, tau_hat) {
    weight_hat <- get_weights(lambda_hat, tau_hat,
                              centered_kernel_mat_at_sampled,
                              centered_kernel_mat_at_grid,
                              centered_kernel_self_grid,
                              x_grid)

    # Ensure weight_hat is a vector
    weight_hat <- as.vector(weight_hat)

    # Calculate the norm difference error
    Norm_Diff_Err <- t(weight_hat - weight_true) %*%
      centered_kernel_mat_at_sampled %*%
      (weight_hat - weight_true)

    return(Norm_Diff_Err)
  }, grid$lambda_hat, grid$tau_hat)

  # Combine results into a data frame
  results <- data.frame(lambda_hat = grid$lambda_hat, tau_hat = grid$tau_hat,
                        Norm_Diff_Err = err)

  return(results)
}
