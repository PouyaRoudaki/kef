#' Leave-One-Out Cross-Validation (LOOCV) Weight Error Grid
#'
#' @param weight_hat A vector of estimated weights.
#' @param centered_kernel_mat_at_sampled A centered kernel matrix at sampled points.
#' @param centered_kernel_mat_at_grid A centered kernel matrix at grid points.
#' @param centered_kernel_self_grid A centered kernel matrix for the grid points themselves.
#' @param x_grid A grid of x values where the function is evaluated.
#' @param lambda_hat_grid A vector of lambda hat values.
#' @param tau_hat_grid A vector of tau hat values.
#'
#' @return A data frame containing lambda_hat, tau_hat, and the corresponding LOOCV error (loocv_err).
#' @export
#'
#' @examples
#' # Example usage (assuming appropriate variables are defined):
#' # result <- loocv_weight_error_grid(weight_hat, centered_kernel_mat_at_sampled, centered_kernel_mat_at_grid, centered_kernel_self_grid, x_grid, lambda_hat_grid, tau_hat_grid)
loocv_weight_error_grid <- function(weight_hat,
                                    centered_kernel_mat_at_sampled,
                                    centered_kernel_mat_at_grid,
                                    centered_kernel_self_grid,
                                    sampled_x,
                                    x_grid,
                                    lambda_hat_grid,
                                    tau_hat_grid,
                                    type_of_p_is_prob=TRUE,
                                    type_of_q_is_prob=TRUE,
                                    method_of_p_calculation="ordinary") {

  # Create a grid of lambda_hat and tau_hat values
  grid <- expand.grid(lambda_hat = lambda_hat_grid, tau_hat = tau_hat_grid)

  # Compute LOOCV error for each combination in the grid using mapply
  err <- mapply(function(lambda_hat, tau_hat) {

    # Calculate LOOCV error for a given (lambda_hat, tau_hat) pair
    loocv_err <- sum(sapply(1:nrow(centered_kernel_mat_at_sampled), function(i) {

      # Remove the i-th sample from the matrices for LOOCV
      temp_centered_kernel_mat_at_sampled <- centered_kernel_mat_at_sampled[-i, -i]
      temp_centered_kernel_mat_at_grid <- centered_kernel_mat_at_grid[-i, ]

      # Get the LOOCV estimated weights by excluding the i-th sample
      w_hat_loocv <- as.numeric(get_weights(lambda_hat = lambda_hat,
                                            tau_hat = tau_hat,
                                            temp_centered_kernel_mat_at_sampled,
                                            temp_centered_kernel_mat_at_grid,
                                            centered_kernel_self_grid,
                                            sampled_x = sampled_x,
                                            x_grid = x_grid,
                                            type_of_p_is_prob=type_of_p_is_prob,
                                            type_of_q_is_prob=type_of_q_is_prob,
                                            method_of_p_calculation=method_of_p_calculation))

      # Calculate the error for the left-out sample
      one_out_err <- t(w_hat_loocv - weight_hat[-i]) %*%
        temp_centered_kernel_mat_at_sampled %*%
        (w_hat_loocv - weight_hat[-i])

      return(one_out_err)
    }))

    return(loocv_err)
  }, grid$lambda_hat, grid$tau_hat)

  # Combine the results into a data frame
  results <- data.frame(lambda_hat = grid$lambda_hat, tau_hat = grid$tau_hat,
                        loocv_err = err)

  return(results)
}
