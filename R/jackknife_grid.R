#' Calculate Jackknife Error for a Given Pair of Parameters
#'
#' Computes the jackknife error
#' given specific values of \code{lambda_hat} and \code{tau_hat}.
#'
#' @param lambda_hat Numeric. Regularization parameter for weight estimation.
#' @param tau_hat Numeric. Regularization parameter for controlling smoothness.
#' @param centered_kernel_mat_at_sampled Centered kernel matrix at sampled points.
#' @param centered_kernel_mat_at_grid Centered kernel matrix at grid points.
#' @param centered_kernel_self_grid Centered kernel matrix for the grid points themselves.
#' @param sampled_x Sampled x values.
#' @param x_grid Grid of x values where the function is evaluated.
#' @param type_of_p_is_prob Logical. If TRUE, p is treated as probabilities.
#' @param type_of_q_is_prob Logical. If TRUE, q is treated as probabilities.
#' @param method_of_p_calculation Character. Method for p calculation.
#'
#' @return Numeric. The calculated jackknife error for the given parameter pair.
#' @export
#'
#' @examples
#' # Example usage:
#' # calc_jackknife_error(0.1, 0.5, weight_hat, ...)
calc_jackknife_error <- function(lambda_hat,
                                 tau_hat,
                                 centered_kernel_mat_at_sampled,
                                 centered_kernel_mat_at_grid,
                                 centered_kernel_self_grid,
                                 sampled_x,
                                 x_grid,
                                 type_of_p_is_prob,
                                 type_of_q_is_prob,
                                 method_of_p_calculation) {

  weights_hat <- get_weights(lambda_hat =lambda_hat,
                             tau_hat = tau_hat,
                             centered_kernel_mat_at_sampled,
                             centered_kernel_mat_at_grid,
                             centered_kernel_self_grid,
                             sampled_x = sampled_x,
                             x_grid = x_grid,
                             type_of_p_is_prob,
                             type_of_q_is_prob,
                             method_of_p_calculation)

  jackknife_err <- sum(sapply(1:nrow(centered_kernel_mat_at_sampled), function(i) {
    # Remove the i-th sample
    temp_centered_kernel_mat_at_sampled <- centered_kernel_mat_at_sampled[-i, -i]
    temp_centered_kernel_mat_at_grid <- centered_kernel_mat_at_grid[-i, ]
    temp_sampled_x <- sampled_x[-i]

    # Get jackknife estimated weights
    w_hat_jackknife <- as.numeric(get_weights(lambda_hat = lambda_hat,
                                              tau_hat = tau_hat,
                                              temp_centered_kernel_mat_at_sampled,
                                              temp_centered_kernel_mat_at_grid,
                                              centered_kernel_self_grid,
                                              sampled_x = temp_sampled_x,
                                              x_grid = x_grid,
                                              type_of_p_is_prob = type_of_p_is_prob,
                                              type_of_q_is_prob = type_of_q_is_prob,
                                              method_of_p_calculation = method_of_p_calculation))

    # Compute jackknife error for the i-th sample
    one_out_err <- t(w_hat_jackknife - weights_hat[-i]) %*%
      temp_centered_kernel_mat_at_sampled %*%
      (w_hat_jackknife - weights_hat[-i])

    return(one_out_err)
  }))

  return(jackknife_err)
}



#' Jackknife Weight Error Grid
#'
#' Computes jackknife error for a grid of lambda_hat and tau_hat values.
#'
#' @param centered_kernel_mat_at_sampled Centered kernel matrix at sampled points.
#' @param centered_kernel_mat_at_grid Centered kernel matrix at grid points.
#' @param centered_kernel_self_grid Centered kernel matrix for the grid points themselves.
#' @param sampled_x Sampled x values.
#' @param x_grid Grid of x values where the function is evaluated.
#' @param lambda_hat_grid A vector of lambda hat values.
#' @param tau_hat_grid A vector of tau hat values.
#' @param type_of_p_is_prob Logical. If TRUE, p is treated as probabilities.
#' @param type_of_q_is_prob Logical. If TRUE, q is treated as probabilities.
#' @param method_of_p_calculation Character. Method for p calculation.
#'
#' @return A data frame containing lambda_hat, tau_hat, and the corresponding Jackknife error (Jackknife_err).
#' @export
#'
#' @examples
#' # Example usage:
#' # jackknife_weight_error_grid(...)
jackknife_weight_error_grid <- function(centered_kernel_mat_at_sampled,
                                        centered_kernel_mat_at_grid,
                                        centered_kernel_self_grid,
                                        sampled_x,
                                        x_grid,
                                        lambda_hat_grid,
                                        tau_hat_grid,
                                        type_of_p_is_prob = TRUE,
                                        type_of_q_is_prob = TRUE,
                                        method_of_p_calculation = "ordinary") {
  # Create a grid of lambda_hat and tau_hat values
  grid <- expand.grid(lambda_hat = lambda_hat_grid, tau_hat = tau_hat_grid)

  # Create a cluster
  num_cores <- parallel::detectCores() - 1
  cl <- parallel::makeCluster(num_cores)


  #func_code <- parallel::clusterEvalQ(cl, deparse(get_weights))
  #print(func_code)
  # Export necessary objects and functions to the cluster
  parallel::clusterExport(cl,
                          varlist = c("get_weights", "calc_jackknife_error",
                                      "centered_kernel_mat_at_sampled",
                                      "centered_kernel_mat_at_grid",
                                      "centered_kernel_self_grid",
                                      "sampled_x", "x_grid",
                                      "type_of_p_is_prob", "type_of_q_is_prob",
                                      "method_of_p_calculation"),
                          envir = environment())

  # Source or reload the updated function in each worker to avoid stale versions
  parallel::clusterEvalQ(cl, devtools::load_all())


  # Compute jackknife error for each combination in the grid using parallelization
  err <- unlist(parallel::parLapply(cl, seq_len(nrow(grid)), function(idx) {
    lambda_hat <- grid$lambda_hat[idx]
    tau_hat <- grid$tau_hat[idx]

    calc_jackknife_error(lambda_hat, tau_hat,
                         centered_kernel_mat_at_sampled,
                         centered_kernel_mat_at_grid,
                         centered_kernel_self_grid,
                         sampled_x, x_grid,
                         type_of_p_is_prob,
                         type_of_q_is_prob,
                         method_of_p_calculation)
  }))







  # Stop the cluster
  parallel::stopCluster(cl)

  # Combine the results into a data frame
  results <- data.frame(lambda_hat = grid$lambda_hat, tau_hat = grid$tau_hat,
                        jackknife_err = err)

  return(results)
}

