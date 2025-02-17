#' Compute Marginal Log Likelihood Using Rcpp
#'
#' This function computes the marginal log likelihood using a Monte Carlo approach
#' with optional parallel computing.
#'
#' @param centered_kernel_mat_at_sampled A matrix of centered kernel values at sampled points.
#' @param sampled_x A numeric vector of sampled points.
#' @param min_x The minimum x value.
#' @param max_x The maximum x value.
#' @param p_vec A probability vector (default: uniform distribution).
#' @param lambda A scalar for the lambda hyperparameter.
#' @param tau A scalar for the tau hyperparameter.
#' @param std_rnorm_matrix A matrix of standard normal random values for Monte Carlo sampling.
#' @param MC_iterations The number of Monte Carlo iterations.
#' @param parallel_computing Logical; if TRUE, enables parallel computing.
#'
#' @return The computed log of the marginal likelihood.
#' @export
marginal_log_likelihood <- function(centered_kernel_mat_at_sampled,
                                    sampled_x,
                                    min_x,
                                    max_x,
                                    p_vec = rep(1, nrow(centered_kernel_mat_at_sampled)),
                                    lambda,
                                    tau,
                                    std_rnorm_matrix,
                                    MC_iterations,
                                    parallel_computing = TRUE) {

  # Call the C++ function using `.Call()`
  .Call("_kef_marginal_log_likelihood",
        centered_kernel_mat_at_sampled,
        sampled_x,
        min_x,
        max_x,
        p_vec,
        lambda,
        tau,
        std_rnorm_matrix,
        MC_iterations,
        parallel_computing)
}


#' Compute Marginal Likelihood Over a Grid of Hyperparameters (Parallel)
#'
#' This function computes the marginal likelihood for each combination of lambda and tau
#' in the provided hyperparameter grid using parallel computing.
#'
#' @param centered_kernel_mat_at_sampled A matrix of centered kernel values at sampled points.
#' @param min_x The minimum x value.
#' @param max_x The maximum x value.
#' @param sampled_x A numeric vector of sampled points.
#' @param hyperparam_grid A dataframe containing pairs of lambda and tau values.
#' @param initial_lambda The initial lambda value.
#' @param initial_w The initial weight vector.
#' @param MC_iterations The number of Monte Carlo iterations.
#' @param max_iterations The maximum number of iterations.
#' @param parallel_computing Logical; if TRUE, enables parallel computing.
#'
#' @return A dataframe containing lambda, tau, and their corresponding marginal log likelihoods.
#' @export
compute_marginal_likelihood_grid_parallel <- function(centered_kernel_mat_at_sampled,
                                                      min_x,
                                                      max_x,
                                                      sampled_x,
                                                      hyperparam_grid,
                                                      initial_lambda = 1,
                                                      initial_w = rep(0, length(sampled_x)),
                                                      MC_iterations,
                                                      max_iterations = 5,
                                                      parallel_computing = TRUE) {

  # Ensure hyperparam_grid is a matrix
  if (!is.matrix(hyperparam_grid) && !is.data.frame(hyperparam_grid)) {
    stop("Error: 'hyperparam_grid' must be a matrix or a dataframe.")
  }

  # Convert dataframe to matrix if needed
  if (is.data.frame(hyperparam_grid)) {
    hyperparam_grid <- as.matrix(hyperparam_grid)
  }

  # Call the C++ function using `.Call()`
  results <- .Call("_kef_compute_marginal_likelihood_grid_parallel",
                   centered_kernel_mat_at_sampled,
                   min_x,
                   max_x,
                   sampled_x,
                   hyperparam_grid,
                   initial_lambda,
                   initial_w,
                   MC_iterations,
                   max_iterations,
                   parallel_computing)

  # Convert the output matrix to a dataframe
  results_df <- as.data.frame(results)
  colnames(results_df) <- c("lambda", "tau", "marginal_log_likelihood")

  return(results_df)
}
