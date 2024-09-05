#' Title: Compute Marginal Likelihood
#'
#' This function calculates the marginal likelihood for a given set of parameters.
#'
#' @param centered_kernel_mat_at_sampled A matrix of centered kernel values at sampled points.
#' @param centered_kernel_mat_at_grid A matrix of centered kernel values at grid points.
#' @param centered_kernel_self_grid A matrix of centered kernel self-values at grid points.
#' @param x_grid A grid of x values where the function is evaluated.
#' @param p_vec A vector of initial probabilities (default is a vector of 1's).
#' @param lambda_hat A scalar parameter, lambda, to be used in the likelihood computation.
#' @param tau_hat A scalar parameter, tau, to be used in the likelihood computation.
#' @param MC_iterations Number of Monte Carlo iterations to perform.
#'
#' @return The log of the marginal likelihood.
#' @export
#'
marginal_likelihood <- function(centered_kernel_mat_at_sampled,
                                centered_kernel_mat_at_grid,
                                centered_kernel_self_grid,
                                sampled_x,
                                x_grid,
                                p_vec = rep(1,nrow(centered_kernel_mat_at_sampled)),
                                lambda_hat,
                                tau_hat,
                                MC_iterations,
                                type_of_p_is_prob = FALSE,
                                type_of_q_is_prob = FALSE,
                                method_of_p_calculation = "ordinary"){

  # Get the number of sampled points
  n <- nrow(centered_kernel_mat_at_sampled)

  # Replicate p_vec for Monte Carlo iterations and convert it into a matrix
  p_matrix <- matrix(rep(p_vec, each = MC_iterations), nrow = MC_iterations, byrow = FALSE)

  # Sample weights w_i from a normal distribution N(0, p(x_i)/tau)
  w_sampled <- matrix(rnorm(MC_iterations * n, mean = 0, sd = sqrt(p_matrix / tau_hat)),
                      nrow = MC_iterations, ncol = n)

  # Calculate the probability for each set of sampled weights using the custom 'get_dens_or_prob' function
  probabilities_for_given_weights <- apply(w_sampled, 1, function(w_vec) {
    get_dens_or_prob(centered_kernel_mat_at_sampled, centered_kernel_mat_at_grid,
              centered_kernel_self_grid,
              sampled_x,
              x_grid,
              lambda_hat,
              w_vec,
              type_of_p_is_prob = type_of_p_is_prob,
              type_of_q_is_prob = type_of_q_is_prob,
              method_of_p_calculation = method_of_p_calculation)
  })

  # Extract the probabilities for the sampled x values and combine them into a matrix
  prob_sampled_x_matrix <- do.call(rbind, lapply(probabilities_for_given_weights, `[[`, "sampled_x"))

  # Compute the likelihood vector by taking the mean across Monte Carlo iterations
  likelihood_vector <- colMeans(prob_sampled_x_matrix, na.rm = TRUE)

  # Compute the log of the marginal likelihood by summing the log of the likelihood vector
  marginal_log_likelihood <- sum(log(likelihood_vector))

  # Print a separator (seems to be for debugging purposes)
  #print("-")

  # Return the computed marginal log likelihood
  return(marginal_log_likelihood)
}

#' Title: Compute Marginal Likelihood over a Grid of Parameters
#'
#' This function computes the marginal likelihood for each combination of lambda and tau
#' in the provided grids.
#'
#' @param centered_kernel_mat_at_sampled A matrix of centered kernel values at sampled points.
#' @param centered_kernel_mat_at_grid A matrix of centered kernel values at grid points.
#' @param centered_kernel_self_grid A matrix of centered kernel self-values at grid points.
#' @param x_grid A grid of x values where the function is evaluated.
#' @param p_vec A vector of initial probabilities (default is a vector of 1's).
#' @param lambda_grid A grid of lambda values to be evaluated.
#' @param tau_grid A grid of tau values to be evaluated.
#' @param MC_iterations Number of Monte Carlo iterations to perform.
#'
#' @return A data frame containing lambda, tau, and their corresponding marginal log likelihoods.
#' @export
#'
compute_marginal_likelihood_grid <- function(centered_kernel_mat_at_sampled,
                                             centered_kernel_mat_at_grid,
                                             centered_kernel_self_grid,
                                             sampled_x,
                                             x_grid,
                                             p_vec = rep(1, nrow(centered_kernel_mat_at_sampled)),
                                             lambda_grid,
                                             tau_grid,
                                             MC_iterations,
                                             type_of_p_is_prob = FALSE,
                                             type_of_q_is_prob = FALSE,
                                             method_of_p_calculation = "ordinary") {
  # Create a grid of all combinations of lambda_hat and tau_hat
  grid <- expand.grid(lambda_hat = lambda_grid, tau_hat = tau_grid)

  # Compute the marginal likelihood for each combination in the grid using mapply
  marginal_log_likelihoods <- mapply(function(lambda_hat, tau_hat) {
    marginal_likelihood(centered_kernel_mat_at_sampled,
                        centered_kernel_mat_at_grid,
                        centered_kernel_self_grid,
                        sampled_x,
                        x_grid,
                        p_vec,
                        lambda_hat,
                        tau_hat,
                        MC_iterations,
                        type_of_p_is_prob,
                        type_of_q_is_prob,
                        method_of_p_calculation)
  }, grid$lambda_hat, grid$tau_hat)

  # Combine the results into a data frame
  results <- data.frame(lambda_hat = grid$lambda_hat, tau_hat = grid$tau_hat, marginal_log_likelihood = marginal_log_likelihoods)

  # Return the results
  return(results)
}
