#' Title
#'
#' @param centered_kernel_mat_at_sampled
#' @param centered_kernel_mat_at_grid
#' @param centerd_kernel_self_grid
#' @param x_grid
#' @param p_old_vec
#' @param lambda_hat
#' @param tau_hat
#' @param MC_iterations
#'
#' @return
#' @export
#'
#' @examples
marginal_likelihood <- function(centered_kernel_mat_at_sampled,
                                centered_kernel_mat_at_grid,
                                centerd_kernel_self_grid,
                                x_grid,
                                p_old_vec = rep(1,nrow(centered_kernel_mat_at_sampled)),
                                lambda_hat,
                                tau_hat,
                                MC_iterations){

  # initials
  n <- nrow(centered_kernel_mat_at_sampled)

  # Take number of iteration in Monte Carlo approximation
  #MC_iterations <- 1000

  # Take probability of p_old(x_i)
  #p_old_vec <- rep(1,n)
  p_old_matrix <- matrix(rep(p_old_vec, each = MC_iterations), nrow = MC_iterations, byrow = F)

  # Take weights w_i from Normal distribution N(0, p(x_i)/\tau)
  w_sampled <- matrix(rnorm(MC_iterations * n, mean = 0, sd = sqrt(p_old_matrix / tau_hat)),
                      nrow = MC_iterations, ncol = n)


  # Apply get_probs function to each column of w_sampled and get the probability
  probabilities_for_given_weights <- apply(w_sampled, 1, function(w_vec) {
    get_probs(centered_kernel_mat_at_sampled, centered_kernel_mat_at_grid,
              centerd_kernel_self_grid, x_grid, lambda_hat, w_vec)
  })


  # Extract the sampled_x vectors and combine them into a matrix
  prob_sampled_x_matrix <- do.call(rbind, lapply(probabilities_for_given_weights, `[[`, "sampled_x"))

  # Likelihood vector
  likelihood_vector <- colMeans(prob_sampled_x_matrix,na.rm = T)

  # Marginal likelihood
  marginal_log_likelihood <- sum(log(likelihood_vector))

  print("-")
  # Return
  return(marginal_log_likelihood)


}

# Function to compute marginal likelihood for a grid of lambda_hat and tau_hat values
#' Title
#'
#' @param centered_kernel_mat_at_sampled
#' @param centered_kernel_mat_at_grid
#' @param centerd_kernel_self_grid
#' @param x_grid
#' @param p_old_vec
#' @param lambda_grid
#' @param tau_grid
#' @param MC_iterations
#'
#' @return
#' @export
#'
#' @examples
compute_marginal_likelihood_grid <- function(centered_kernel_mat_at_sampled,
                                             centered_kernel_mat_at_grid,
                                             centerd_kernel_self_grid,
                                             x_grid,
                                             p_old_vec = rep(1, nrow(centered_kernel_mat_at_sampled)),
                                             lambda_grid,
                                             tau_grid,
                                             MC_iterations) {
  # Create a grid of lambda_hat and tau_hat values
  grid <- expand.grid(lambda_hat = lambda_grid, tau_hat = tau_grid)

  # Compute marginal likelihood for each combination in the grid using mapply
  marginal_log_likelihoods <- mapply(function(lambda_hat, tau_hat) {
    marginal_likelihood(centered_kernel_mat_at_sampled,
                        centered_kernel_mat_at_grid,
                        centerd_kernel_self_grid,
                        x_grid,
                        p_old_vec,
                        lambda_hat,
                        tau_hat,
                        MC_iterations)
  }, grid$lambda_hat, grid$tau_hat)

  # Combine results into a data frame
  results <- data.frame(lambda_hat = grid$lambda_hat, tau_hat = grid$tau_hat, marginal_log_likelihood = marginal_log_likelihoods)

  return(results)
}

# Define grids for lambda_hat
#log_min_lambda <- log10(0.04)
#log_max_lambda <- log10(400)

# Number of intervals
#num_intervals <- 10

# Calculate the step size
#step_size <- (log_max_lambda - log_min_lambda) / num_intervals

# Generate the log scale values
#log_lambda_values <- seq(log_min_lambda, log_max_lambda, by = step_size)

# Convert back to the original scale
#lambda_grid <- 10^log_lambda_values

# Define grids for tau_hat
#log_min_tau <- log10(0.0001)
#log_max_tau <- log10(10)

# Number of intervals
#num_intervals <- 100

# Calculate the step size
#step_size <- (log_max_tau - log_min_tau) / num_intervals

# Generate the log scale values
#log_tau_values <- seq(log_min_tau, log_max_tau, by = step_size)

# Convert back to the original scale
#tau_grid <- 10^log_tau_values




# Compute marginal likelihood over the grid
#marginal_likelihood_wrt_tau_lambda <- compute_marginal_likelihood_grid(centered_kernel_mat_at_sampled,
#                                            centered_kernel_mat_at_grid,
#                                            centerd_kernel_self_grid,
#                                            x_grid,
#                                            p_old_vec= rep(1, nrow(centered_kernel_mat_at_sampled)),
#                                            lambda_grid,
#                                            tau_grid,
#                                            MC_iterations = 1000)


save.image(file = "marginal_likelihoods.RData")

#install.packages("plot3D")
#install.packages("rgl")
#install.packages("akima")
library(rgl)
library(akima)


marginal_likelihood_wrt_tau_lambda_1 <- na.omit(marginal_likelihood_wrt_tau_lambda)
# Interpolate to create a grid
interp_result <- with(marginal_likelihood_wrt_tau_lambda_1,
                      interp(lambda_hat,
                             tau_hat,
                             marginal_log_likelihood,
                             duplicate = "mean"))

# Handle NAs in the interpolated grid
interp_result$z[is.na(interp_result$z)] <- mean(interp_result$z, na.rm = TRUE)

# Plot the surface
persp3d(interp_result$x, interp_result$y, interp_result$z, col = "lightblue")

# Add points
points3d(marginal_likelihood_wrt_tau_lambda_1$lambda_hat,
         marginal_likelihood_wrt_tau_lambda_1$tau_hat,
         marginal_likelihood_wrt_tau_lambda_1$marginal_log_likelihood, col = "red", size = 5)


