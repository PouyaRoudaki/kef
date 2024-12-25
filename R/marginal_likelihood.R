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
                                sampled_x,
                                min_x,
                                max_x,
                                p_vec = rep(1,nrow(centered_kernel_mat_at_sampled)),
                                lambda_hat,
                                tau_hat,
                                MC_iterations){

  # Get the number of sampled points
  n <- nrow(centered_kernel_mat_at_sampled)

  # Sample weights w_i from a normal distribution N(0, p(x_i)/tau)
  w_sampled <- matrix(rnorm(MC_iterations * n, mean = 0, sd = sqrt(p_vec / tau_hat)),
                      byrow = TRUE,
                      nrow = MC_iterations, ncol = n)

  #data <- as.data.frame(cbind(rep(1:length(sampled_x),each = MC_iterations),
  #                            as.numeric(w_sampled)))

  #data <- data %>% mutate(sampled_x = sampled_x[V1])

  #data1 <- data %>%
  #  group_by(V1) %>%
  #  summarise(var(V2))

  #print(plot(sampled_x, sqrt(data1$`var(V2)`)))

  #print(plot(sampled_x,w_sampled))
  # Calculate the probability for each set of sampled weights using the custom 'get_dens_wo_grid' function
  probabilities_for_given_weights <- apply(w_sampled, 1, function(w_vector) {
    get_dens_wo_grid(centered_kernel_mat_at_sampled,
                     min_x,
                     max_x,
                     sampled_x,
                     lambda_hat,
                     w_vector)
  })

  #print(dim(probabilities_for_given_weights))

  sample_mid_points <- get_middle_points_grid(min_x, sampled_x, max_x)
  base_measure_weights <- sample_mid_points[-1] - sample_mid_points[-length(sample_mid_points)]

  prob_inverse_base <- probabilities_for_given_weights * exp(-base_measure_weights)


  # Compute the likelihood vector by taking the mean across Monte Carlo iterations
  likelihood_vector <- rowMeans(prob_inverse_base, na.rm = TRUE)

  #print(length(likelihood_vector))
  #print(plot(sampled_x,likelihood_vector))

  # Create a data frame with your data
  plot_data <- data.frame(x = sampled_x, likelihood = likelihood_vector)

  # Plot using ggplot
  plot <- ggplot(plot_data, aes(x = x, y = likelihood)) +
    geom_line() +  # Line plot
    geom_point() +  # Add points to emphasize sampled values
    labs(title = paste("Likelihood vs. Sampled X, lambda = ", lambda_hat),
         x = "Sampled X",
         y = "Likelihood") +
    theme_bw()  # Use a clean theme

  print(plot)

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
                                             min_x,
                                             max_x,
                                             sampled_x,
                                             lambda_grid,
                                             tau_hat,
                                             initial_lambda_hat = 1,
                                             MC_iterations) {


  # Compute the marginal likelihood for each combination in the grid using mapply
  counter <- 1

  lambda_hat <- initial_lambda_hat
  while (counter <= 30) {

    w_vec <- as.numeric(get_weights_wo_grid(
      lambda_hat,
      tau_hat,
      centered_kernel_mat_at_sampled,
      sampled_x,
      min_x,
      max_x))

    p_vec <- get_dens_wo_grid(centered_kernel_mat_at_sampled,
                              min_x,
                              max_x,
                              sampled_x,
                              lambda_hat,
                              w_vec)

    max_marginal_log_likelihood <- -Inf

    for (lambda_hat in lambda_grid) {

      cat("-")

      #print(length(p_vec))
      marginal_log_likelihood <- marginal_likelihood(
      centered_kernel_mat_at_sampled,
                        sampled_x,
                        min_x,
                        max_x,
                        p_vec,
                        lambda_hat,
                        tau_hat,
                        MC_iterations)

      #print(lambda_hat)
      #print(marginal_log_likelihood)
      #print(max_marginal_log_likelihood)


      if(marginal_log_likelihood > max_marginal_log_likelihood){
        max_marginal_log_likelihood <- marginal_log_likelihood
        max_likelihood_lambda <- lambda_hat
      }


    }

    lambda_hat <- max_likelihood_lambda
    # Add an extra newline after the loop
    cat("\n")
    print(paste0("Iteration = ", counter,
                   ", lambda_hat = ", lambda_hat,
                   ", tau_hat = ", tau_hat,
                   ", max_marginal_log_likelihood = ",max_marginal_log_likelihood))

    counter <- counter +1
  }
}
