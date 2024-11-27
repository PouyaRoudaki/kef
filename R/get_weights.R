#' Estimate Weights Using the Newton-Raphson Method
#'
#' This function estimates the weight vector using an iterative Newton-Raphson method.
#' The weights are updated based on the provided kernel matrices, regularization parameters,
#' and grid points.
#'
#' @param lambda_hat A scalar representing the estimated lambda parameter.
#' @param tau_hat A scalar representing the estimated tau parameter.
#' @param centered_kernel_mat_at_sampled A square matrix (n x n) representing the centered kernel matrix
#'        evaluated at the sampled points.
#' @param centered_kernel_mat_at_grid A matrix (n x m) representing the centered kernel matrix evaluated
#'        at the grid points, where n is the number of sampled points and m is the number of grid points.
#' @param centered_kernel_self_grid A vector of length m representing the diagonal of the centered kernel
#'        matrix evaluated at the grid points.
#' @param x_grid A vector representing the grid points.
#'
#' @return A numeric vector of length n representing the estimated weight vector.
#' @export
#'
#' @examples
#' # Example usage:
#' lambda_hat <- 1
#' tau_hat <- 0.5
#' n <- 10
#' m <- 20
#' centered_kernel_mat_at_sampled <- matrix(runif(n * n), n, n)
#' centered_kernel_mat_at_grid <- matrix(runif(n * m), n, m)
#' centered_kernel_self_grid <- runif(m)
#' x_grid <- seq(0, 1, length.out = m)
#' weights <- get_weights(lambda_hat, tau_hat, centered_kernel_mat_at_sampled,
#'                        centered_kernel_mat_at_grid, centered_kernel_self_grid, x_grid)
get_weights <- function(lambda_hat, tau_hat, centered_kernel_mat_at_sampled,
                        centered_kernel_mat_at_grid, centered_kernel_self_grid,
                        sampled_x, x_grid,
                        type_of_p_is_prob=TRUE,
                        type_of_q_is_prob=TRUE,
                        method_of_p_calculation="ordinary") {

  max_iteration <- 2000  # Maximum number of iterations for the Newton-Raphson method
  NRstepsize <- 0.1  # Step size for the Newton-Raphson update
  n <- nrow(centered_kernel_mat_at_sampled)  # Number of sampled points
  weight_hat_vec <- rep(0, n)  # Initialize the weight vector with zeros
  #s <- rep(1000, n)
  #weight_hat_change <- rep(1000, n)
  #counter <- 1
  #while ((norm(s, p = 2) > 10^(-10)) & (norm(weight_hat_change, p = 2) > 10^(-10))) {
  for (i in 1:max_iteration) {
    # Calculate probabilities for sampled and grid points
    probs <- get_dens_or_prob(centered_kernel_mat_at_sampled, centered_kernel_mat_at_grid,
                       centered_kernel_self_grid,
                       sampled_x,
                       x_grid,
                       lambda_hat,
                       weight_hat_vec,
                       type_of_p_is_prob,
                       type_of_q_is_prob,
                       method_of_p_calculation)

    prob_sampled_x <- probs$sampled_x
    prob_grid_x <- probs$grid_x

    #print(probs)
    #print(length(probs$grid_x))
    #print(length(prob_grid_x))
    #print(dim(centered_kernel_mat_at_grid))
    # Compute the gradient (s)
    s <- lambda_hat * (colSums(centered_kernel_mat_at_sampled) -
                         n * prob_grid_x %*% t(centered_kernel_mat_at_grid)) -
      tau_hat * weight_hat_vec / prob_sampled_x

    # Compute the inverse of the Hessian matrix
    if(type_of_q_is_prob == TRUE){
      Hessian_inv <- solve(lambda_hat^2 * n * (centered_kernel_mat_at_grid %*%
                             diag(prob_grid_x) %*% t(centered_kernel_mat_at_grid) -
                             (centered_kernel_mat_at_grid %*% prob_grid_x) %*%
                             t(centered_kernel_mat_at_grid %*% prob_grid_x)  ) +
                             diag(tau_hat / prob_sampled_x))
    }else{
      Hessian_inv <- solve(lambda_hat^2 * n * centered_kernel_mat_at_grid %*%
                             diag(prob_grid_x) %*% t(centered_kernel_mat_at_grid) +
                             diag(tau_hat / prob_sampled_x))
    }

    #print(s)
    #print(Hessian_inv)
    # Update the weight vector using the Newton-Raphson method
    #weight_hat_change <- NRstepsize * s %*% Hessian_inv
    weight_hat_vec <- weight_hat_vec + NRstepsize * s %*% Hessian_inv

    # Print progress every 10% of the iterations or at the first iteration
    if (i %% round(max_iteration / 10) == 0 || i == 1) {
      print(paste("Iteration", i, ": ||s||_2 =", format_number(norm(s, p = 2))))
    }
    #counter = counter + 1
  }

  return(weight_hat_vec)
}


#dim(centered_kernel_mat_at_grid %*% diag(prob_grid_x) %*%  t(centered_kernel_mat_at_grid))
#dim(diag(prob_grid_x))
#dim(centered_kernel_mat_at_grid)
# Example usage
# Define the weights for the mixture distribution
#mixture_weights <- c(1/2, 1/6, 1/6, 1/6)

# Define the parameters for the normal distributions
# First distribution: N(0, 1)
#means <- c(0, -1, 0, 1)
#sds <- c(1, 0.1, 0.1, 0.1)


#sampled_x <- sort(normal_mixture(1000, means, sds, mixture_weights))
#x_grid <-  seq(-4,4,length.out =2000)
#centering_grid <- sort(runif(1000, min = -4, max =4))

#centered_kernel_mat_at_sampled <- centered_kernel_matrix(first_vec_kernel = sampled_x,
#                                                         second_vec_kernel = sampled_x,
#                                                         centering_grid = centering_grid,
#                                                         hurst_coef = 0.5)
#centered_kernel_mat_at_grid <- centered_kernel_matrix(first_vec_kernel = sampled_x,
#                                                         second_vec_kernel = x_grid,
#                                                         centering_grid = centering_grid,
#                                                         hurst_coef = 0.5)
#centered_kernel_self_grid <- diag(centered_kernel_matrix(first_vec_kernel = x_grid,
#                                                        second_vec_kernel = x_grid,
#                                                        centering_grid = centering_grid,
#                                                        hurst_coef = 0.5))




# Save the entire global environment to a file
#save.image(file = "my_environment.RData")
#lambda_hat <- 10
#tau_hat <- 10

#weights_hat <- get_weights(lambda_hat =lambda_hat, tau_hat = tau_hat,
#                      centered_kernel_mat_at_sampled, centered_kernel_mat_at_grid,
#                      centered_kernel_self_grid, x_grid = x_grid)

#plot(sampled_x,weights_hat[1,])


#probs <- get_probs(centered_kernel_mat_at_sampled, centered_kernel_mat_at_grid,
#                   centered_kernel_self_grid, x_grid, lambda_hat, weights_hat)





# Plot the histogram of the samples
#hist(sampled_x, breaks = 50, probability = TRUE, main = paste("Density estimation", "lambda =", lambda_hat, "tau =", tau_hat), xlab = "Value")

# Add density lines for each of the component distributions
#curve((mixture_weights[1] * dnorm(x, mean = means[1], sd = sds[1]) +
#         mixture_weights[2] * dnorm(x, mean = means[2], sd = sds[2]) +
#         mixture_weights[3] * dnorm(x, mean = means[3], sd = sds[3]) +
#         mixture_weights[4] * dnorm(x, mean = means[4], sd = sds[4])) , add = TRUE, col = "red",lwd = 2)


# Perform KDE with adaptive bandwidth
#kde_result <- kde(sampled_x)

# Overlay the KDE plot on the histogram
#lines(kde_result$eval.points, kde_result$estimate, col = "blue", lwd  =2)

#lines(x_grid, probs$grid_x, col = "orange", lwd  =2)

