#' Estimate Weights Using the Newton-Raphson Method
#'
#' This function estimates the weight vector using an iterative Newton-Raphson method.
#' The method updates weights based on the provided kernel matrices, regularization
#' parameters, sampled points, and grid points.
#'
#' @param lambda_hat A scalar representing the estimated lambda parameter, which
#'        controls the contribution of the kernel matrices to the weight estimation.
#' @param tau_hat A scalar representing the estimated tau parameter, which
#'        determines the regularization strength applied to the weights.
#' @param centered_kernel_mat_at_sampled A square matrix (n x n) representing the centered
#'        kernel matrix evaluated at the sampled points, where n is the number of sampled points.
#' @param centered_kernel_mat_at_grid A matrix (n x m) representing the centered kernel
#'        matrix evaluated at the grid points, where n is the number of sampled points
#'        and m is the number of grid points.
#' @param centered_kernel_self_grid A vector of length m representing the diagonal of
#'        the centered kernel matrix evaluated at the grid points.
#' @param sampled_x A vector of sampled points for which the weights are to be estimated.
#' @param x_grid A vector of grid points that define the evaluation domain for the kernel.
#' @param type_of_p_is_prob Logical; if TRUE, treats the "p" function as a probability density.
#' @param type_of_q_is_prob Logical; if TRUE, treats the "q" function as a probability density.
#' @param method_of_p_calculation A string indicating the method used to calculate probabilities
#'        for the "p" function. Default is "ordinary".
#' @param print_trace Logical; if TRUE, prints progress updates during the Newton-Raphson iterations.
#'
#' @return A numeric vector of length n representing the estimated weight vector for the sampled points.
#'
#' @export
get_weights <- function(lambda_hat,
                        tau_hat,
                        centered_kernel_mat_at_sampled,
                        centered_kernel_mat_at_grid,
                        centered_kernel_self_grid,
                        sampled_x,
                        x_grid,
                        type_of_p_is_prob=TRUE,
                        type_of_q_is_prob=TRUE,
                        method_of_p_calculation="ordinary",
                        print_trace = FALSE) {

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
    probs <- get_dens_or_prob(centered_kernel_mat_at_sampled,
                              centered_kernel_mat_at_grid,
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
    if ((i %% round(max_iteration / 10) == 0 || i == 1) & print_trace == TRUE) {
      print(paste("Iteration", i, ": ||s||_2 =", pracma::Norm(s)))
    }
    #counter = counter + 1
  }

  return(weight_hat_vec)
}


