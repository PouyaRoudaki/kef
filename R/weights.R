#' Estimate Weights Using the Barzilai-Borwein Method (Optimized with Rcpp)
#'
#' This function estimates the weight vector using the Barzilai-Borwein method,
#' a numerical approach to solve nonlinear systems of equations. It optimizes
#' performance by leveraging Rcpp for computation and adaptive subsampling for
#' better initial values in the Barzilai-Borwein method.
#'
#' @param lambda_hat A scalar representing the estimated lambda parameter.
#' @param tau_hat A scalar representing the estimated tau parameter.
#' @param centered_kernel_mat_at_sampled A matrix representing the centered kernel at sampled points.
#' @param sampled_x A numeric vector of sampled points.
#' @param min_x The minimum domain value.
#' @param max_x The maximum domain value.
#' @param prior_variance_p_vector A numeric vector for prior variance probabilities. Default is NULL.
#' @param print_trace Logical; if TRUE, prints progress updates.
#' @param init Logical; if TRUE, the weights calculation has a subsampling method for initialization.
#'
#' @return A numeric vector of estimated weights.
#' @export
get_weights_wo_grid_BBsolve <- function(lambda_hat, tau_hat, centered_kernel_mat_at_sampled,
                                        sampled_x, min_x, max_x,
                                        prior_variance_p_vector = NULL,
                                        print_trace = FALSE,
                                        init = FALSE) {

  n <- nrow(centered_kernel_mat_at_sampled)  # Number of sampled points

  # Wrapper for BBsolve using the Rcpp function
  s_function <- function(weight_hat_vec) {
    result <- get_s_function(weight_hat_vec, lambda_hat, tau_hat,
                             centered_kernel_mat_at_sampled, sampled_x, min_x, max_x,
                             prior_variance_p_vector)  # Pass prior_variance_p_vector

    return(as.numeric(result))  # Ensure it's a standard numeric vector
  }

  # Default initial weights for BBsolve.
  initial_weights <- rep(0, n)

  # If the sample size is larger than 500, we randomly select 100 samples
  # without replacement and estimate the weights for this subsample.
  # Using these estimated weights, we interpolate to obtain weights
  # for the entire dataset using fast linear interpolation.
  # Finally, we scale the interpolated weights according to the sample
  # size proportion by multiplying them by 100 / n.
  if(n > 500 && init == TRUE){

    # Take 100 random samples without replacement
    sampled_idx_100 <- sort(sample.int(n = n, size = 100, replace = FALSE))
    sampled_x_100 <- sampled_x[sampled_idx_100]
    centered_kernel_mat_at_sampled_100 <- centered_kernel_mat_at_sampled[sampled_idx_100, sampled_idx_100]

    # BBsolve Wrapper for subsampled data
    s_function_100 <- function(weight_hat_vec_100) {
      result <- get_s_function(weight_hat_vec_100, lambda_hat, tau_hat,
                               centered_kernel_mat_at_sampled_100, sampled_x_100, min_x, max_x,
                               prior_variance_p_vector)  # Pass prior_variance_p_vector

      return(as.numeric(result))  # Ensure it's a standard numeric vector
    }

    # Solve using BBsolve for subsample
    result_100 <- BB::BBsolve(par = rep(0,length(sampled_x_100)),
                          fn = s_function_100,
                          control = list(maxit = 10000,
                                         tol = 1e-6,
                                         trace = T))

    weights_100 <- result_100$par

    plot(sampled_x_100, weights_100)

    # Interpolate from 100 points to full dataset
    weights_full_data <- interp_linear(x = sampled_x_100, y = weights_100, xnew = sampled_x)

    # Scale weights based on sample size proportion
    initial_weights <- weights_full_data * 10 / sqrt(n)

    plot(sampled_x, initial_weights)
  }

  # Solve using BBsolve for full dataset
  result <- BB::BBsolve(par = as.numeric(initial_weights),
                    fn = s_function,
                    control = list(maxit = 10000,
                                   tol = 1e-4,
                                   trace = print_trace))
  if (print_trace) {
    print(result$message)
  }

  return(result$par)  # Return optimized weight vector
}
