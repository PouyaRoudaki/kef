#' Estimate Weights Using the Barzilai-Borwein Method (Optimized with Rcpp)
#'
#' This function estimates the weight vector using a Barzilai-Borwein method.
#'
#' @param lambda_hat A scalar representing the estimated lambda parameter.
#' @param tau_hat A scalar representing the estimated tau parameter.
#' @param centered_kernel_mat_at_sampled A matrix representing the centered kernel at sampled points.
#' @param sampled_x A vector of sampled points.
#' @param min_x The minimum domain value.
#' @param max_x The maximum domain value.
#' @param print_trace Logical; if TRUE, prints progress updates.
#'
#' @return A numeric vector of estimated weights.
#' @export
get_weights_wo_grid_BBsolve <- function(lambda_hat, tau_hat, centered_kernel_mat_at_sampled,
                                        sampled_x, min_x, max_x, print_trace = FALSE) {

  n <- nrow(centered_kernel_mat_at_sampled)  # Number of sampled points

  weight_hat_vec <- rep(0,n)
  min_x <- -3.1
  max_x <- 3.1
  # Wrapper for BBsolve using the Rcpp function
  s_function <- function(weight_hat_vec) {
    result <- get_s_function(weight_hat_vec, lambda_hat, tau_hat,
                             centered_kernel_mat_at_sampled, sampled_x, min_x, max_x)

    return(as.numeric(result))  # Ensure it's a standard numeric vector
  }

  # Solve using BBsolve (still in R)
  result <- BBsolve(par = rep(0, n), fn = s_function, control = list(maxit = 10000,
                                                                     tol = 1e-4,
                                                                     trace = print_trace))

  if (print_trace) {
    print(result$message)
  }

  return(result$par)  # Return optimized weight vector
}
