#' Kernel Density Estimation Using the KEF Method
#'
#' This function estimates density using kernel exponential family (KEF) methods.
#' It constructs a kernel matrix and computes weights using the Barzilai-Borwein
#' optimization method. If the sample size is large, it avoids computational complexity
#' by returning density estimates only at sample points.
#'
#' @param sample A numeric vector representing the observed sample points.
#' @param grid A numeric vector representing the evaluation grid.
#' @param lambda A numeric scalar for the lambda parameter (regularization term).
#' @param tau A numeric scalar for the tau parameter.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{weights} - The estimated weight vector.
#'   \item \code{probs_sample} - Density estimates at sample points.
#'   \item \code{probs_grid} - Density estimates at grid points (only if \code{length(sample) < 1000}).
#'   \item \code{time} - The function's execution time.
#' }
#'
#' @details
#' If the sample size is \code{>= 1000}, the function **only returns densities at sample points**
#' to avoid computational complexity. A warning is issued in this case.
#'
#' @export
#'
#' @examples
#' sample <- rnorm(500)
#' grid <- seq(-3, 3, length.out = 100)
#' lambda <- 0.1
#' tau <- 0.5
#' result <- kef(sample, grid, lambda, tau)
#' plot(grid, result$probs_grid, type = "l", main = "Estimated Density", ylab = "Density")

kef <- function(sample, grid, lambda, tau, with_base = TRUE, with_prob_in_var = TRUE) {

  # Start timer
  start_time <- Sys.time()

  # Compute the centered kernel matrix at sampled points
  centered_kernel_mat_at_sampled <- centered_kernel_matrix(
    first_vec_kernel = sample,
    second_vec_kernel = sample,
    centering_grid = grid,
    hurst_coef = 0.5
  )

  # Check if density should be computed on the grid
  density_only_sample <- length(sample) >= 1000
  if (density_only_sample) {
    warning("Your sample size is large. To avoid computational complexity,
             the output includes only density estimates at sample points.")
  }

  # Compute kernel matrices for the grid if needed
  if (!density_only_sample) {
    centered_kernel_mat_at_grid <- centered_kernel_matrix(
      first_vec_kernel = sample,
      second_vec_kernel = grid,
      centering_grid = grid,
      hurst_coef = 0.5
    )

    centered_kernel_self_grid <- diag(centered_kernel_matrix(
      first_vec_kernel = grid,
      second_vec_kernel = grid,
      centering_grid = grid,
      hurst_coef = 0.5
    ))
  }

  # Estimate the weight vector using the Barzilai-Borwein optimization method
  weights_hat_wo_grid <- get_weights_wo_grid_BBsolve(
    lambda_hat = lambda,
    tau_hat = tau,
    centered_kernel_mat_at_sampled = centered_kernel_mat_at_sampled,
    sampled_x = sample,
    min_x = min(grid),
    max_x = max(grid),
    with_base = with_base,
    with_prob_in_var = with_prob_in_var
  )


  # Compute density estimates based on whether grid evaluation is required
  if (!density_only_sample) {
    probs <- get_dens_or_prob(
      centered_kernel_mat_at_sampled,
      centered_kernel_mat_at_grid,
      centered_kernel_self_grid,
      sample,
      grid,
      lambda_hat = lambda,
      as.vector(weights_hat_wo_grid),
      type_of_p_is_prob = FALSE,
      type_of_q_is_prob = FALSE
    )


    # Store results including grid estimates
    result_list <- list(
      weights = as.vector(weights_hat_wo_grid),
      probs_sample = as.vector(probs$sampled_x),
      probs_grid = as.vector(probs$grid_x)
    )
  } else {
    probs <- get_dens_wo_grid(
      centered_kernel_mat_at_sampled,
      min_x = min(grid),
      max_x = max(grid),
      sampled_x = sample,
      lambda_hat = lambda,
      weight_hat_vec = as.vector(weights_hat_wo_grid)
    )


    # Store results (no grid estimates)
    result_list <- list(
      weights = as.vector(weights_hat_wo_grid),
      probs_sample = as.vector(probs)
    )
  }

  # End timer
  end_time <- Sys.time()

  # Compute total execution time
  total_time <- difftime(end_time, start_time, units = "secs")

  # Store the execution time.
  result_list$time <- total_time

  #return(invisible(result_list))  # Avoids unnecessary console output
  return(result_list)
}

