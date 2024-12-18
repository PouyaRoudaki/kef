
#' Compute LOOCV Errors on a Grid (Parallelized)
#'
#' This function calculates LOOCV weight errors for a given grid of lambda and tau values using parallelization.
#'
#' @param centered_kernel_mat_at_sampled A square matrix (n x n) representing the centered
#'        kernel matrix evaluated at the sampled points, where n is the number of sampled points.
#' @param sampled_x A vector of sampled points for which the weights are to be estimated.
#' @param min_x A scalar representing the minimum value of the domain.
#' @param max_x A scalar representing the maximum value of the domain.
#' @param lambda_hat_grid A numeric vector of lambda values for the grid.
#' @param tau_hat_grid A numeric vector of tau values for the grid.
#' @param cloud_computing Logical; if `TRUE`, all available cores are used for parallelization. Default is `FALSE`.
#'
#' @return A data frame containing:
#'   - `lambda_hat`: The lambda values from the grid.
#'   - `tau_hat`: The tau values from the grid.
#'   - `loocv_err_mean`: The computed jackknife errors mean.
#'   - `loocv_err_se`: The computed jackknife errors standard error.
#' @export
#'
#' @examples
#' # Example usage:
#' lambda_grid <- seq(0.1, 1, by = 0.1)
#' tau_grid <- seq(0.1, 1, by = 0.1)
#' sampled_data <- rnorm(100)
#' grid_points <- seq(-3, 3, length.out = 50)
#'
#' # Assume kernel matrices are precomputed:
#' centered_kernel_sampled <- matrix(runif(10000), ncol = 100)
#' centered_kernel_grid <- matrix(runif(5000), ncol = 50)
#' self_centered_grid <- matrix(runif(2500), ncol = 50)
#'
#' result <- loocv_error_grid_inner_parallelized(
#'   centered_kernel_mat_at_sampled = centered_kernel_sampled,
#'   sampled_x = sampled_data,
#'   min_x = min(grid_points),
#'   max_x = max(grid_points),
#'   lambda_hat_grid = lambda_grid,
#'   tau_hat_grid = tau_grid,
#'   cloud_computing = FALSE
#' )
loocv_error_grid_inner_parallelized <- function(centered_kernel_mat_at_sampled,
                                                           sampled_x,
                                                           min_x,
                                                           max_x,
                                                           lambda_hat_grid,
                                                           tau_hat_grid,
                                                           cloud_computing = FALSE) {
  # Create a grid of lambda_hat and tau_hat values
  grid <- expand.grid(lambda_hat = lambda_hat_grid, tau_hat = tau_hat_grid)

  # Initialize results
  err_mean <- numeric(nrow(grid))
  err_se <- numeric(nrow(grid))

  # Number of cores for inner parallelization
  num_cores_inner <- if (cloud_computing) parallel::detectCores() else max(1, parallel::detectCores() - 2)

  for (outer_index in seq_len(nrow(grid))) {

    lambda_hat <- grid$lambda_hat[outer_index]
    tau_hat <- grid$tau_hat[outer_index]
    print(paste0("outer_index = ",outer_index, " lambda = ",lambda_hat, " tau = ", tau_hat))

    # Create a cluster for the inner loop
    cl_inner <- parallel::makeCluster(num_cores_inner)

    # Export necessary variables to the cluster
    parallel::clusterExport(cl_inner,
                            varlist = c("get_weights_wo_grid",
                                        "centered_kernel_mat_at_sampled",
                                        "sampled_x", "min_x","max_x",
                                        "lambda_hat", "tau_hat"),
                            envir = environment())

    # Source or reload the updated function in each worker
    parallel::clusterEvalQ(cl_inner, {
      library(kef)
      library(pracma)
    })

    # Perform leave-one-out computations in parallel with tryCatch
    loocv_err_vec <- unlist(parallel::parLapply(cl_inner, seq_len(nrow(centered_kernel_mat_at_sampled)), function(i) {
      tryCatch({
        temp_centered_kernel_mat_at_sampled <- centered_kernel_mat_at_sampled[-i, -i]
        temp_sampled_x <- sampled_x[-i]

        w_wo_i <- as.numeric(get_weights_wo_grid(lambda_hat = lambda_hat,
                                                 tau_hat = tau_hat,
                                                 temp_centered_kernel_mat_at_sampled,
                                                 temp_sampled_x,
                                                 min_x,
                                                 max_x))

        dens_wo_i <- get_dens_wo_grid(temp_centered_kernel_mat_at_sampled,
                                      min_x,
                                      max_x,
                                      temp_sampled_x,
                                      lambda_hat,
                                      w_wo_i)

        prob_wo_i <- dens_wo_i / sum(dens_wo_i)

        one_out_err <- centered_kernel_mat_at_sampled[i, i] -
          2 * lambda_hat * prob_wo_i %*% centered_kernel_mat_at_sampled[-i, i] +
          lambda_hat^2 * prob_wo_i %*% centered_kernel_mat_at_sampled[-i, -i] %*% prob_wo_i

        return(one_out_err)
      }, error = function(e) {
        message(sprintf("Non-invertible Hessian for iteration = %d, lambda_hat = %f, tau_hat = %f: %s", i, lambda_hat, tau_hat, e$message))
        return(NA) # Return NA if an error occurs
      })
    }))


    # Store the error for the current outer index
    err_mean[outer_index] <- mean(loocv_err_vec)

    # Store the error for the current outer index
    err_se[outer_index] <- sqrt(var(loocv_err_vec))


    # Stop the cluster
    parallel::stopCluster(cl_inner)
  }

  # Combine the results into a data frame
  results <- data.frame(lambda_hat = grid$lambda_hat, tau_hat = grid$tau_hat,
                        loocv_err_mean = err_mean,
                        loocv_err_se = err_se)
  return(results)
}
