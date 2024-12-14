  #' Calculate Jackknife Error for a Given Pair of Parameters
  #'
  #' Computes the jackknife error
  #' given specific values of \code{lambda_hat} and \code{tau_hat}.
  #'
  #' @param lambda_hat Numeric. Regularization parameter for weight estimation.
  #' @param tau_hat Numeric. Regularization parameter for controlling smoothness.
  #' @param centered_kernel_mat_at_sampled A square matrix (n x n) representing the centered
  #'        kernel matrix evaluated at the sampled points, where n is the number of sampled points.
  #' @param sampled_x A vector of sampled points for which the weights are to be estimated.
  #' @param min_x A scalar representing the minimum value of the domain.
  #' @param max_x A scalar representing the maximum value of the domain.
  #'
  #' @return Numeric. The calculated jackknife error for the given parameter pair.
  #' @export
  #'
  #' @examples
  #' # Example usage:
  #' # calc_jackknife_error(0.1, 0.5, weight_hat, ...)
  calc_jackknife_error <- function(lambda_hat,
                                   tau_hat,
                                   centered_kernel_mat_at_sampled,
                                   sampled_x,
                                   min_x,
                                   max_x) {

    weights_hat <- get_weights_wo_grid(lambda_hat,
                               tau_hat,
                               centered_kernel_mat_at_sampled,
                               sampled_x,
                               min_x,
                               max_x)

    jackknife_err <- sum(sapply(1:nrow(centered_kernel_mat_at_sampled), function(i) {
      # Remove the i-th sample
      temp_centered_kernel_mat_at_sampled <- centered_kernel_mat_at_sampled[-i, -i]
      #temp_centered_kernel_mat_at_grid <- centered_kernel_mat_at_grid[-i, ]
      temp_sampled_x <- sampled_x[-i]

      # Get jackknife estimated weights
      w_hat_jackknife <- as.numeric(get_weights(lambda_hat = lambda_hat,
                                                tau_hat = tau_hat,
                                                temp_centered_kernel_mat_at_sampled,
                                                sampled_x = temp_sampled_x,
                                                min_x,
                                                max_x))

      # Compute jackknife error for the i-th sample
      one_out_err <- t(w_hat_jackknife - weights_hat[-i]) %*%
        temp_centered_kernel_mat_at_sampled %*%
        (w_hat_jackknife - weights_hat[-i])

      return(one_out_err)
    }))

    return(jackknife_err)
  }



  #' Jackknife Weight Error Grid
  #'
  #' Computes jackknife error for a grid of lambda_hat and tau_hat values.
  #'
  #' @param centered_kernel_mat_at_sampled A square matrix (n x n) representing the centered
  #'        kernel matrix evaluated at the sampled points, where n is the number of sampled points.
  #' @param sampled_x A vector of sampled points for which the weights are to be estimated.
  #' @param min_x A scalar representing the minimum value of the domain.
  #' @param max_x A scalar representing the maximum value of the domain.
  #' @param lambda_hat_grid A v
  #' @param tau_hat_grid A vector of tau hat values.
  #'
  #' @return A data frame containing lambda_hat, tau_hat, and the corresponding Jackknife error (Jackknife_err).
  #' @export
  #'
  #' @examples
  #' # Example usage:
  #' # jackknife_weight_error_grid(...)
  jackknife_weight_error_grid <- function(centered_kernel_mat_at_sampled,
                                          sampled_x,
                                          min_x,
                                          max_x,
                                          lambda_hat_grid,
                                          tau_hat_grid) {
    # Create a grid of lambda_hat and tau_hat values
    grid <- expand.grid(lambda_hat = lambda_hat_grid, tau_hat = tau_hat_grid)

    # Create a cluster
    num_cores <- parallel::detectCores() - 1
    cl <- parallel::makeCluster(num_cores)


    #func_code <- parallel::clusterEvalQ(cl, deparse(get_weights))
    #print(func_code)
    # Export necessary objects and functions to the cluster
    parallel::clusterExport(cl,
                            varlist = c("get_weights_wo_grid", "calc_jackknife_error",
                                        "centered_kernel_mat_at_sampled",
                                        "sampled_x",
                                        "min_x",
                                        "max_x"),
                            envir = environment())

    # Source or reload the updated function in each worker to avoid stale versions
    parallel::clusterEvalQ(cl, devtools::load_all())

    # Compute jackknife error for each combination in the grid using parallelization
    err <- unlist(parallel::parLapply(cl, seq_len(nrow(grid)), function(idx) {
      lambda_hat <- grid$lambda_hat[idx]
      tau_hat <- grid$tau_hat[idx]

      calc_jackknife_error(lambda_hat,
                           tau_hat,
                           centered_kernel_mat_at_sampled,
                           sampled_x,
                           min_x,
                           max_x)
    }))

    # Stop the cluster
    parallel::stopCluster(cl)

    # Combine the results into a data frame
    results <- data.frame(lambda_hat = grid$lambda_hat, tau_hat = grid$tau_hat,
                          jackknife_err = err)
    return(results)
  }


#' Parallelized Calculation of Jackknife Error
#'
#' This function computes the jackknife error for a given set of inputs, with parallelization
#' applied to the leave-one-out computations to improve performance.
#'
#' @param lambda_hat Estimated lambda parameter.
#' @param tau_hat Estimated tau parameter.
#' @param centered_kernel_mat_at_sampled A matrix of kernel values for the sampled points, centered.
#' @param centered_kernel_mat_at_grid A matrix of kernel values between sampled points and grid points, centered.
#' @param centered_kernel_self_grid A vector of self-kernel values for grid points.
#' @param sampled_x A vector of sampled input values.
#' @param x_grid A vector of grid points where calculations are performed.
#' @param type_of_p_is_prob Logical, whether the "p" kernel is interpreted as a probability.
#' @param type_of_q_is_prob Logical, whether the "q" kernel is interpreted as a probability.
#' @param method_of_p_calculation A string indicating the method to compute "p".
#'
#' @return The computed jackknife error as a numeric value.
#' @export
#'
#' @examples
#' # Example usage of calc_jackknife_error_nested_parallel
#' jackknife_error <- calc_jackknife_error_nested_parallel(lambda_hat, tau_hat, kernel_sampled, kernel_grid,
#'                                                  kernel_self_grid, sampled_x, x_grid, TRUE, TRUE, "method")
calc_jackknife_error_nested_parallel <- function(lambda_hat,
                                          tau_hat,
                                          centered_kernel_mat_at_sampled,
                                          centered_kernel_mat_at_grid,
                                          centered_kernel_self_grid,
                                          sampled_x,
                                          x_grid,
                                          type_of_p_is_prob,
                                          type_of_q_is_prob,
                                          method_of_p_calculation,
                                          outer_index) {

  # Step 1: Compute the weights for the entire dataset
  weights_hat <- get_weights(lambda_hat = lambda_hat,
                             tau_hat = tau_hat,
                             centered_kernel_mat_at_sampled,
                             centered_kernel_mat_at_grid,
                             centered_kernel_self_grid,
                             sampled_x = sampled_x,
                             x_grid = x_grid,
                             type_of_p_is_prob,
                             type_of_q_is_prob,
                             method_of_p_calculation)

  # Step 2: Setup parallelization for leave-one-out computations
  num_cores_inner <- max(1, parallel::detectCores() - 2)  # Leave some cores free for the system
  cl_inner <- parallel::makeCluster(num_cores_inner)  # Create a cluster of cores



  # Export necessary variables to the cluster
  parallel::clusterExport(cl_inner,
                          varlist = c("get_weights","calc_jackknife_error_nested_parallel",
                                      "lambda_hat", "tau_hat",
                                      "centered_kernel_mat_at_sampled",
                                      "centered_kernel_mat_at_grid",
                                      "centered_kernel_self_grid",
                                      "sampled_x", "x_grid", "weights_hat",
                                      "type_of_p_is_prob", "type_of_q_is_prob",
                                      "method_of_p_calculation"),
                          envir = environment())
  # Source or reload the updated function in each worker to avoid stale versions
  parallel::clusterEvalQ(cl_inner, {
    library(kef)
    library(pracma)
  })


  jackknife_err <- sum(unlist(parallel::parLapply(cl_inner, 1:nrow(centered_kernel_mat_at_sampled), function(i) {
    # Define a log file for this worker
    #log_file <- paste0("log_outer_", outer_index, "_inner_", i, ".txt")

    # Open a sink to write logs
    #sink(log_file, append = TRUE)

    # Log details
    #cat("Processing leave-one-out index:", i, "\n")
    #cat("Connected to outer index:", outer_index, "\n")
    #cat("Process ID:", Sys.getpid(), "\n")
    #cat("Using get_weights function:\n")
    #print(get_weights)  # Log the function definition

    # Step 4.1: Remove the i-th sample to create a leave-one-out dataset
    temp_centered_kernel_mat_at_sampled <- centered_kernel_mat_at_sampled[-i, -i]
    temp_centered_kernel_mat_at_grid <- centered_kernel_mat_at_grid[-i, ]
    temp_sampled_x <- sampled_x[-i]

    # Step 4.2: Recompute the weights for the leave-one-out dataset
    w_hat_jackknife <- as.numeric(get_weights(lambda_hat = lambda_hat,
                                              tau_hat = tau_hat,
                                              temp_centered_kernel_mat_at_sampled,
                                              temp_centered_kernel_mat_at_grid,
                                              centered_kernel_self_grid,
                                              sampled_x = temp_sampled_x,
                                              x_grid = x_grid,
                                              type_of_p_is_prob = type_of_p_is_prob,
                                              type_of_q_is_prob = type_of_q_is_prob,
                                              method_of_p_calculation = method_of_p_calculation
                                              ))

    # Step 4.3: Compute the leave-one-out error for the i-th sample
    one_out_err <- t(w_hat_jackknife - weights_hat[-i]) %*%
      temp_centered_kernel_mat_at_sampled %*%
      (w_hat_jackknife - weights_hat[-i])

    # Close the sink
    #sink()

    return(one_out_err)
  })))

  # Step 5: Stop the cluster after computations are complete
  parallel::stopCluster(cl_inner)

  # Step 6: Return the computed jackknife error
  return(jackknife_err)
}




#' Jackknife Weight Error Grid
#'
#' Computes jackknife error for a grid of lambda_hat and tau_hat values.
#'
#' @param centered_kernel_mat_at_sampled Centered kernel matrix at sampled points.
#' @param centered_kernel_mat_at_grid Centered kernel matrix at grid points.
#' @param centered_kernel_self_grid Centered kernel matrix for the grid points themselves.
#' @param sampled_x Sampled x values.
#' @param x_grid Grid of x values where the function is evaluated.
#' @param lambda_hat_grid A vector of lambda hat values.
#' @param tau_hat_grid A vector of tau hat values.
#' @param type_of_p_is_prob Logical. If TRUE, p is treated as probabilities.
#' @param type_of_q_is_prob Logical. If TRUE, q is treated as probabilities.
#' @param method_of_p_calculation Character. Method for p calculation.
#'
#' @return A data frame containing lambda_hat, tau_hat, and the corresponding Jackknife error (Jackknife_err).
#' @export
#'
#' @examples
#' # Example usage:
#' # jackknife_weight_error_grid(...)
jackknife_weight_error_grid_nested_parallel <- function(centered_kernel_mat_at_sampled,
                                        centered_kernel_mat_at_grid,
                                        centered_kernel_self_grid,
                                        sampled_x,
                                        x_grid,
                                        lambda_hat_grid,
                                        tau_hat_grid,
                                        type_of_p_is_prob = TRUE,
                                        type_of_q_is_prob = TRUE,
                                        method_of_p_calculation = "ordinary",
                                        cloud_computing = FALSE,
                                        os = "Windows") {
  # Create a grid of lambda_hat and tau_hat values
  grid <- expand.grid(lambda_hat = lambda_hat_grid, tau_hat = tau_hat_grid)

  # Create a cluster
  if(cloud_computing == FALSE){
    num_cores <- max(1, parallel::detectCores() - 2)
  }else{
    num_cores <- parallel::detectCores()
  }


  if(os == "Windows"){
    cl_outer <- parallel::makeCluster(num_cores)


    #func_code <- parallel::clusterEvalQ(cl, deparse(get_weights))
    #print(func_code)
    # Export necessary objects and functions to the cluster
    parallel::clusterExport(cl_outer,
                            varlist = c("get_weights","calc_jackknife_error_nested_parallel",
                                        "centered_kernel_mat_at_sampled",
                                        "centered_kernel_mat_at_grid",
                                        "centered_kernel_self_grid",
                                        "sampled_x", "x_grid",
                                        "type_of_p_is_prob", "type_of_q_is_prob",
                                        "method_of_p_calculation"),
                            envir = environment())

    # Source or reload the updated function in each worker to avoid stale versions
    parallel::clusterEvalQ(cl_outer, {
      library(kef)
      library(pracma)
    })

  # Compute jackknife error for each combination in the grid using parallelization
  err <- unlist(parallel::parLapply(cl_outer, seq_len(nrow(grid)), function(outer_index) {
    # Log the processing index and process ID
    #log_file <- paste0("log_outer_", outer_index, ".txt")
    #sink(log_file, append = TRUE)
    #cat("Starting computation for outer index:", outer_index, "\n")
    #cat("Process ID:", Sys.getpid(), "\n\n")
    #sink()

    lambda_hat <- grid$lambda_hat[outer_index]
    tau_hat <- grid$tau_hat[outer_index]

    calc_jackknife_error_nested_parallel(
      lambda_hat = lambda_hat,
      tau_hat = tau_hat,
      centered_kernel_mat_at_sampled = centered_kernel_mat_at_sampled,
      centered_kernel_mat_at_grid = centered_kernel_mat_at_grid,
      centered_kernel_self_grid = centered_kernel_self_grid,
      sampled_x = sampled_x,
      x_grid = x_grid,
      type_of_p_is_prob = type_of_p_is_prob,
      type_of_q_is_prob = type_of_q_is_prob,
      method_of_p_calculation = method_of_p_calculation,
      outer_index = outer_index  # Pass the outer thread index
    )

    # Stop the cluster
    parallel::stopCluster(cl_outer)
  }))
  }else{
    # Define a function to compute jackknife error for a single index
    compute_error <- function(outer_index) {

      log_file <- paste0("log_process_", Sys.getpid(), ".txt")
      cat("Processing index:", outer_index, "on process:", Sys.getpid(), "\n", file = log_file, append = TRUE)

      lambda_hat <- grid$lambda_hat[outer_index]
      tau_hat <- grid$tau_hat[outer_index]

      calc_jackknife_error_nested_parallel(lambda_hat, tau_hat,
                           centered_kernel_mat_at_sampled,
                           centered_kernel_mat_at_grid,
                           centered_kernel_self_grid,
                           sampled_x, x_grid,
                           type_of_p_is_prob,
                           type_of_q_is_prob,
                           method_of_p_calculation,
                           outer_index)
    }

    # Use parallel::mclapply to compute jackknife errors
    err <- unlist(parallel::mclapply(seq_len(nrow(grid)), compute_error, mc.cores = num_cores))
  }

  # Combine the results into a data frame
  results <- data.frame(lambda_hat = grid$lambda_hat, tau_hat = grid$tau_hat,
                        jackknife_err = err)
  return(results)
}


#' Compute Jackknife Weight Errors on a Grid (Parallelized)
#'
#' This function calculates jackknife weight errors for a given grid of lambda and tau values using parallelization.
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
#'   - `jackknife_err`: The computed jackknife errors.
#'
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
#' result <- jackknife_weight_error_grid_inner_parallelized(
#'   centered_kernel_mat_at_sampled = centered_kernel_sampled,
#'   sampled_x = sampled_data,
#'   min_x = min(grid_points),
#'   max_x = max(grid_points),
#'   lambda_hat_grid = lambda_grid,
#'   tau_hat_grid = tau_grid,
#'   cloud_computing = FALSE
#' )
jackknife_weight_error_grid_inner_parallelized <- function(centered_kernel_mat_at_sampled,
                                                           sampled_x,
                                                           min_x,
                                                           max_x,
                                                           lambda_hat_grid,
                                                           tau_hat_grid,
                                                           cloud_computing = FALSE) {
  # Create a grid of lambda_hat and tau_hat values
  grid <- expand.grid(lambda_hat = lambda_hat_grid, tau_hat = tau_hat_grid)

  # Initialize results
  err <- numeric(nrow(grid))

  # Number of cores for inner parallelization
  num_cores_inner <- if (cloud_computing) parallel::detectCores() else max(1, parallel::detectCores() - 2)

  for (outer_index in seq_len(nrow(grid))) {
    print(paste0("outer_index = ",outer_index))
    lambda_hat <- grid$lambda_hat[outer_index]
    tau_hat <- grid$tau_hat[outer_index]

    # Compute the weights for the entire dataset
    weights_hat <- get_weights_wo_grid(lambda_hat = lambda_hat,
                               tau_hat = tau_hat,
                               centered_kernel_mat_at_sampled,
                               sampled_x = sampled_x,
                               min_x = min_x,
                               max_x = max_x)

    # Create a cluster for the inner loop
    cl_inner <- parallel::makeCluster(num_cores_inner)

    # Export necessary variables to the cluster
    parallel::clusterExport(cl_inner,
                            varlist = c("get_weights_wo_grid","weights_hat",
                                        "centered_kernel_mat_at_sampled",
                                        "sampled_x", "min_x","max_x",
                                        "lambda_hat", "tau_hat"),
                            envir = environment())

    # Source or reload the updated function in each worker
    parallel::clusterEvalQ(cl_inner, {
      library(kef)
      library(pracma)
    })

    # Perform leave-one-out computations in parallel
    jackknife_err <- sum(unlist(parallel::parLapply(cl_inner, seq_len(nrow(centered_kernel_mat_at_sampled)), function(i) {

      temp_centered_kernel_mat_at_sampled <- centered_kernel_mat_at_sampled[-i, -i]
      temp_sampled_x <- sampled_x[-i]

      w_hat_jackknife <- as.numeric(get_weights_wo_grid(lambda_hat = lambda_hat,
                                                tau_hat = tau_hat,
                                                temp_centered_kernel_mat_at_sampled,
                                                temp_sampled_x,
                                                min_x,
                                                max_x))

      one_out_err <- t(w_hat_jackknife - weights_hat[-i]) %*%
        temp_centered_kernel_mat_at_sampled %*%
        (w_hat_jackknife - weights_hat[-i])
      return(one_out_err)
    })))

    # Store the error for the current outer index
    err[outer_index] <- jackknife_err

    # Stop the cluster
    parallel::stopCluster(cl_inner)
  }

  # Combine the results into a data frame
  results <- data.frame(lambda_hat = grid$lambda_hat, tau_hat = grid$tau_hat,
                        jackknife_err = err)
  return(results)
}
