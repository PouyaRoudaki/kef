#' Compute Marginal Log Likelihood Using Rcpp
#'
#' This function computes the marginal log likelihood using a Monte Carlo approach
#' with optional parallel computing.
#'
#' @param centered_kernel_mat_at_sampled A matrix of centered kernel values at sampled points.
#' @param sampled_x A numeric vector of sampled points.
#' @param min_x The minimum x value.
#' @param max_x The maximum x value.
#' @param p_vec A probability vector (default: uniform distribution).
#' @param lambda A scalar for the lambda hyperparameter.
#' @param tau A scalar for the tau hyperparameter.
#' @param std_rnorm_matrix A matrix of standard normal random values for Monte Carlo sampling.
#' @param MC_iterations The number of Monte Carlo iterations.
#' @param parallel_computing Logical; if TRUE, enables parallel computing.
#'
#' @return The computed log of the marginal likelihood.
#' @export
marginal_log_likelihood <- function(centered_kernel_mat_at_sampled,
                                    sampled_x,
                                    min_x,
                                    max_x,
                                    p_vec = rep(1, nrow(centered_kernel_mat_at_sampled)),
                                    lambda,
                                    tau,
                                    std_rnorm_matrix,
                                    MC_iterations,
                                    parallel_computing = TRUE) {

  # Call the C++ function using `.Call()`
  .Call("_kef_marginal_log_likelihood",
        centered_kernel_mat_at_sampled,
        sampled_x,
        min_x,
        max_x,
        p_vec,
        lambda,
        tau,
        std_rnorm_matrix,
        MC_iterations,
        parallel_computing)
}


#' Compute Marginal Likelihood Over a Grid of Hyperparameters (Parallel)
#'
#' This function computes the marginal likelihood for each combination of lambda and tau
#' in the provided hyperparameter grid using parallel computing.
#'
#' @param centered_kernel_mat_at_sampled A matrix of centered kernel values at sampled points.
#' @param min_x The minimum x value.
#' @param max_x The maximum x value.
#' @param sampled_x A numeric vector of sampled points.
#' @param hyperparam_grid A dataframe containing pairs of lambda and tau values.
#' @param initial_lambda The initial lambda value.
#' @param initial_w The initial weight vector.
#' @param MC_iterations The number of Monte Carlo iterations.
#' @param max_iterations The maximum number of iterations.
#' @param parallel_computing Logical; if TRUE, enables parallel computing.
#'
#' @return A dataframe containing lambda, tau, and their corresponding marginal log likelihoods.
#' @export
compute_marginal_likelihood_grid_parallel <- function(centered_kernel_mat_at_sampled,
                                                      min_x,
                                                      max_x,
                                                      sampled_x,
                                                      hyperparam_grid,
                                                      initial_lambda = 1,
                                                      initial_w = rep(0, length(sampled_x)),
                                                      MC_iterations,
                                                      max_iterations = 5,
                                                      parallel_computing = TRUE) {

  # Ensure hyperparam_grid is a matrix
  if (!is.matrix(hyperparam_grid) && !is.data.frame(hyperparam_grid)) {
    stop("Error: 'hyperparam_grid' must be a matrix or a dataframe.")
  }

  # Convert dataframe to matrix if needed
  if (is.data.frame(hyperparam_grid)) {
    hyperparam_grid <- as.matrix(hyperparam_grid)
  }

  # Call the C++ function using `.Call()`
  results <- .Call("_kef_compute_marginal_likelihood_grid_parallel",
                   centered_kernel_mat_at_sampled,
                   min_x,
                   max_x,
                   sampled_x,
                   hyperparam_grid,
                   initial_lambda,
                   initial_w,
                   MC_iterations,
                   max_iterations,
                   parallel_computing)

  # Convert the output matrix to a dataframe
  results_df <- as.data.frame(results)
  colnames(results_df) <- c("lambda", "tau", "marginal_log_likelihood")

  return(results_df)
}


#' Optimize Marginal Log-Likelihood with Convergence Check
#'
#' This function optimizes the hyperparameters \eqn{\lambda} and \eqn{\tau} by maximizing
#' the marginal log-likelihood using L-BFGS-B optimization. Instead of a grid search,
#' it efficiently finds the best parameters while checking for convergence.
#'
#' @param centered_kernel_mat_at_sampled The kernel matrix centered at sampled points.
#' @param min_x Minimum value of the sampled domain.
#' @param max_x Maximum value of the sampled domain.
#' @param sampled_x Vector of sampled points.
#' @param initial_lambda Initial value for \eqn{\lambda} (default: 1).
#' @param initial_w Initial weight vector (default: zeros of length `sampled_x`).
#' @param MC_iterations Number of Monte Carlo iterations.
#' @param max.iterations Maximum number of iterations (default: 50).
#' @param tol Convergence tolerance (default: `1e-4`). If the relative change in
#' \eqn{\lambda} and \eqn{\tau} falls below this threshold, optimization stops early.
#' @param parallel_computing Boolean flag indicating whether parallelization should be
#' used for optimization (default: `TRUE`).
#'
#' @return A list containing:
#' \describe{
#'   \item{lambda}{Optimized value of \eqn{\lambda}.}
#'   \item{tau}{Optimized value of \eqn{\tau}.}
#'   \item{max_marginal_log_likelihood}{Maximum marginal log-likelihood value.}
#'   \item{converged}{Logical value indicating whether the optimization converged
#'   before reaching `max.iterations`.}
#' }
#'
#' @details
#' - Uses the **L-BFGS-B** optimization method to efficiently find the best
#'   \eqn{\lambda} and \eqn{\tau} values.
#' - The optimization stops early if the relative change in \eqn{\lambda} and
#'   \eqn{\tau} is smaller than `tol`.
#' - If convergence is not achieved within `max.iterations`, the function
#'   prints a warning and returns the best estimate.
#'
#' @examples
#' \dontrun{
#' result <- optimize_marginal_log_likelihood(
#'   centered_kernel_mat_at_sampled, min_x, max_x, sampled_x,
#'   MC_iterations = 10000, max.iterations = 50, tol = 1e-4
#' )
#' print(result)
#' }
#'
#' @export
optimize_marginal_log_likelihood <- function(centered_kernel_mat_at_sampled,
                                               min_x,
                                               max_x,
                                               sampled_x,
                                               initial_lambda = 1,
                                               initial_w = rep(0, length(sampled_x)),
                                               MC_iterations = 10000,
                                               max.iterations = 10,
                                               tol = 1e-4,  # Convergence tolerance
                                               parallel_computing = TRUE) {
  converged <- FALSE
  t <- 1
  n <- length(sampled_x)

  # Generate matrix with each row independently sampled from Normal(0,1)
  set.seed(110)
  std_rnorm_matrix <- matrix(rnorm(MC_iterations * n, mean = 0, sd = 1),
                             nrow = MC_iterations,
                             ncol = n)

  lambda <- initial_lambda
  tau <- exp(log(lambda) - 4.5*log(10) + (.Machine$double.xmin))
  w_vec <- initial_w


  dens_vec <- get_dens_wo_grid(centered_kernel_mat_at_sampled,
                                          min_x,
                                          max_x,
                                          sampled_x,
                                          lambda,
                                          w_vec)

  p_vec <- dens_vec / sum(dens_vec)


  repeat {

    cat(paste0("Iteration: ", t, "\n"))

    # Store old values for convergence check
    lambda_old <- lambda
    tau_old <- tau

    # Define objective function for optimization
    objective_function <- function(params) {
      log_lambda <- params[1]  # Optimizing log(lambda)
      theta <- params[2]  # Free parameter for tau reparameterization
      log_tau <- log_lambda - 5*log(10) + exp(theta)  # Enforcing constraint 9.671
      lambda <- exp(log_lambda)
      tau <- exp(log_tau)



      - marginal_log_likelihood(
        centered_kernel_mat_at_sampled,
        sampled_x,
        min_x,
        max_x,
        p_vec,
        lambda,
        tau,
        std_rnorm_matrix,
        MC_iterations,
        parallel_computing)
    }

    cat(paste0("Initial lambda: ", lambda,", Initial tau: ", exp(log(lambda) - 5*log(10) + (.Machine$double.xmin)), "\n"))

    # Optimization using L-BFGS-B (bounded optimization)
    opt_result <- optim(
      par = c(log(lambda), -.Machine$double.xmax ),  # log(.Machine$double.xmin) Initial values for log(lambda) and theta
      fn = objective_function,
      method = "L-BFGS-B",
      lower = c(log(1e-1),  -Inf),
      upper = c(log(1e2), Inf)
    )

    # Retrieve optimal lambda and tau
    log_lambda <- opt_result$par[1]
    theta <- opt_result$par[2]
    log_tau <- log_lambda - 5*log(10) + exp(theta)
    lambda <- exp(log_lambda)
    tau <- exp(log_tau)

    cat(paste0("Optimized lambda: ", lambda, ", tau: ", tau, ", MLL: ", -opt_result$value,
               ", The ratio: ", lambda^2/tau ,"\n"))

    #cat(paste0("lambda_old: ",lambda_old,"tau_old: ",tau_old,".\n"))
    # Convergence check: Stop if parameters don't change significantly
    delta_lambda <- abs(lambda - lambda_old) / lambda_old
    delta_tau <- abs(tau - tau_old) / tau_old

    #cat(paste0("delta_lambda: ",delta_lambda,"delta_tau: ",delta_tau,".\n"))

    if (delta_lambda < tol & delta_tau < tol) {
      converged <- TRUE
      cat("✔ Convergence reached: Estimated Lambda and Tau by Maximum Marginal Likelihood are stable.\n")
      break
    }

    # Update weights
    w_vec <- get_weights_wo_grid_BBsolve(lambda_hat = lambda,
                                     tau_hat = tau,
                                     centered_kernel_mat_at_sampled = centered_kernel_mat_at_sampled,
                                     sampled_x = sampled_x,
                                     min_x = min_x,
                                     max_x = max_x,
                                     prior_variance_p_vector = p_vec,
                                     print_trace = FALSE)

    # Update density and p_vec
    dens_vec <- get_dens_wo_grid(centered_kernel_mat_at_sampled,
                                            min_x,
                                            max_x,
                                            sampled_x,
                                            lambda,
                                            w_vec)

    p_vec <- dens_vec / sum(dens_vec)

    t <- t + 1

    if (t > max.iterations) {
      cat("Error: Max iterations reached without full convergence.\n")
      break
    }
  }

  return(list(lambda = lambda, tau = tau, max_marginal_log_likelihood = -opt_result$value, converged = converged))
}


#' Optimize Marginal Log-Likelihood with Convergence Check
#'
#' This function optimizes the hyperparameters \eqn{\lambda} and \eqn{\tau} by maximizing
#' the marginal log-likelihood using L-BFGS-B optimization. Instead of a grid search,
#' it efficiently finds the best parameters while checking for convergence.
#'
#' @param centered_kernel_mat_at_sampled The kernel matrix centered at sampled points.
#' @param min_x Minimum value of the sampled domain.
#' @param max_x Maximum value of the sampled domain.
#' @param sampled_x Vector of sampled points.
#' @param initial_lambda Initial value for \eqn{\lambda} (default: 1).
#' @param initial_tau Initial value for \eqn{\tau} (default: 1).
#' @param initial_w Initial weight vector (default: zeros of length `sampled_x`).
#' @param MC_iterations Number of Monte Carlo iterations.
#' @param max.iterations Maximum number of iterations (default: 5).
#' @param tol Convergence tolerance (default: `1e-4`). If the relative change in
#' \eqn{\lambda} and \eqn{\tau} falls below this threshold, optimization stops early.
#' @param parallel_computing Boolean flag indicating whether parallelization should be
#' used for optimization (default: `TRUE`).
#' @param seed An integer that controls the randomness.
#'
#' @return A list containing:
#' \describe{
#'   \item{lambda}{Optimized value of \eqn{\lambda}.}
#'   \item{tau}{Optimized value of \eqn{\tau}.}
#'   \item{max_marginal_log_likelihood}{Maximum marginal log-likelihood value.}
#'   \item{converged}{Logical value indicating whether the optimization converged
#'   before reaching `max.iterations`.}
#' }
#'
#' @details
#' - Uses the **L-BFGS-B** optimization method to efficiently find the best
#'   \eqn{\lambda} and \eqn{\tau} values.
#' - The optimization stops early if the relative change in \eqn{\lambda} and
#'   \eqn{\tau} is smaller than `tol`.
#' - If convergence is not achieved within `max.iterations`, the function
#'   prints a warning and returns the best estimate.
#'
#' @examples
#' \dontrun{
#' result <- optimize_marginal_log_likelihood(
#'   centered_kernel_mat_at_sampled, min_x, max_x, sampled_x,
#'   MC_iterations = 10000, max.iterations = 50, tol = 1e-4
#' )
#' print(result)
#' }
#'
#' @export
optimize_marginal_log_likelihood_new <- function(centered_kernel_mat_at_sampled,
                                             min_x,
                                             max_x,
                                             sampled_x,
                                             initial_lambda = 1,
                                             initial_tau = 1,
                                             initial_w = rep(0, length(sampled_x)),
                                             MC_iterations = 10000,
                                             max.iterations = 5,
                                             tol = 1e-4,  # Convergence tolerance
                                             parallel_computing = TRUE,
                                             seed = 1) {

  converged <- FALSE
  t <- 1
  n <- length(sampled_x)

  #set.seed(13)
  # Generate matrix with each row independently sampled from Normal(0,1)
  set.seed(seed)
  std_rnorm_matrix <- matrix(rnorm(MC_iterations * n, mean = 0, sd = 1),
                             nrow = MC_iterations,
                             ncol = n)

  lambda <- initial_lambda
  tau <- initial_tau
  w_vec <- initial_w


  dens_vec <- get_dens_wo_grid(centered_kernel_mat_at_sampled,
                               min_x,
                               max_x,
                               sampled_x,
                               lambda,
                               w_vec)

  p_vec <- dens_vec / sum(dens_vec)

  p_vec_init <- as.numeric(p_vec)

  repeat {

    cat(paste0("Iteration: ", t, "\n"))

    # Store old values for convergence check
    lambda_old <- lambda
    tau_old <- tau

    # Define objective function for optimization
    objective_function <- function(params) {
      #log_lambda <- params[1]  # Optimizing log(lambda)
      #theta <- params[2]  # Free parameter for tau reparameterization
      #log_tau <- log_lambda - 4.2*log(10) + exp(theta)  # Enforcing constraint 9.671
      #log_tau <- params[2]
      #lambda <- exp(log_lambda)
      #tau <- exp(log_tau)
      lambda <- params[1]
      tau <- params[2]


      - marginal_log_likelihood(
        centered_kernel_mat_at_sampled,
        sampled_x,
        min_x,
        max_x,
        p_vec,
        lambda,
        tau,
        std_rnorm_matrix,
        MC_iterations,
        parallel_computing)
    }

    cat(paste0("Initial lambda: ", lambda,", Initial tau: ", tau, "\n"))

    # Optimization using L-BFGS-B (bounded optimization)
    opt_result <- optim(
      par = c(lambda, tau),  # log(.Machine$double.xmin) Initial values for log(lambda) and theta
      fn = objective_function,
      method = "L-BFGS-B",
      lower = c(1e-1, 1e-6),
      upper = c(1e2, 1e2)
    )

    # Retrieve optimal lambda and tau
    #log_lambda <- opt_result$par[1]
    #theta <- opt_result$par[2]
    #log_tau <- log_lambda - 4.2*log(10) + exp(theta)
    #log_tau <- opt_result$par[2]
    #lambda <- exp(log_lambda)
    #tau <- exp(log_tau)
    lambda <- opt_result$par[1]
    tau <- opt_result$par[2]


    cat(paste0("Optimized lambda: ", lambda, ", tau: ", tau, ", MLL: ", -opt_result$value,
               ", The ratio: ", lambda^2/tau ,"\n"))

    #cat(paste0("lambda_old: ",lambda_old,"tau_old: ",tau_old,".\n"))
    # Convergence check: Stop if parameters don't change significantly
    delta_lambda <- abs(lambda - lambda_old) / lambda_old
    delta_tau <- abs(tau - tau_old) / tau_old

    #cat(paste0("delta_lambda: ",delta_lambda,"delta_tau: ",delta_tau,".\n"))

    if (delta_lambda < tol & delta_tau < tol) {
      converged <- TRUE
      cat("✔ Convergence reached: Estimated Lambda and Tau by Maximum Marginal Likelihood are stable.\n")
      break
    }

    # Update weights
    w_vec <- get_weights_wo_grid_BBsolve(lambda_hat = lambda,
                                         tau_hat = tau,
                                         centered_kernel_mat_at_sampled = centered_kernel_mat_at_sampled,
                                         sampled_x = sampled_x,
                                         min_x = min_x,
                                         max_x = max_x,
                                         prior_variance_p_vector = p_vec,
                                         print_trace = FALSE)

    # Update density and p_vec
    dens_vec <- get_dens_wo_grid(centered_kernel_mat_at_sampled,
                                 min_x,
                                 max_x,
                                 sampled_x,
                                 lambda,
                                 w_vec)

    p_vec <- dens_vec / sum(dens_vec)

    t <- t + 1

    if (t > max.iterations) {
      cat("Error: Max iterations reached without full convergence.\n")
      break
    }
  }

  return(list(lambda = lambda,
              tau = tau,
              max_marginal_log_likelihood = -opt_result$value,
              converged = converged,
              std_rnorm_matrix = std_rnorm_matrix,
              p_vec = p_vec_init))
}



#' Optimize Marginal Log-Likelihood with Convergence Check
#'
#' This function optimizes the hyperparameters \eqn{\lambda} and \eqn{\tau} by maximizing
#' the marginal log-likelihood using L-BFGS-B optimization. Instead of a grid search,
#' it efficiently finds the best parameters while checking for convergence.
#'
#' @param centered_kernel_mat_at_sampled The kernel matrix centered at sampled points.
#' @param min_x Minimum value of the sampled domain.
#' @param max_x Maximum value of the sampled domain.
#' @param sampled_x Vector of sampled points.
#' @param initial_lambda Initial value for \eqn{\lambda} (default: 1).
#' @param initial_tau Initial value for \eqn{\tau} (default: 1).
#' @param initial_w Initial weight vector (default: zeros of length `sampled_x`).
#' @param MC_iterations Number of Monte Carlo iterations.
#' @param max.iterations Maximum number of iterations (default: 5).
#' @param tol Convergence tolerance (default: `1e-4`). If the relative change in
#' \eqn{\lambda} and \eqn{\tau} falls below this threshold, optimization stops early.
#' @param parallel_computing Boolean flag indicating whether parallelization should be
#' used for optimization (default: `TRUE`).
#' @param seed An integer that controls the randomness.
#'
#' @return A list containing:
#' \describe{
#'   \item{lambda}{Optimized value of \eqn{\lambda}.}
#'   \item{tau}{Optimized value of \eqn{\tau}.}
#'   \item{max_marginal_log_likelihood}{Maximum marginal log-likelihood value.}
#'   \item{converged}{Logical value indicating whether the optimization converged
#'   before reaching `max.iterations`.}
#' }
#'
#' @details
#' - Uses the **L-BFGS-B** optimization method to efficiently find the best
#'   \eqn{\lambda} and \eqn{\tau} values.
#' - The optimization stops early if the relative change in \eqn{\lambda} and
#'   \eqn{\tau} is smaller than `tol`.
#' - If convergence is not achieved within `max.iterations`, the function
#'   prints a warning and returns the best estimate.
#'
#' @examples
#' \dontrun{
#' result <- optimize_marginal_log_likelihood(
#'   centered_kernel_mat_at_sampled, min_x, max_x, sampled_x,
#'   MC_iterations = 10000, max.iterations = 50, tol = 1e-4
#' )
#' print(result)
#' }
#'
#' @export
optimize_marginal_log_likelihood_init_p <- function(centered_kernel_mat_at_sampled,
                                                 min_x,
                                                 max_x,
                                                 sampled_x,
                                                 initial_lambda = 1,
                                                 initial_tau = 1/1350,
                                                 MC_iterations = 10000,
                                                 max.iterations = 5,
                                                 tol = 1e-4,  # Convergence tolerance
                                                 parallel_computing = TRUE,
                                                 seed = 1) {

  converged <- FALSE
  t <- 1
  n <- length(sampled_x)



  #set.seed(13)
  # Generate matrix with each row independently sampled from Normal(0,1)
  set.seed(seed)
  std_rnorm_matrix <- matrix(rnorm(MC_iterations * n, mean = 0, sd = 1),
                             nrow = MC_iterations,
                             ncol = n)

  lambda <- initial_lambda
  tau <- initial_tau

  grid <- seq(from = min_x, to = max_x,length.out = 4*n)

  dens_vec <- kef(sampled_x,grid,lambda = lambda,tau = tau)$probs_sample

  p_vec <- dens_vec / sum(dens_vec)

  p_vec_init <- as.numeric(p_vec)

  repeat {

    cat(paste0("Iteration: ", t, "\n"))

    # Store old values for convergence check
    lambda_old <- lambda
    tau_old <- tau

    # Define objective function for optimization
    objective_function <- function(params) {
      #log_lambda <- params[1]  # Optimizing log(lambda)
      #theta <- params[2]  # Free parameter for tau reparameterization
      #log_tau <- log_lambda - 4.2*log(10) + exp(theta)  # Enforcing constraint 9.671
      #log_tau <- params[2]
      #lambda <- exp(log_lambda)
      #tau <- exp(log_tau)
      lambda <- params[1]
      tau <- params[2]


      - marginal_log_likelihood(
        centered_kernel_mat_at_sampled,
        sampled_x,
        min_x,
        max_x,
        p_vec,
        lambda,
        tau,
        std_rnorm_matrix,
        MC_iterations,
        parallel_computing)
    }

    cat(paste0("Initial lambda: ", lambda,", Initial tau: ", tau, "\n"))

    # Optimization using L-BFGS-B (bounded optimization)
    opt_result <- optim(
      par = c(lambda, tau),  # log(.Machine$double.xmin) Initial values for log(lambda) and theta
      fn = objective_function,
      method = "L-BFGS-B",
      lower = c(1e-1, 1e-6),
      upper = c(1e2, 1e2)
    )

    # Retrieve optimal lambda and tau
    #log_lambda <- opt_result$par[1]
    #theta <- opt_result$par[2]
    #log_tau <- log_lambda - 4.2*log(10) + exp(theta)
    #log_tau <- opt_result$par[2]
    #lambda <- exp(log_lambda)
    #tau <- exp(log_tau)
    lambda <- opt_result$par[1]
    tau <- opt_result$par[2]


    cat(paste0("Optimized lambda: ", lambda, ", tau: ", tau, ", MLL: ", -opt_result$value,
               ", The ratio: ", lambda^2/tau ,"\n"))

    #cat(paste0("lambda_old: ",lambda_old,"tau_old: ",tau_old,".\n"))
    # Convergence check: Stop if parameters don't change significantly
    delta_lambda <- abs(lambda - lambda_old) / lambda_old
    delta_tau <- abs(tau - tau_old) / tau_old

    #cat(paste0("delta_lambda: ",delta_lambda,"delta_tau: ",delta_tau,".\n"))

    if (delta_lambda < tol & delta_tau < tol) {
      converged <- TRUE
      cat("✔ Convergence reached: Estimated Lambda and Tau by Maximum Marginal Likelihood are stable.\n")
      break
    }

    # Update weights
    w_vec <- get_weights_wo_grid_BBsolve(lambda_hat = lambda,
                                         tau_hat = tau,
                                         centered_kernel_mat_at_sampled = centered_kernel_mat_at_sampled,
                                         sampled_x = sampled_x,
                                         min_x = min_x,
                                         max_x = max_x,
                                         prior_variance_p_vector = p_vec,
                                         print_trace = FALSE)

    # Update density and p_vec
    dens_vec <- get_dens_wo_grid(centered_kernel_mat_at_sampled,
                                 min_x,
                                 max_x,
                                 sampled_x,
                                 lambda,
                                 w_vec)

    p_vec <- dens_vec / sum(dens_vec)

    t <- t + 1

    if (t > max.iterations) {
      cat("Error: Max iterations reached without full convergence.\n")
      break
    }
  }

  return(list(lambda = lambda,
              tau = tau,
              max_marginal_log_likelihood = -opt_result$value,
              converged = converged,
              std_rnorm_matrix = std_rnorm_matrix,
              p_vec = p_vec_init))
}



#' Optimize Marginal Log-Likelihood with Convergence Check
#'
#' This function optimizes the hyperparameters \eqn{\lambda} and \eqn{\tau} by maximizing
#' the marginal log-likelihood using L-BFGS-B optimization. Instead of a grid search,
#' it efficiently finds the best parameters while checking for convergence.
#'
#' @param centered_kernel_mat_at_sampled The kernel matrix centered at sampled points.
#' @param min_x Minimum value of the sampled domain.
#' @param max_x Maximum value of the sampled domain.
#' @param sampled_x Vector of sampled points.
#' @param initial_lambda Initial value for \eqn{\lambda} (default: 1).
#' @param initial_tau Initial value for \eqn{\tau} (default: 1).
#' @param initial_w Initial weight vector (default: zeros of length `sampled_x`).
#' @param MC_iterations Number of Monte Carlo iterations.
#' @param max.iterations Maximum number of iterations (default: 50).
#' @param tol Convergence tolerance (default: `1e-4`). If the relative change in
#' \eqn{\lambda} and \eqn{\tau} falls below this threshold, optimization stops early.
#' @param parallel_computing Boolean flag indicating whether parallelization should be
#' used for optimization (default: `TRUE`).
#'
#' @return A list containing:
#' \describe{
#'   \item{lambda}{Optimized value of \eqn{\lambda}.}
#'   \item{tau}{Optimized value of \eqn{\tau}.}
#'   \item{max_marginal_log_likelihood}{Maximum marginal log-likelihood value.}
#'   \item{converged}{Logical value indicating whether the optimization converged
#'   before reaching `max.iterations`.}
#' }
#'
#' @details
#' - Uses the **L-BFGS-B** optimization method to efficiently find the best
#'   \eqn{\lambda} and \eqn{\tau} values.
#' - The optimization stops early if the relative change in \eqn{\lambda} and
#'   \eqn{\tau} is smaller than `tol`.
#' - If convergence is not achieved within `max.iterations`, the function
#'   prints a warning and returns the best estimate.
#'
#' @examples
#' \dontrun{
#' result <- optimize_marginal_log_likelihood(
#'   centered_kernel_mat_at_sampled, min_x, max_x, sampled_x,
#'   MC_iterations = 10000, max.iterations = 50, tol = 1e-4
#' )
#' print(result)
#' }
#'
#' @export


optimize_marginal_log_likelihood_nloptr <- function(centered_kernel_mat_at_sampled,
                                                    min_x,
                                                    max_x,
                                                    sampled_x,
                                                    initial_lambda = 0.1,
                                                    initial_tau = 0.1,
                                                    initial_w = rep(0, length(sampled_x)),
                                                    MC_iterations = 10000,
                                                    max.iterations = 10,
                                                    tol = 1e-4,  # Convergence tolerance
                                                    parallel_computing = TRUE) {

  converged <- FALSE
  t <- 1
  n <- length(sampled_x)

  # Generate matrix with each row independently sampled from Normal(0,1)
  set.seed(110)
  std_rnorm_matrix <- matrix(rnorm(MC_iterations * n, mean = 0, sd = 1),
                             nrow = MC_iterations,
                             ncol = n)

  lambda <- initial_lambda
  tau <- initial_tau
  w_vec <- initial_w

  dens_vec <- get_dens_wo_grid(centered_kernel_mat_at_sampled,
                               min_x,
                               max_x,
                               sampled_x,
                               lambda,
                               w_vec)

  p_vec <- dens_vec / sum(dens_vec)

  repeat {
    cat(paste0("Iteration: ", t, "\n"))

    lambda_old <- lambda
    tau_old <- tau

    # Define objective function
    objective_function <- function(params) {
      lambda <- params[1]
      tau <- params[2]

      -marginal_log_likelihood(
        centered_kernel_mat_at_sampled,
        sampled_x,
        min_x,
        max_x,
        p_vec,
        lambda,
        tau,
        std_rnorm_matrix,
        MC_iterations,
        parallel_computing)
    }

    # Define constraint function: lambda^2 / tau - 10^4.2 < 0
    constraint_function <- function(params) {
      lambda <- params[1]
      tau <- params[2]
      return(log(lambda) - log(tau) - log(10^(4.2)))  # log(lambda/tau) < log(10^(4.2))
    }


    cat(paste0("Initial lambda: ", 0.1, ", Initial tau: ", 0.1, "\n"))

    # Optimization using nloptr
    opt_result <- nloptr::nloptr(
      x0 = c(0.1, 0.1),
      eval_f = objective_function,
      eval_g_ineq = function(x) c(constraint_function(x)),  # Nonlinear inequality constraint
      lb = c(1e-2, 1e-6),
      ub = c(1e1, 1e1),
      opts = list("algorithm" = "NLOPT_LN_AUGLAG", "xtol_rel" = 1e-6)
    )

    lambda <- opt_result$solution[1]
    tau <- opt_result$solution[2]
    #lambda <- exp(log_lambda)
    #tau <- exp(log_tau)

    cat(paste0("Optimized lambda: ", lambda, ", tau: ", tau, ", MLL: ", -opt_result$objective,
               ", The ratio: ", lambda^2 / tau, "\n"))

    delta_lambda <- abs(lambda - lambda_old) / lambda_old
    delta_tau <- abs(tau - tau_old) / tau_old

    if (delta_lambda < tol & delta_tau < tol) {
      converged <- TRUE
      cat("✔ Convergence reached: Estimated Lambda and Tau by Maximum Marginal Likelihood are stable.\n")
      break
    }

    # Update weights
    w_vec <- get_weights_wo_grid_BBsolve(lambda_hat = lambda,
                                         tau_hat = tau,
                                         centered_kernel_mat_at_sampled = centered_kernel_mat_at_sampled,
                                         sampled_x = sampled_x,
                                         min_x = min_x,
                                         max_x = max_x,
                                         prior_variance_p_vector = p_vec,
                                         print_trace = FALSE)

    dens_vec <- get_dens_wo_grid(centered_kernel_mat_at_sampled,
                                 min_x,
                                 max_x,
                                 sampled_x,
                                 lambda,
                                 w_vec)

    p_vec <- dens_vec / sum(dens_vec)

    t <- t + 1

    if (t > max.iterations) {
      cat("Error: Max iterations reached without full convergence.\n")
      break
    }
  }

  return(list(lambda = lambda, tau = tau, max_marginal_log_likelihood = -opt_result$objective, converged = converged))
}
