# Load required packages
packages <- c("ks", "quantreg", "spatstat", "BB", "pracma", "akima",
              "tidyverse", "dplyr", "ggplot2", "parallel", "doParallel", "foreach", "kef")

lapply(packages, require, character.only = TRUE)

# Set up parallel backend
num_cores <- detectCores()  # Use all but one core
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Define the weights for the mixture distribution
mixture_weights = c(1/2, 1/10, 1/10, 1/10, 1/10, 1/10)

# Define the parameters for the normal distributions
# First distribution: N(0, 1)
means = c(0, -1, -0.5, 0, 0.5, 1)
sds = c(1, 0.1, 0.1, 0.1, 0.1, 0.1)


# Define parameters
min_x <- -3
max_x <- 3
#min_x <- 0
#max_x <- 1

n_grid <- 10000
grid <- seq(min_x, max_x, length.out = n_grid)

# Mixture Normal Example
params_mixture <- list(means = means, sds = sds, weights = mixture_weights)
density_characterization_mixture <- list(type = "mixture_normal", parameters = params_mixture)
true_density_grid <- true_density_function(grid, density_characterization_mixture)

# Uniform Example
#params_uniform <- list(min = min_x, max = max_x)
#density_characterization_uniform <- list(type = "uniform", parameters = params_uniform)
#true_density_grid <- true_density_function(grid, density_characterization_uniform)

n_iter <- 5000
result_list <- vector("list", n_iter)

# Parallel execution of iterations
result_list <- foreach(i = 1:n_iter, .packages = c("quantreg","ks", "pracma", "foreach", "kef")) %dopar% {
  cat("Iteration:", i, "\n")

  n_sample <- 100
  set.seed(i)
  sample <- sort(normal_mixture(n_sample, means, sds, mixture_weights))
  #sample <- sort(runif(n_sample,min = min_x, max = max_x))

  # Compute kernel matrices
  centered_kernel_mat_at_sampled <- centered_kernel_matrix(sample, sample, grid, 0.5)
  centered_kernel_mat_at_grid <- centered_kernel_matrix(sample, grid, grid, 0.5)
  centered_kernel_self_grid <- diag(centered_kernel_matrix(grid, grid, grid, 0.5))

  # Compute true weights
  w_true <- get_true_weights(sample, grid, true_density_grid)

  # Kernel Density Estimation (KDE)
  start_time <- Sys.time()
  estimated_density_grid <- kde(x = sample, h = hpi(sample), eval.points = grid)$estimate
  pi_time <- difftime(Sys.time(), start_time, units = "secs")
  pi_ise <- l2_ise(grid, true_density_grid, estimated_density_grid)

  # Square Cross Validation KDE
  start_time <- Sys.time()
  estimated_density_grid <- kde(x = sample, h = hscv(sample), eval.points = grid)$estimate
  scv_time <- difftime(Sys.time(), start_time, units = "secs")
  scv_ise <- l2_ise(grid, true_density_grid, estimated_density_grid)

  # Least Square Cross Validation KDE
  start_time <- Sys.time()
  estimated_density_grid <- kde(x = sample, h = hlscv(sample), eval.points = grid)$estimate
  lscv_time <- difftime(Sys.time(), start_time, units = "secs")
  lscv_ise <- l2_ise(grid, true_density_grid, estimated_density_grid)

  # Normal Scale KDE
  start_time <- Sys.time()
  estimated_density_grid <- kde(x = sample, h = hns(sample), eval.points = grid)$estimate
  ns_time <- difftime(Sys.time(), start_time, units = "secs")
  ns_ise <- l2_ise(grid, true_density_grid, estimated_density_grid)

  # Adaptive KDE (Adhoc)
  start_time <- Sys.time()
  estimated_density_grid <- akj(x = sample, z = grid)$dens
  adaptive_adhoc_time <- difftime(Sys.time(), start_time, units = "secs")
  adaptive_adhoc_ise <- l2_ise(grid, true_density_grid, estimated_density_grid)

  # KEF: Rule of Thumb
  kef_rot_lambda <- 1
  kef_rot_tau <- 1 / 1350
  kef_rot_ratio <- (kef_rot_lambda)^2 / kef_rot_tau
  kef_rot <- kef(sample, grid, lambda = kef_rot_lambda, tau = kef_rot_tau)
  kef_rot_time <- kef_rot$time
  estimated_density_grid <- kef_rot$probs_grid
  kef_rot_ise <- l2_ise(grid, true_density_grid, estimated_density_grid)
  kef_rot_se <- rkhs_se(w_hat_vec = kef_rot$weights, w_vec = w_true, kernel_matrix_at_samples = centered_kernel_mat_at_sampled)

  # KEF: Marginal Log Likelihood (MLL) Optimization
  start_time_mll <- Sys.time()
  optimized_mll <- optimize_marginal_log_likelihood(
    centered_kernel_mat_at_sampled, min_x, max_x, sample,
    initial_lambda = 1, initial_tau = 1, initial_w = rep(0, length(sample)),
    MC_iterations = 100000, max.iterations = 10, tol = 1e-1,
    parallel_computing = TRUE, seed = 4
  )
  kef_mll <- kef(sample, grid, lambda = optimized_mll$lambda, tau = optimized_mll$tau)
  kef_mll_time <- as.numeric(kef_mll$time) + as.numeric(difftime(Sys.time(), start_time_mll, units = "secs"))
  estimated_density_grid <- kef_mll$probs_grid
  kef_mll_ise <- l2_ise(grid, true_density_grid, estimated_density_grid)
  kef_mll_se <- rkhs_se(w_hat_vec = kef_mll$weights, w_vec = w_true, kernel_matrix_at_samples = centered_kernel_mat_at_sampled)

  kef_mll_lambda <- optimized_mll$lambda
  kef_mll_tau <- optimized_mll$tau
  kef_mll_ratio <- (kef_mll_lambda)^2 / kef_mll_tau

  # Store results in a dataframe
  result_iter_i <- data.frame(
    method = c("fixed_pi", "fixed_scv", "fixed_lscv", "fixed_ns",
               "adaptive_adhoc", "kef_rot", "kef_mll"),
    time = c(pi_time, scv_time, lscv_time, ns_time, adaptive_adhoc_time, kef_rot_time, kef_mll_time),
    ISE = c(pi_ise, scv_ise, lscv_ise, ns_ise, adaptive_adhoc_ise, kef_rot_ise, kef_mll_ise),
    RKHS_SE = c(NA, NA, NA, NA, NA, kef_rot_se, kef_mll_se),
    lambda = c(NA, NA, NA, NA, NA, kef_rot_lambda, kef_mll_lambda),
    tau = c(NA, NA, NA, NA, NA, kef_rot_tau, kef_mll_tau),
    ratio = c(NA, NA, NA, NA, NA, kef_rot_ratio, kef_mll_ratio)
  )

  return(result_iter_i)
}

# Stop parallel cluster
stopCluster(cl)

# Combine results into one data frame
df_combined <- bind_rows(result_list, .id = "iteration")

# Compute summary statistics
result_summary <- df_combined %>%
  group_by(method) %>%
  summarise(
    time_mean = mean(time, na.rm = TRUE),
    time_se = sd(time, na.rm = TRUE),
    MISE = mean(ISE, na.rm = TRUE),
    ISE_se = sd(ISE, na.rm = TRUE),
    RKHS_MSE = mean(RKHS_SE, na.rm = TRUE),
    RKHS_SE_se = sd(RKHS_SE, na.rm = TRUE),
    lambda_log_mean = mean(log10(lambda), na.rm = TRUE),
    lambda_log_se = sd(log10(lambda), na.rm = TRUE),
    tau_log_mean = mean(log10(tau), na.rm = TRUE),
    tau_log_se = sd(log10(tau), na.rm = TRUE),
    ratio_log_mean = mean(log10(ratio), na.rm = TRUE),
    ratio_log_se = sd(log10(ratio), na.rm = TRUE)
  )

# Save results
write.csv(result_summary, file = "summary_claw1.csv")

# View the summary
#View(result_summary)
