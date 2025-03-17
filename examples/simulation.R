packages <- c("ks", "quantreg","spatstat","BB","pracma", "akima", "tidyverse", "dplyr", "ggplot2", "parallel", "doParallel","foreach")


#install.packages(packages, dependencies = TRUE)
lapply(packages, require, character.only = TRUE)

# Specify the Density and its domain: ----

#The true density is **mixture normal with 3 peaks**.

# Define the weights for the mixture distribution
#mixture_weights <- c(1/2,1/6,1/6,1/6)

# Define the parameters for the normal distributions
# First distribution: N(0, 1)
#means <- c(0,-1, 0,1)
#sds <- c(1,0.1,0.1,0.1)

# Define the weights for the mixture distribution
mixture_weights = c(1/2, 1/10, 1/10, 1/10, 1/10, 1/10)

# Define the parameters for the normal distributions
# First distribution: N(0, 1)
means = c(0, -1, -0.5, 0, 0.5, 1)
sds = c(1, 0.1, 0.1, 0.1, 0.1, 0.1)

# Define the weights for the mixture distribution
#mixture_weights = c(1/2, 1/2)

# Define the parameters for the normal distributions
# First distribution: N(0, 1)
#means = c(0, 1.5)
#sds = c(0.5, 0.1)


# Define the domain of density and the corresponding grid
min_x <- -3
max_x <- 3

#min_x <- 0
#max_x <- 1

n_grid <- 10000

grid <-  seq(min_x,max_x,length.out = n_grid)

# Uniform Example
#params_uniform <- list(min = min_x, max = max_x)
#density_characterization_uniform <- list(type = "uniform", parameters = params_uniform)
#true_density_grid <- true_density_function(grid, density_characterization_uniform)

# Mixture Normal Example
params_mixture <- list(means = means, sds = sds, weights = mixture_weights)
density_characterization_mixture <- list(type = "mixture_normal", parameters = params_mixture)
true_density_grid <- true_density_function(grid, density_characterization_mixture)


# Iterative procedure ----
n_iter <- 1000
result_list <- list()

for (i in 1:n_iter) {
  cat("------------------------------------------------------------------------\n")
  cat(paste0("Iter ",i,":\n"))
  cat("------------------------------------------------------------------------\n")
  ## Take a Sample ----
  n_sample <- 100

  cat(paste0("✔ Take a Sample :\n"))
  set.seed(i)
  #sample <- sort(normal_mixture(n_sample, means, sds, mixture_weights))
  sample <- sort(runif(n_sample,min = min_x, max = max_x))

  ## Find the Cenetered Kernel Matrices ----


  centered_kernel_mat_at_sampled <- centered_kernel_matrix(first_vec_kernel = sample,                                       second_vec_kernel = sample,
                                                           centering_grid = grid,
                                                           hurst_coef = 0.5)
  centered_kernel_mat_at_grid <- centered_kernel_matrix(first_vec_kernel= sample,
                                                        second_vec_kernel = grid,
                                                        centering_grid = grid,
                                                        hurst_coef = 0.5)
  centered_kernel_self_grid <- diag(centered_kernel_matrix(first_vec_kernel = grid,
                                                           second_vec_kernel = grid,
                                                           centering_grid = grid,
                                                           hurst_coef = 0.5))


  ## Find the approximated true weights: ----

  w_true <- get_true_weights(sample,grid,true_density_grid)
  cat(paste0("✔ Find approximated true weights :\n"))

  ## Estimate the density and kef's weights: ----
  ### Fixed Bandwidth Density Estimation ----
  #### Fixed bandwidth: plug-in ----

  start_time <- Sys.time()
  estimated_density_grid <- kde(x = sample,h = hpi(sample), eval.points = grid)$estimate
  end_time <- Sys.time()
  # Compute total execution time
  pi_time <- difftime(end_time, start_time, units = "secs")
  pi_ise <- l2_ise(grid,true_density_grid,estimated_density_grid)
  cat(paste0("✔ Fixed plug-in :\n"))
  #### Fixed bandwidth: Square Cross Validation: ----

  start_time <- Sys.time()
  estimated_density_grid <- kde(x = sample,h = hscv(sample), eval.points = grid)$estimate
  end_time <- Sys.time()
  # Compute total execution time
  scv_time <- difftime(end_time, start_time, units = "secs")
  scv_ise <- l2_ise(grid,true_density_grid,estimated_density_grid)
  cat(paste0("✔ Fixed scv:\n"))
  #### Fixed bandwidth: Least Square Cross Validation: ----

  start_time <- Sys.time()
  estimated_density_grid <- kde(x = sample,h = hlscv(sample), eval.points = grid)$estimate
  end_time <- Sys.time()
  # Compute total execution time
  lscv_time <- difftime(end_time, start_time, units = "secs")
  lscv_ise <- l2_ise(grid,true_density_grid,estimated_density_grid)
  cat(paste0("✔ Fixed lscv :\n"))
  #### Fixed bandwidth: Normal Scale ----

  start_time <- Sys.time()
  estimated_density_grid <- kde(x = sample,h = hns(sample), eval.points = grid)$estimate
  end_time <- Sys.time()
  # Compute total execution time
  ns_time <- difftime(end_time, start_time, units = "secs")
  ns_ise <- l2_ise(grid,true_density_grid,estimated_density_grid)
  cat(paste0("✔ Fixed ns :\n"))
  ### Kernel Adaptive Density Estimation ----
  #### Adaptive bandwidth: Adhoc ----

  start_time <- Sys.time()
  estimated_density_grid <- akj(x = sample, z = grid)$dens
  end_time <- Sys.time()
  # Compute total execution time
  adaptive_adhoc_time <- difftime(end_time, start_time, units = "secs")
  adaptive_adhoc_ise <- l2_ise(grid,true_density_grid,estimated_density_grid)
  cat(paste0("✔ Adaptive adhoc:\n"))
  #### Adaptive bandwidth: CV optim ----

  adaptive_cv_optim <- function(sample, grid, alpha_init, kappa_init, lower_bounds, upper_bounds) {

    # Define the objective function to minimize (cv_error)
    objective_function <- function(params) {
      alpha <- params[1]
      kappa <- params[2]

      # Adaptive kernel density estimation
      adaptive_result <- akj(sample, grid, alpha = alpha, kappa = kappa)$dens

      # Roughness (L2 norm)
      roughness <- pracma::trapz(grid, (adaptive_result)^2)

      # Optimized Cross Product Term using proper indexing
      leave_one_out_densities <- vapply(seq_along(sample), function(i) {
        # Remove only the i-th occurrence, not all instances of sample[i]
        sample_wo_i <- sample[-i]

        # Compute Adaptive kernel density estimation without the i-th sample point
        adaptive_result_wo_i <- akj(x = sample_wo_i, z = sample[i], alpha = alpha, kappa = kappa)$dens

        # Extract density values
        densities_wo_i <- adaptive_result_wo_i

        return(densities_wo_i)
      }, numeric(1))  # Ensures output is numeric

      # Cross-validation error calculation
      cross_product_sum <- sum(leave_one_out_densities, na.rm = TRUE)
      cv_error <- roughness - (2 / length(sample)) * cross_product_sum

      return(cv_error)
    }

    # Optimization using Nelder-Mead or BFGS
    result <- optim(
      par = c(alpha_init, kappa_init),   # Initial values
      fn = objective_function,           # Function to minimize
      method = "L-BFGS-B",               # Bounded optimization
      lower = lower_bounds,              # Lower bounds for alpha and kappa
      upper = upper_bounds               # Upper bounds for alpha and kappa
    )

    # Extract best parameters
    best_alpha <- result$par[1]
    best_kappa <- result$par[2]
    min_cv_error <- result$value

    # Return results
    return(list(best_alpha = best_alpha, best_kappa = best_kappa, min_cv_error = min_cv_error))
  }



  start_time_adap <- Sys.time()

  alpha_init <- 0.5  # Initial guess for alpha
  kappa_init <- 0.5  # Initial guess for kappa
  lower_bounds <- c(0.01, 0.01)  # Lower bounds for alpha and kappa
  upper_bounds <- c(5, 5)  # Upper bounds for alpha and kappa

  # Adaptive hyperparam selection
  adaptive_cv_optim <- adaptive_cv_optim(sample, grid, alpha_init, kappa_init, lower_bounds, upper_bounds)

  # Adaptive best estimation
  estimated_density_grid <- akj(x = sample, z = grid, alpha = adaptive_cv_optim$best_alpha, kappa = adaptive_cv_optim$best_kappa)$dens
  end_time_adap <- Sys.time()
  adaptive_cv_time <- difftime(end_time_adap, start_time_adap, units = "secs")
  adaptive_cv_optim_ise <- l2_ise(grid,true_density_grid,estimated_density_grid)
  cat(paste0("✔ Adaptive optim mise:\n"))
  ###KEF ----

  #### KEF: Rule of thumb ----
  kef_rot_lambda <- 1
  kef_rot_tau <- 1/1350
  kef_rot_ratio <- (kef_rot_lambda)^2/kef_rot_tau
  kef_rot <- kef(sample, grid, lambda = kef_rot_lambda, tau = kef_rot_tau)
  kef_rot_time <- kef_rot$time
  estimated_density_grid <- kef_rot$probs_grid
  kef_rot_ise <- l2_ise(grid,true_density_grid,estimated_density_grid)
  kef_rot_se <- rkhs_se(w_hat_vec = kef_rot$weights,
                        w_vec = w_true,
                        kernel_matrix_at_samples = centered_kernel_mat_at_sampled)

  cat(paste0("✔ kef rot:\n"))
  #### KEF: MLL ----

  start_time_mll <- Sys.time()
  optimized_mll <- optimize_marginal_log_likelihood_new(centered_kernel_mat_at_sampled,
                                                        min_x = min_x,
                                                        max_x = max_x,
                                                        sample,
                                                        initial_lambda = 1,
                                                        initial_w = rep(0, length(sample)),
                                                        MC_iterations = 10000,
                                                        max.iterations = 10,
                                                        tol = 1e-3,
                                                        parallel_computing = TRUE,
                                                        seed = 4)
  end_time_mll <- Sys.time()
  kef_mll <- kef(sample, grid, lambda = optimized_mll$lambda, tau = optimized_mll$tau)
  kef_mll_time <- as.numeric(kef_mll$time) + as.numeric(difftime(end_time_mll, start_time_mll, units = "secs"))
  estimated_density_grid <- kef_mll$probs_grid
  kef_mll_ise <- l2_ise(grid,true_density_grid,estimated_density_grid)
  kef_mll_se <- rkhs_se(w_hat_vec = kef_mll$weights,
                        w_vec = w_true,
                        kernel_matrix_at_samples = centered_kernel_mat_at_sampled)

  kef_mll_lambda <- optimized_mll$lambda
  kef_mll_tau <- optimized_mll$tau
  kef_mll_ratio <- (kef_mll_lambda)^2/kef_mll_tau

  cat(paste0("✔ kef mll:\n"))
  #### KEF: Jackknife optimised MISE ----

  kef_jackknife_optim <- function(sample, grid, lambda_init, tau_init, lower_bounds, upper_bounds,loss_metric = "MISE") {

    if (loss_metric == "MISE"){
      # Define the objective function to minimize (cv_error)
      objective_function_mise <- function(params) {
        lambda <- params[1]
        tau <- params[2]

        # KEF density estimation
        kef_result <- kef(sample, grid, lambda, tau)

        # Roughness (L2 norm)
        roughness <- pracma::trapz(grid, (kef_result$probs_grid)^2)

        # Optimized Cross Product Term using proper indexing
        leave_one_out_densities <- vapply(seq_along(sample), function(i) {
          # Remove only the i-th occurrence, not all instances of sample[i]
          sample_wo_i <- sample[-i]

          # Compute KEF density estimation without the i-th sample point
          kef_result_wo_i <- kef(sample_wo_i, grid, lambda, tau)

          # Extract density values
          densities_wo_i <- kef_result_wo_i$probs_sample[i]

          return(densities_wo_i)
        }, numeric(1))  # Ensures output is numeric

        # Cross-validation error calculation
        cross_product_sum <- sum(leave_one_out_densities, na.rm = TRUE)
        cv_error <- roughness - (2 / length(sample)) * cross_product_sum

        return(cv_error)
      }

      # Optimization using Nelder-Mead or BFGS
      result <- optim(
        par = c(lambda_init, tau_init),   # Initial values
        fn = objective_function_mise,      # Function to minimize
        method = "L-BFGS-B",               # Bounded optimization
        lower = lower_bounds,              # Lower bounds for alpha and kappa
        upper = upper_bounds               # Upper bounds for alpha and kappa
      )

      # Extract best parameters
      best_lambda <- result$par[1]
      best_tau <- result$par[2]
      min_cv_error <- result$value

      # Return results
      return(list(best_lambda = best_lambda, best_tau = best_tau, min_cv_error = min_cv_error))
    } else if (loss_metric == "MSE"){
      # Define the objective function to minimize (cv_error)
      objective_function_mse <- function(params) {
        lambda <- params[1]
        tau <- params[2]

        # KEF density estimation using whole data
        kef_result <- kef(sample, grid, lambda, tau)

        # Estimated weights using full data
        w_estimated_full <- kef_result$weights

        # Optimized the Jackknife RKHS Mean Square Error term
        jackknife_error<- mean(vapply(seq_along(sample), function(i) {
          # Remove only the i-th occurrence, not all instances of sample[i]
          sample_wo_i <- sample[-i]

          # Compute KEF density estimation without the i-th sample point
          kef_result_wo_i <- kef(sample_wo_i, grid, lambda, tau)

          # Extract estimated weights with leave one out data
          w_estimated_wo_i <- kef_result_wo_i$weights

          # Estimate the Jackknife error
          jackknife_error_val <- rkhs_se(w_estimated_wo_i, w_estimated_full[-i],
                                         kernel_matrix_at_samples[-i,-i])

          return(jackknife_error_val)
        }, numeric(length(sample))))  # Ensures output is numeric

        return(jackknife_error)
      }

      # Optimization using Nelder-Mead or BFGS
      result <- optim(
        par = c(lambda_init, tau_init),   # Initial values
        fn = objective_function_mse,       # Function to minimize
        method = "L-BFGS-B",               # Bounded optimization
        lower = lower_bounds,              # Lower bounds for alpha and kappa
        upper = upper_bounds               # Upper bounds for alpha and kappa
      )

      # Extract best parameters
      best_lambda <- result$par[1]
      best_tau <- result$par[2]
      min_jackknife_error <- result$value

      # Return results
      return(list(best_lambda = best_lambda, best_tau = best_tau, min_jackknife_error = min_jackknife_error))
    } else {
      stop("loss_metric should be `MISE` or `MSE`")
    }
  }



  start_time_mise <- Sys.time()

  #lambda_init <- 0.1  # Initial guess for lambda
  #tau_init <- 0.01/1350  # Initial guess for tau
  #lower_bounds <- c(1e-2, 1e-6)  # Lower bounds for lambda and tau
  #upper_bounds <- c(30, 100)  # Upper bounds for lambda and tau

  # kef hyperparam selection
  #mise_jk_optim <- kef_jackknife_optim(sample, grid, lambda_init, tau_init, lower_bounds, upper_bounds,loss_metric = "MISE")

  #end_time_mise <- Sys.time()

  #kef_mise <- kef(sample, grid, lambda = mise_jk_optim$lambda, tau = mise_jk_optim$tau)
  #kef_mise_jk_time <- as.numeric(kef_mise$time) + as.numeric(difftime(end_time_mise, start_time_mise, units = "secs"))
  #estimated_density_grid <- kef_mise$probs_grid
  #kef_mise_jk_ise <- l2_ise(grid,true_density_grid,estimated_density_grid)
  #kef_mise_jk_se <- rkhs_se(w_hat_vec = kef_mise$weights,
  #                          w_vec = w_true,
  #                          kernel_matrix_at_samples = centered_kernel_mat_at_sampled)
  #cat(paste0("✔ kef optim mise:\n"))

  #### KEF: Jackknife optimised MSE ----


  #start_time_mse <- Sys.time()

  #lambda_init <- 0.1  # Initial guess for lambda
  #tau_init <- 0.01/1350  # Initial guess for tau
  #lower_bounds <- c(1e-2, 1e-6)  # Lower bounds for lambda and tau
  #upper_bounds <- c(30, 100)  # Upper bounds for lambda and tau

  # kef hyperparam selection
  #mse_jk_optim <- kef_jackknife_optim(sample, grid, lambda_init, tau_init, lower_bounds, upper_bounds,loss_metric = "MSE")

  #end_time_mse <- Sys.time()

  #kef_mse <- kef(sample, grid, lambda = mse_jk_optim$lambda, tau = mse_jk_optim$tau)
  #kef_mse_jk_time <- as.numeric(kef_mse$time) + as.numeric(difftime(end_time_mse, start_time_mse, units = "secs"))
  #estimated_density_grid <- kef_mse$probs_grid
  #kef_mse_jk_ise <- l2_ise(grid,true_density_grid,estimated_density_grid)
  #kef_mse_jk_se <- rkhs_se(w_hat_vec = kef_mse$weights,
  #                         w_vec = w_true,
  #                         kernel_matrix_at_samples = centered_kernel_mat_at_sampled)

  #cat(paste0("✔ kef optim mse:\n"))
  ## Final Result ----

  #Store the final result here.

  #result_iter_i <- data.frame(matrix(data = NA, nrow = 10, ncol = 4))
  #colnames(result_iter_i) <- c("method","time","ISE","RKHS_SE")
  #result_iter_i$method <- c("fixed_pi", "fixed_scv","fixed_lscv","fixed_ns",
  #                          "adaptive_adhoc","adaptive_cv_optim",
  #                          "kef_rot","kef_mll", "kef_jk_mise_optim", "kef_jk_mse_optim")

  #result_iter_i$time <- c(as.numeric(pi_time),as.numeric(scv_time),
  #                        as.numeric(lscv_time),as.numeric(ns_time),
  #                        as.numeric(adaptive_adhoc_time),
  #                        as.numeric(adaptive_cv_time),
  #                        as.numeric(kef_rot_time),
  #                        as.numeric(kef_mll_time),
  #                        as.numeric(kef_mise_jk_time),
  #                        as.numeric(kef_mse_jk_time))

  #result_iter_i$ISE <- c(pi_ise,scv_ise,lscv_ise,ns_ise,adaptive_adhoc_ise,
  #                       adaptive_cv_optim_ise,kef_rot_ise,kef_mll_ise,
  #                       kef_mise_jk_ise,kef_mse_jk_ise)

  #result_iter_i$SE <- c(NA,NA,NA,NA,NA,
  #                      NA,kef_rot_se,kef_mll_se,
  #                      kef_mise_jk_se,kef_mse_jk_se)


  result_iter_i <- data.frame(matrix(data = NA, nrow = 8, ncol = 7))
  colnames(result_iter_i) <- c("method","time","ISE","RKHS_SE","lambda","tau","ratio")
  result_iter_i$method <- c("fixed_pi", "fixed_scv","fixed_lscv","fixed_ns",
                            "adaptive_adhoc","adaptive_cv_optim",
                            "kef_rot","kef_mll")

  result_iter_i$time <- c(as.numeric(pi_time),as.numeric(scv_time),
                          as.numeric(lscv_time),as.numeric(ns_time),
                          as.numeric(adaptive_adhoc_time),
                          as.numeric(adaptive_cv_time),
                          as.numeric(kef_rot_time),
                          as.numeric(kef_mll_time))

  result_iter_i$ISE <- c(pi_ise,scv_ise,lscv_ise,ns_ise,adaptive_adhoc_ise,
                         adaptive_cv_optim_ise,kef_rot_ise,kef_mll_ise)

  result_iter_i$RKHS_SE <- c(NA,NA,NA,NA,NA,
                        NA,kef_rot_se,kef_mll_se)

  result_iter_i$lambda <- c(NA,NA,NA,NA,NA,
                             NA,kef_rot_lambda,kef_mll_lambda)

  result_iter_i$tau <- c(NA,NA,NA,NA,NA,
                             NA,kef_rot_tau,kef_mll_tau)

  result_iter_i$ratio <- c(NA,NA,NA,NA,NA,
                             NA,kef_rot_ratio,kef_mll_ratio)

  result_list[[i]] <- result_iter_i
  cat(paste0("✔ store the result:\n"))
}

#mix_norm_3_peaks <- result_list
claw <- result_list
#bimodal <- result_list
#uniform <- result_list

#save(mix_norm_3_peaks, file = "mix_norm_3_peaks.RData")
save(claw, file = "claw.RData")
#save(bimodal, file = "bimodal.RData")
#save(uniform, file = "uniform.RData")

#View(mix_norm_3_peaks[[i]])


# Assuming your list of data frames is named `df_list`
library(dplyr)

# Combine all data frames into one with an identifier for iterations
df_combined <- bind_rows(result_list, .id = "iteration")

# Compute mean and standard deviation for first three numeric columns
result_summary <- df_combined %>%
  group_by(method) %>%
  summarise(
    time_mean = mean(time, na.rm = TRUE),
    time_se = sd(time, na.rm = TRUE),
    MISE = mean(ISE, na.rm = TRUE),
    ISE_se = sd(ISE, na.rm = TRUE),
    RKHS_MSE = mean(RKHS_SE, na.rm = TRUE),
    RKHS_SE_se = sd(RKHS_SE, na.rm = TRUE),

    # Compute mean and standard deviation of log values for "lambda", "tau", and "ratio"
    lambda_log_mean = mean(log10(lambda), na.rm = TRUE),
    lambda_log_se = sd(log10(lambda), na.rm = TRUE),
    tau_log_mean = mean(log10(tau), na.rm = TRUE),
    tau_log_se = sd(log10(tau), na.rm = TRUE),
    ratio_log_mean = mean(log10(ratio), na.rm = TRUE),
    ratio_log_se = sd(log10(ratio), na.rm = TRUE)
  )

# View the final dataframe
View(result_summary)

# Store the summary dataframe
#summary_mix_norm_3_peaks <- result_summary
summary_claw <- result_summary
#summary_bimodal <- result_summary
#summary_uniform <- result_summary

#write.csv(summary_mix_norm_3_peaks, file = "summary_mix_norm_3_peaks.csv")
write.csv(summary_claw, file = "summary_claw.csv")
#write.csv(summary_bimodal, file = "summary_bimodal.csv")
#write.csv(summary_uniform, file = "summary_uniform.csv")
