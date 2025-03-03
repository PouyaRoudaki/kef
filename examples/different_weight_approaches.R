library(quantreg)
library(ggplot2)
library(dplyr)
library(MASS)
library(BB)
library(ks)
# Example usage

set.seed(7)
# Define the weights for the mixture distribution
mixture_weights <- c(1/2,1/6,1/6,1/6)

# Define the parameters for the normal distributions
# First distribution: N(0, 1)
means <- c(0,-1, 0,1)
sds <- c(1,0.1,0.1,0.1)

sampled_x <- sort(normal_mixture(100, means, sds, mixture_weights))
x_grid <-  seq(-3.1,3.1,length.out = 400)
# For 10000 it takes 1446.795 secs

#lambda_hat <- 1
#tau_hat <- 0.00316227766016838

#lambda_hat <- 1
#tau_hat <- 1/1350

centered_kernel_mat_at_sampled <- centered_kernel_matrix(first_vec_kernel = sampled_x,
                                                         second_vec_kernel = sampled_x,
                                                         centering_grid = x_grid,
                                                         hurst_coef = 0.5)

optimal_kef <- optimize_marginal_log_likelihood(centered_kernel_mat_at_sampled,
                                 min_x = -3.1,
                                 max_x = 3.1,
                                 sampled_x,
                                 initial_lambda = 1,
                                 initial_w = rep(0, length(sampled_x)),
                                 MC_iterations = 10000,
                                 max.iterations = 3,
                                 tol = 1e-4,
                                 parallel_computing = TRUE
                                  )

kef_res <- kef(sampled_x,grid = x_grid,lambda = optimal_kef$lambda, tau = optimal_kef$tau)

kef_res$time

kef_df <- data.frame(grid = sampled_x, kef_pdf = kef_res$probs_sample)
#kef_df <- data.frame(grid = sampled_x, kef_pdf = probs)

# Define a matrix of normal densities for each mean and standard deviation
density_matrix <- sapply(seq_along(means), function(i) {
  dnorm(x_grid, mean = means[i], sd = sds[i])
})

# Define a matrix of normal densities for each mean and standard deviation
density_matrix_sampled <- sapply(seq_along(means), function(i) {
  dnorm(sampled_x, mean = means[i], sd = sds[i])
})

# Calculate the true density by taking the weighted sum of the columns
true_density <- density_matrix %*% mixture_weights

# Calculate the true density by taking the weighted sum of the columns
true_density_sampled <- density_matrix_sampled %*% mixture_weights

true_density_df <- data.frame(grid = x_grid, true_pdf = true_density)
true_density_df_sampled <- data.frame(grid = sampled_x, true_pdf = true_density_sampled)

# Perform the adaptive KDE
#kde_adaptive <- akj(sampled_x,sampled_x,kappa = 0.35,alpha = 0.9)
#kde_adaptive <- akj(sampled_x,sampled_x)
#kde_adaptive_df <- data.frame(grid = sampled_x, kde_adaptive_pdf = kde_adaptive$dens)

# Perform the adaptive KDE
kde_adaptive <- densityAdaptiveKernel(sampled_x)

kde_adaptive_df <- data.frame(grid = kde_adaptive$x, kde_adaptive_pdf = kde_adaptive$y)

kde_fixed <- kde(sampled_x,eval.points = sampled_x)
kde_fixed_df <- data.frame(grid = sampled_x, kde_fixed_pdf = kde_fixed$estimate)



ggplot() +
  geom_histogram(aes(x = sampled_x, y = ..density..), bins = 40, fill = 'gray', alpha = 1, color = 'black') +
  geom_line(data = true_density_df_sampled, aes(x = grid, y = true_pdf, color = 'True Density'), linewidth = 1) +
  #geom_point(data = true_density_df_sampled, aes(x = grid, y = weights_var, color = 'Weights Var'), size = 1) +
  geom_line(data = kde_adaptive_df, aes(x = grid, y = kde_adaptive_pdf, color = 'KDE Adaptive'), linewidth = 1) +
  geom_line(data = kde_fixed_df, aes(x = grid, y = kde_fixed_pdf, color = 'KDE Fixed'), linewidth = 1) +
  geom_line(data = kef_df, aes(x = grid, y = kef_pdf, color = 'KEF'), linewidth = 1) +
  scale_color_manual(name = "Type of density", values = c('True Density' = 'red', 'KDE Adaptive' = 'blue','KDE Fixed' = 'limegreen', 'KEF' = 'orange')) +
  ggtitle(paste('Histogram and Kernel Density Estimate for lambda_hat =',
                format(lambda_hat,digits = 3,scientific = T),'and tau_hat =',format(tau_hat,digits = 3,scientific = T))) +
  xlab('Value') +
  ylab('Density') +
  theme_bw() +
  theme(legend.position = "bottom")

################################################################################

library(quantreg)
library(ggplot2)
library(dplyr)
library(MASS)
library(BB)
library(ks)
# Example usage

set.seed(7)

sampled_x <- sort(runif(1000,0,1))
x_grid <-  seq(0,1,length.out = 4000)
# For 10000 it takes 1446.795 secs

centered_kernel_mat_at_sampled <- centered_kernel_matrix(first_vec_kernel = sampled_x,
                                                         second_vec_kernel = sampled_x,
                                                         centering_grid = x_grid,
                                                         hurst_coef = 0.5)

optimal_kef <- optimize_marginal_log_likelihood(centered_kernel_mat_at_sampled,
                                                min_x = 0,
                                                max_x = 1,
                                                sampled_x,
                                                initial_lambda = 1,
                                                initial_w = rep(0, length(sampled_x)),
                                                MC_iterations = 10000,
                                                max.iterations = 3,
                                                censoring = FALSE)

kef_res <- kef(sampled_x,grid = x_grid,lambda = optimal_kef$lambda, tau = optimal_kef$tau)

kef_res$time

kef_df <- data.frame(grid = sampled_x, kef_pdf = kef_res$probs_sample)
#kef_df <- data.frame(grid = sampled_x, kef_pdf = probs)



true_density_df_sampled <- data.frame(grid = x_grid, true_pdf =1)

# Perform the adaptive KDE
#kde_adaptive <- akj(sampled_x,sampled_x,kappa = 0.35,alpha = 0.9)
#kde_adaptive <- akj(sampled_x,sampled_x)
#kde_adaptive_df <- data.frame(grid = sampled_x, kde_adaptive_pdf = kde_adaptive$dens)

# Perform the adaptive KDE
kde_adaptive <- densityAdaptiveKernel(sampled_x)

kde_adaptive_df <- data.frame(grid = kde_adaptive$x, kde_adaptive_pdf = kde_adaptive$y)

kde_fixed <- kde(sampled_x,eval.points = sampled_x)
kde_fixed_df <- data.frame(grid = sampled_x, kde_fixed_pdf = kde_fixed$estimate)



ggplot() +
  geom_histogram(aes(x = sampled_x, y = ..density..), bins = 40, fill = 'gray', alpha = 1, color = 'black') +
  geom_line(data = true_density_df_sampled, aes(x = grid, y = true_pdf, color = 'True Density'), linewidth = 1) +
  #geom_point(data = true_density_df_sampled, aes(x = grid, y = weights_var, color = 'Weights Var'), size = 1) +
  geom_line(data = kde_adaptive_df, aes(x = grid, y = kde_adaptive_pdf, color = 'KDE Adaptive'), linewidth = 1) +
  geom_line(data = kde_fixed_df, aes(x = grid, y = kde_fixed_pdf, color = 'KDE Fixed'), linewidth = 1) +
  geom_line(data = kef_df, aes(x = grid, y = kef_pdf, color = 'KEF'), linewidth = 1) +
  scale_color_manual(name = "Type of density", values = c('True Density' = 'red', 'KDE Adaptive' = 'blue','KDE Fixed' = 'limegreen', 'KEF' = 'orange')) +
  ggtitle(paste('Histogram and Kernel Density Estimate for lambda_hat =',
                format(lambda_hat,digits = 3,scientific = T),'and tau_hat =',format(tau_hat,digits = 3,scientific = T))) +
  xlab('Value') +
  ylab('Density') +
  theme_bw() +
  theme(legend.position = "bottom")


