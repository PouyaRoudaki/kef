library(quantreg)
library(ggplot2)
library(dplyr)
library(MASS)
library(BB)
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


lambda_hat <- 10
tau_hat <- 0.003981072

kefde <- kef(sampled_x,grid = x_grid,lambda = lambda_hat, tau = tau_hat)


# Convert the data to a data frame for use with ggplot2
plot_data <- data.frame(sampled_x = sampled_x, weights_hat = kefde$weights)

# Create the ggplot
p <- ggplot(plot_data) +
  geom_line(aes(x = sampled_x, y = weights_hat), color = "blue") +
  #geom_line(aes(x = sampled_x, y = weights_hat_init), color = "red") +
  labs(x = "Sampled x",
       y = "Weights Hat")+
  ggtitle(paste('Weights Hat vs Sampled x for lambda_hat =',
                format(lambda_hat,digits = 3,scientific = T),'and tau_hat =',format(tau_hat,digits = 3,scientific = T))) +
  theme_bw()

print(p)



kef_df <- data.frame(grid = sampled_x, kef_pdf = kefde$probs_sample)
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
kde_adaptive <- akj(sampled_x,sampled_x,kappa = 0.35,alpha = 0.9)
kde_adaptive_df <- data.frame(grid = sampled_x, kde_adaptive_pdf = kde_adaptive$dens)



ggplot() +
  geom_histogram(aes(x = sampled_x, y = ..density..), bins = 60, fill = 'gray', alpha = 1, color = 'black') +
  geom_line(data = true_density_df_sampled, aes(x = grid, y = true_pdf, color = 'True Density'), linewidth = 1) +
  #geom_point(data = true_density_df_sampled, aes(x = grid, y = weights_var, color = 'Weights Var'), size = 1) +
  geom_line(data = kde_adaptive_df, aes(x = grid, y = kde_adaptive_pdf, color = 'KDE Adaptive'), linewidth = 1) +
  geom_line(data = kef_df, aes(x = grid, y = kef_pdf, color = 'KEF'), linewidth = 1) +
  scale_color_manual(name = "Type of density", values = c('True Density' = 'red', 'KDE Adaptive' = 'blue', 'KEF' = 'orange')) +
  ggtitle(paste('Histogram and Kernel Density Estimate for lambda_hat =',
                format(lambda_hat,digits = 3,scientific = T),'and tau_hat =',format(tau_hat,digits = 3,scientific = T))) +
  xlab('Value') +
  ylab('Density') +
  theme_bw() +
  theme(legend.position = "bottom")


########################################################################

lambda_grid <- 10^(seq(-1,1,by=0.1))
tau_grid <- 10^(seq(-4,1,by=0.1))

# Create a data frame with all possible combinations of lambda and tau
grid <- expand.grid(lambda = lambda_grid, tau = tau_grid)

# Calculate log10(tau) and log10(lambda)
grid$log10_tau <- log10(grid$tau)
grid$log10_lambda <- log10(grid$lambda)

# Filter the grid based on the condition
filtered_grid <- grid %>% filter(log10_tau >= log10_lambda - 3 )

centered_kernel_mat_at_sampled <- centered_kernel_matrix(first_vec_kernel = sampled_x,
                                                         second_vec_kernel = sampled_x,
                                                         centering_grid = x_grid,
                                                         hurst_coef = 0.5) # add this to the maximum marginal likelihood function

lst_df <- compute_marginal_likelihood_grid_parallel(centered_kernel_mat_at_sampled,
                                                    min_x = -3.1,
                                                    max_x = 3.1,
                                                    sampled_x,
                                                    hyperparam_grid = filtered_grid,
                                                    initial_lambda = 1,
                                                    initial_w = rep(0, length(sampled_x)),
                                                    MC_iterations = 1000,
                                                    max_iterations = 5
)

library(parallel)
library(parallelly)
library(foreach)
lst_df_R <- compute_marginal_likelihood_grid_parallel_R(centered_kernel_mat_at_sampled,
         -3.1,
         3.1,
         sampled_x,
         hyperparam_grid = filtered_grid,
         initial_lambda = 1,
         initial_w = rep(0, length(sampled_x)),
         MC_iterations = 1000,
         max.iterations = 10,
         censoring = F)
