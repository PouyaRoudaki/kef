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


lambda_hat <- 1
tau_hat <- 1/1350

kef_res <- kef(sampled_x,grid = x_grid,lambda = lambda_hat, tau = tau_hat)


centered_kernel_mat_at_sampled <- centered_kernel_matrix(first_vec_kernel = sampled_x,
                                                         second_vec_kernel = sampled_x,
                                                         centering_grid = x_grid,
                                                         hurst_coef = 0.5)

centered_kernel_mat_at_grid <- centered_kernel_matrix(first_vec_kernel = sampled_x,
                                                      second_vec_kernel = x_grid,
                                                      centering_grid = x_grid,
                                                      hurst_coef = 0.5)
centered_kernel_self_grid <- diag(centered_kernel_matrix(first_vec_kernel = x_grid,
                                                         second_vec_kernel = x_grid,
                                                         centering_grid = x_grid,
                                                         hurst_coef = 0.5))


########################################################################

lambda_grid <- 10^(seq(-1,1,by=0.1))
tau_grid <- 10^(seq(-4,1,by=0.1))

# Create a data frame with all possible combinations of lambda and tau
grid <- expand.grid(lambda = lambda_grid, tau = tau_grid)

# Calculate log10(tau) and log10(lambda)
grid$log10_tau <- log10(grid$tau)
grid$log10_lambda <- log10(grid$lambda)

# Filter the grid based on the condition
filtered_grid <- grid %>% filter(log10_tau >= log10_lambda - 4.2 )

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
