# library
library(ggplot2)

# Define the weights for the mixture distribution
mixture_weights <- c(1/2,1/6,1/6,1/6)

# Define the parameters for the normal distributions
# First distribution: N(0, 1)
means <- c(0,-1, 0,1)
sds <- c(1,0.1,0.1,0.1)

set.seed(7)

sample <- sort(normal_mixture(100, means, sds, mixture_weights))
grid <-  seq(-3.1,3.1,length.out = 10000)

centered_kernel_mat_at_sampled <- centered_kernel_matrix(first_vec_kernel = sample,
                                                         second_vec_kernel = sample,
                                                         centering_grid = grid,
                                                         hurst_coef = 0.5)
centered_kernel_mat_at_grid <- centered_kernel_matrix(first_vec_kernel = sample,
                                                      second_vec_kernel = grid,
                                                      centering_grid = grid,
                                                      hurst_coef = 0.5)
centered_kernel_self_grid <- diag(centered_kernel_matrix(first_vec_kernel = grid,
                                                         second_vec_kernel = grid,
                                                         centering_grid = grid,
                                                         hurst_coef = 0.5))

# Mixture Normal Example
params_mixture <- list(means = means, sds = sds, weights = mixture_weights)
density_characterization_mixture <- list(type = "mixture_normal", parameters = params_mixture)
true_density_mixture <- true_density_function(grid, density_characterization_mixture)
plot(grid, true_density_mixture, type = "l", main = "Mixture Normal Density", xlab = "x", ylab = "Density")


w_vec <- get_true_weights(sample,grid,true_density_mixture)

plot(sample,w_vec, type = "l", main = "Mixture Normal weights", xlab = "x", ylab = "Weights")


approx_density <- get_dens_or_prob(
  centered_kernel_mat_at_sampled,
  centered_kernel_mat_at_grid,
  centered_kernel_self_grid,
  sample, grid,
  lambda_hat = 1,
  weight_hat_vec = w_vec,
  type_of_p_is_prob = FALSE,
  type_of_q_is_prob = FALSE
)$grid_x


# Create the dataframe
true_density_df <- data.frame(
  x = as.numeric(grid),
  true = as.numeric(true_density_mixture),
  approx = as.numeric(approx_density)
)


# Plot using ggplot2
ggplot(true_density_df, aes(x = x)) +
  geom_line(aes(y = true, color = "True Density"), size = 1) +
  geom_line(aes(y = approx, color = "Approx Density"), linetype = "dashed", size = 1) +
  labs(title = "True vs Approximate Density",
       x = "x",
       y = "Density",
       color = "Legend") +
  theme_bw() +
  scale_color_manual(values = c("True Density" = "blue", "Approx Density" = "red"))
