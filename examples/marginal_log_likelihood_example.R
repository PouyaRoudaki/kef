library(quantreg)
library(ggplot2)
library(plotly)
library(akima)
library(dplyr)
library(MASS)
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
# centering_grid <- sampled_x This doesn't work because using this centering grid the kernel mean embedding is zero.
#centering_grid <- runif(min = -3.1,max = 3.1,n = 4000)

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


# Save the entire global environment to a file
#save.image(file = "my_environment.RData")
#lambda_hat <- 30
#tau_hat <- 0.4
# lambda^2/tau
#> 0.45^2/0.00015
#[1] 1350
#> 45^2/1.5
#[1] 1350
#> 4.5^2/0.015
#[1] 1350
#> 1^2/0.000333333

lambda_hat <- 45
tau_hat <- 1.5

lambda_hat <- 4.5
tau_hat <- 0.015

lambda_hat <- 1
tau_hat <- 1/1350

lambda_hat <- 0.45
tau_hat <- 0.00015

lambda_hat <- sqrt(1350)
tau_hat <- 1


lambda_hat <- sqrt(13.5)
tau_hat <- 0.01

lambda_hat <- 3
tau_hat <- 9/(1350)

lambda_hat <- 1
tau_hat <- 1/1350


weights_hat_wo_grid <- get_weights_wo_grid(lambda_hat =lambda_hat,
                                           tau_hat = tau_hat,
                                           centered_kernel_mat_at_sampled,
                                           sampled_x = sampled_x,
                                           min_x = min(x_grid),
                                           max_x = max(x_grid),
                                           print_trace = T
)



# Convert the data to a data frame for use with ggplot2
plot_data <- data.frame(sampled_x = sampled_x, weights_hat = as.vector(weights_hat_wo_grid))

# Create the ggplot
p <- ggplot(plot_data, aes(x = sampled_x, y = weights_hat)) +
  geom_point(color  = "black") +
  geom_line(color = "blue") +
  labs(x = "Sampled x",
       y = "Weights Hat")+
  ggtitle(paste('Weights Hat vs Sampled x for lambda_hat =',
                format(lambda_hat,digits = 3,scientific = T),'and tau_hat =',format(tau_hat,digits = 3,scientific = T))) +
  theme_bw()

print(p)


probs <- get_dens_or_prob(centered_kernel_mat_at_sampled,
                          centered_kernel_mat_at_grid,
                          centered_kernel_self_grid,
                          sampled_x,x_grid,
                          lambda_hat, as.vector(weights_hat_wo_grid),
                          type_of_p_is_prob = FALSE,
                          type_of_q_is_prob = FALSE,
                          method_of_p_calculation = "ordinary")





#probs <- get_dens_wo_grid(centered_kernel_mat_at_sampled,
#                          -3.1,3.1,
#                          sampled_x,
#                          lambda_hat, as.vector(weights_hat_wo_grid))


kef_df <- data.frame(grid = x_grid, kef_pdf = probs$grid_x)
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
kde_adaptive <- akj(sampled_x,x_grid,kappa = 0.35,alpha = 0.9)
kde_adaptive_df <- data.frame(grid = x_grid, kde_adaptive_pdf = kde_adaptive$dens)



ggplot() +
  geom_histogram(aes(x = sampled_x, y = ..density..), bins = 60, fill = 'gray', alpha = 1, color = 'black') +
  geom_line(data = true_density_df, aes(x = grid, y = true_pdf, color = 'True Density'), linewidth = 1) +
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

################################################################################

lambda_grid <- as.vector(outer(seq(1,9.5,1), 10^c(-1:1), "*"))
lambda_grid_trimmed <- lambda_grid[lambda_grid<20]

tau_grid <- as.vector(outer(seq(1,9.5,1), 10^c(-4:0), "*"))
tau_grid_trimmed <- tau_grid[tau_grid<10]

#lambda_grid_trimmed <- seq(30,60,by =1)
library(ggplot2)
# Specify the PDF output file
#pdf("output3.pdf")  # Adjust width and height as needed
#dev.off()

res_df <- compute_marginal_likelihood_grid(centered_kernel_mat_at_sampled,
                                 min_x = -6.1,
                                 max_x = 6.1,
                                 sampled_x,
                                 lambda_grid = lambda_grid_trimmed,
                                 tau_grid = tau_grid_trimmed,
                                 initial_lambda = 1,
                                 initial_w = rep(0, length(sampled_x)),
                                 MC_iterations = 1000,
                                 max.iterations = 4)


# Find the row with max MLL for each lambda
max_mll_per_lambda <- res_df %>%
  group_by(tau) %>%
  slice_max(order_by = mll, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(prop1 = lambda^2/tau, prop2 = lambda/tau, prop3 = log10(lambda)/log10(tau),
         suggested_lambda = sqrt(tau*1350))

# Print results
print(max_mll_per_lambda)

best_hyperparams <- res_df %>% arrange(-mll) %>% head(1)

best_weights <- get_weights_wo_grid(lambda_hat =best_hyperparams$lambda,
                                           tau_hat = best_hyperparams$tau,
                                           centered_kernel_mat_at_sampled,
                                           sampled_x = sampled_x,
                                           min_x = min(x_grid),
                                           max_x = max(x_grid),
                                           print_trace = T
)

best_weights <- get_weights_wo_grid(lambda_hat = 0.2,
                                    tau_hat = 0.0004,
                                    centered_kernel_mat_at_sampled,
                                    sampled_x = sampled_x,
                                    min_x = min(x_grid),
                                    max_x = max(x_grid),
                                    print_trace = T
)

plot(sampled_x,as.vector(best_weights))

best_probs <- get_dens_or_prob(centered_kernel_mat_at_sampled,
                          centered_kernel_mat_at_grid,
                          centered_kernel_self_grid,
                          sampled_x,x_grid,
                          lambda_hat, as.vector(best_weights),
                          type_of_p_is_prob = FALSE,
                          type_of_q_is_prob = FALSE,
                          method_of_p_calculation = "ordinary")

plot(sampled_x,as.vector(best_probs$sampled_x))
kef_df <- data.frame(grid = x_grid, kef_pdf = best_probs$grid_x)

ggplot() +
  geom_histogram(aes(x = sampled_x, y = ..density..), bins = 60, fill = 'gray', alpha = 1, color = 'black') +
  geom_line(data = true_density_df, aes(x = grid, y = true_pdf, color = 'True Density'), linewidth = 1) +
  #geom_point(data = true_density_df_sampled, aes(x = grid, y = weights_var, color = 'Weights Var'), size = 1) +
  geom_line(data = kde_adaptive_df, aes(x = grid, y = kde_adaptive_pdf, color = 'KDE Adaptive'), linewidth = 1) +
  geom_line(data = kef_df, aes(x = grid, y = kef_pdf, color = 'KEF'), linewidth = 1) +
  scale_color_manual(name = "Type of density", values = c('True Density' = 'red', 'KDE Adaptive' = 'blue', 'KEF' = 'orange')) +
  ggtitle(paste('Histogram and Kernel Density Estimate for lambda_hat =',
                format(best_hyperparams$lambda,digits = 3,scientific = T),'and tau_hat =',format(best_hyperparams$tau,digits = 3,scientific = T))) +
  xlab('Value') +
  ylab('Density') +
  theme_bw() +
  theme(legend.position = "bottom")

# Close the PDF device
dev.off()

#sample_mid_points <- get_middle_points_grid(-3.1, sampled_x, 3.1)
#base_measure_weights <- sample_mid_points[-1] - sample_mid_points[-length(sample_mid_points)]

#plot(base_measure_weights^(-1),x = sampled_x)

plot_ml <- data.frame(
  lambda = numeric(),
  marginal_log_likelihood = numeric()
)

plot_ml <- rbind(plot_ml, data.frame(
  lambda = 6,
  marginal_log_likelihood = 7
))
