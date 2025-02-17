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

lambda_grid <- 10^(seq(-1,2,by=0.1))
tau_grid <- 10^(seq(-4,1,by=0.1))

# Create a data frame with all possible combinations of lambda and tau
grid <- expand.grid(lambda = lambda_grid, tau = tau_grid)

# Calculate log10(tau) and log10(lambda)
grid$log10_tau <- log10(grid$tau)
grid$log10_lambda <- log10(grid$lambda)

# Filter the grid based on the condition
filtered_grid <- grid %>% filter(log10_tau >= log10_lambda - 4.2 )

# Print the filtered grid
View(filtered_grid)


#lambda_grid_trimmed <- seq(30,60,by =1)
library(ggplot2)
# Specify the PDF output file
#pdf("output3.pdf")  # Adjust width and height as needed
#dev.off()
library(doParallel)
library(foreach)

lst_df <- compute_marginal_likelihood_grid_parallel(centered_kernel_mat_at_sampled,
                                 min_x = -3.1,
                                 max_x = 3.1,
                                 sampled_x,
                                 hyperparam_grid = filtered_grid,
                                 initial_lambda = 1,
                                 initial_w = rep(0, length(sampled_x)),
                                 MC_iterations = 1000,
                                 max_iterations = 4
                                 )

result_df <- lst_df[[1]]

head_df <- result_df %>% arrange(desc(mll)) %>% head(200)



# Fit the linear model
model <- lm(log10(tau) ~ log10(lambda), data = head_df)

# Extract coefficients
coeffs <- coef(model)
intercept <- round(coeffs[1], 4)
slope <- round(coeffs[2], 4)

# Create equation text
eqn <- paste0("log10(tau) = ", intercept, " + ", slope, " * log10(lambda)")

# Plot
ggplot(head_df, aes(x = log10(lambda), y = log10(tau))) +
  geom_point(colour = "blue") +  # Scatter points
  geom_smooth(method = "lm", se = FALSE, colour = "red") +  # Regression line
  annotate("text", x = max(log10(head_df$lambda)), y = max(log10(head_df$tau)),
           label = eqn, hjust = 1, vjust = 1, size = 5, colour = "black") +
  labs(title = "Scatter Plot with Linear Fit",
       x = "log10(lambda)",
       y = "log10(tau)") +
  theme_minimal()


#df_max_mll_lambda <- my_df %>%
#  group_by(lambda) %>%
#  slice_max(mll, n = 5, with_ties = FALSE) %>%
#  ungroup()

# View the result
#View(df_max_mll_tau)
# Find the row with max MLL for each lambda
max_mll_per_lambda <- result_df %>%
  group_by(lambda) %>%
  slice_max(order_by = mll, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(prop1 = lambda^2/tau, prop2 = lambda/tau, prop3 = log10(lambda)/log10(tau), prop4 = log10(lambda)/(log10(tau)+log10(1350)),
         suggested_lambda = sqrt(tau*1350))

# Print results
print(max_mll_per_lambda)

best_hyperparams <- result_df %>% arrange(-mll) %>% head(1)
print(best_hyperparams)

lambda_hat  <- best_hyperparams$lambda
tau_hat <- best_hyperparams$tau
#lambda_hat <- 0.1
#tau_hat <- 0.0003
best_weights <- get_weights_wo_grid(lambda_hat = lambda_hat,
                                    tau_hat = tau_hat,
                                    centered_kernel_mat_at_sampled,
                                    sampled_x = sampled_x,
                                    min_x = min(x_grid),
                                    max_x = max(x_grid),
                                    print_trace = T
)

#plot(sampled_x,as.vector(best_weights))

best_probs <- get_dens_or_prob(centered_kernel_mat_at_sampled,
                          centered_kernel_mat_at_grid,
                          centered_kernel_self_grid,
                          sampled_x,x_grid,
                          lambda_hat, as.vector(best_weights),
                          type_of_p_is_prob = FALSE,
                          type_of_q_is_prob = FALSE,
                          method_of_p_calculation = "ordinary")

#plot(sampled_x,as.vector(best_probs$sampled_x))
kef_df <- data.frame(grid = x_grid, kef_pdf = best_probs$grid_x)

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

# Close the PDF device
lst_df_not_censored <- compute_marginal_likelihood_grid_parallel(centered_kernel_mat_at_sampled,
                                                    min_x = -3.1,
                                                    max_x = 3.1,
                                                    sampled_x,
                                                    hyperparam_grid = filtered_grid,
                                                    initial_lambda = 1,
                                                    initial_w = rep(0, length(sampled_x)),
                                                    MC_iterations = 1000,
                                                    max.iterations = 4,
                                                    censoring = F)

result_df_not_censored <- lst_df_not_censored[[1]]

best_hyperparams <- result_df_not_censored %>% arrange(-mll) %>% head(1)
print(best_hyperparams)

lambda_hat  <- best_hyperparams$lambda
tau_hat <- best_hyperparams$tau
#lambda_hat <- 0.1
#tau_hat <- 0.0003
best_weights <- get_weights_wo_grid(lambda_hat = lambda_hat,
                                    tau_hat = tau_hat,
                                    centered_kernel_mat_at_sampled,
                                    sampled_x = sampled_x,
                                    min_x = min(x_grid),
                                    max_x = max(x_grid),
                                    print_trace = T
)

#plot(sampled_x,as.vector(best_weights))

best_probs <- get_dens_or_prob(centered_kernel_mat_at_sampled,
                               centered_kernel_mat_at_grid,
                               centered_kernel_self_grid,
                               sampled_x,x_grid,
                               lambda_hat, as.vector(best_weights),
                               type_of_p_is_prob = FALSE,
                               type_of_q_is_prob = FALSE,
                               method_of_p_calculation = "ordinary")

#plot(sampled_x,as.vector(best_probs$sampled_x))
kef_df <- data.frame(grid = x_grid, kef_pdf = best_probs$grid_x)

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

sampled_x <- sort(normal_mixture(1000, means, sds, mixture_weights))
x_grid <-  seq(-3.1,3.1,length.out = 4000)
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

lambda_grid <- 10^(seq(-1,2,by=0.1))
tau_grid <- 10^(seq(-4,1,by=0.1))

# Create a data frame with all possible combinations of lambda and tau
grid <- expand.grid(lambda = lambda_grid, tau = tau_grid)

# Calculate log10(tau) and log10(lambda)
grid$log10_tau <- log10(grid$tau)
grid$log10_lambda <- log10(grid$lambda)

# Filter the grid based on the condition
filtered_grid <- grid[grid$log10_tau >= grid$log10_lambda - 4.2, ]

# Print the filtered grid
View(filtered_grid)


#lambda_grid_trimmed <- seq(30,60,by =1)
library(ggplot2)
# Specify the PDF output file
#pdf("output3.pdf")  # Adjust width and height as needed
#dev.off()
library(doParallel)
library(foreach)

lst_df <- compute_marginal_likelihood_grid_parallel(centered_kernel_mat_at_sampled,
                                                    min_x = -3.1,
                                                    max_x = 3.1,
                                                    sampled_x,
                                                    hyperparam_grid = filtered_grid,
                                                    initial_lambda = 1,
                                                    initial_w = rep(0, length(sampled_x)),
                                                    MC_iterations = 5000,
                                                    max.iterations = 5)

result_df <- lst_df[[1]]

head_df <- result_df %>% arrange(desc(mll)) %>% head(200)



# Fit the linear model
model <- lm(log10(tau) ~ log10(lambda), data = head_df)

# Extract coefficients
coeffs <- coef(model)
intercept <- round(coeffs[1], 4)
slope <- round(coeffs[2], 4)

# Create equation text
eqn <- paste0("log10(tau) = ", intercept, " + ", slope, " * log10(lambda)")

# Plot
ggplot(head_df, aes(x = log10(lambda), y = log10(tau))) +
  geom_point(colour = "blue") +  # Scatter points
  geom_smooth(method = "lm", se = FALSE, colour = "red") +  # Regression line
  annotate("text", x = max(log10(head_df$lambda)), y = max(log10(head_df$tau)),
           label = eqn, hjust = 1, vjust = 1, size = 5, colour = "black") +
  labs(title = "Scatter Plot with Linear Fit",
       x = "log10(lambda)",
       y = "log10(tau)") +
  theme_minimal()


#df_max_mll_lambda <- my_df %>%
#  group_by(lambda) %>%
#  slice_max(mll, n = 5, with_ties = FALSE) %>%
#  ungroup()

# View the result
#View(df_max_mll_tau)
# Find the row with max MLL for each lambda
max_mll_per_lambda <- result_df %>%
  group_by(lambda) %>%
  slice_max(order_by = mll, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(prop1 = lambda^2/tau, prop2 = lambda/tau, prop3 = log10(lambda)/log10(tau), prop4 = log10(lambda)/(log10(tau)+log10(1350)),
         suggested_lambda = sqrt(tau*1350))

# Print results
print(max_mll_per_lambda)

best_hyperparams <- result_df %>% arrange(-mll) %>% head(1)
print(best_hyperparams)

lambda_hat  <- best_hyperparams$lambda
tau_hat <- best_hyperparams$tau
#lambda_hat <- 0.1
#tau_hat <- 0.0003
best_weights <- get_weights_wo_grid(lambda_hat = lambda_hat,
                                    tau_hat = tau_hat,
                                    centered_kernel_mat_at_sampled,
                                    sampled_x = sampled_x,
                                    min_x = min(x_grid),
                                    max_x = max(x_grid),
                                    print_trace = T
)

#plot(sampled_x,as.vector(best_weights))

best_probs <- get_dens_or_prob(centered_kernel_mat_at_sampled,
                               centered_kernel_mat_at_grid,
                               centered_kernel_self_grid,
                               sampled_x,x_grid,
                               lambda_hat, as.vector(best_weights),
                               type_of_p_is_prob = FALSE,
                               type_of_q_is_prob = FALSE,
                               method_of_p_calculation = "ordinary")

#plot(sampled_x,as.vector(best_probs$sampled_x))
kef_df <- data.frame(grid = x_grid, kef_pdf = best_probs$grid_x)

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

# Close the PDF device
lst_df_not_censored <- compute_marginal_likelihood_grid_parallel(centered_kernel_mat_at_sampled,
                                                                 min_x = -3.1,
                                                                 max_x = 3.1,
                                                                 sampled_x,
                                                                 hyperparam_grid = filtered_grid,
                                                                 initial_lambda = 1,
                                                                 initial_w = rep(0, length(sampled_x)),
                                                                 MC_iterations = 5000,
                                                                 max.iterations = 4,
                                                                 censoring = F)

result_df_not_censored <- lst_df_not_censored[[1]]

best_hyperparams <- result_df_not_censored %>% arrange(-mll) %>% head(1)
print(best_hyperparams)

lambda_hat  <- best_hyperparams$lambda
tau_hat <- best_hyperparams$tau
#lambda_hat <- 0.1
#tau_hat <- 0.0003
best_weights <- get_weights_wo_grid(lambda_hat = lambda_hat,
                                    tau_hat = tau_hat,
                                    centered_kernel_mat_at_sampled,
                                    sampled_x = sampled_x,
                                    min_x = min(x_grid),
                                    max_x = max(x_grid),
                                    print_trace = T
)

#plot(sampled_x,as.vector(best_weights))

best_probs <- get_dens_or_prob(centered_kernel_mat_at_sampled,
                               centered_kernel_mat_at_grid,
                               centered_kernel_self_grid,
                               sampled_x,x_grid,
                               lambda_hat, as.vector(best_weights),
                               type_of_p_is_prob = FALSE,
                               type_of_q_is_prob = FALSE,
                               method_of_p_calculation = "ordinary")

#plot(sampled_x,as.vector(best_probs$sampled_x))
kef_df <- data.frame(grid = x_grid, kef_pdf = best_probs$grid_x)

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
