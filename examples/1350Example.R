#lambda_hat * (colSums(centered_kernel_mat_at_sampled) -
#                n * (prob_sampled_base) %*% t(centered_kernel_mat_at_sampled)) -
#  tau_hat * weight_hat_vec / prob_sampled

# for the case that lambda_hat = 1, tau_hat = 1/1350, n = 1000

#pdf("difficult_density_1.pdf")

library(quantreg)
library(ggplot2)
library(dplyr)
library(MASS)
# Example usage

set.seed(7)
# Define the weights for the mixture distribution
mixture_weights <- c(1/10,1/10,1/2,1/10,1/10,1/10)

# Define the parameters for the normal distributions
# First distribution: N(0, 1)
means <- c(0, -2,-1, 0,1,2)
sds <- c(1,0.1,0.1,0.1,0.1,0.1)

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

lambda_hat <- 1
tau_hat <- 1/1350



weights_hat_wo_grid <- get_weights_wo_grid(lambda_hat =lambda_hat,
                                           tau_hat = tau_hat,
                                           centered_kernel_mat_at_sampled,
                                           sampled_x = sampled_x,
                                           min_x = min(x_grid),
                                           max_x = max(x_grid),
                                           max_iteration = 1000,
                                           NRstepsize = 0.1,
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


# Precompute formatted values
lambda_hat_val <- format(lambda_hat, digits = 3, scientific = TRUE)
tau_hat_val <- format(tau_hat, digits = 3, scientific = TRUE)

# Plot
ggplot() +
  geom_histogram(aes(x = sampled_x, y = ..density..), bins = 60, fill = 'gray', alpha = 1, color = 'black') +
  geom_line(data = true_density_df, aes(x = grid, y = true_pdf, color = 'True Density'), linewidth = 1) +
  geom_line(data = kde_adaptive_df, aes(x = grid, y = kde_adaptive_pdf, color = 'KDE Adaptive'), linewidth = 1) +
  geom_line(data = kef_df, aes(x = grid, y = kef_pdf, color = 'KEF'), linewidth = 1) +
  scale_color_manual(name = "Type of density", values = c('True Density' = 'red', 'KDE Adaptive' = 'blue', 'KEF' = 'orange')) +
  ggtitle(
    bquote(
      "Histogram and Kernel Density Estimate for " ~
        lambda == .(lambda_hat_val) ~ "and" ~ tau == .(tau_hat_val)
    )
  ) +
  xlab('Value') +
  ylab('Density') +
  theme_bw() +
  theme(legend.position = "bottom")


dens <- get_dens_wo_grid(centered_kernel_mat_at_sampled,
                         min_x = -3.1,
                         max_x = 3.1,
                         sampled_x,
                         lambda_hat,
                         as.vector(weights_hat_wo_grid))

# Find the base measure of samples
sample_mid_points <- get_middle_points_grid(-3.1, sampled_x, 3.1)
base_measure_weights <- sample_mid_points[-1] - sample_mid_points[-length(sample_mid_points)]

dens_sampled_base <- dens * base_measure_weights

prob_sampled_base <- dens_sampled_base / sum(dens_sampled_base)
prob_sampled <- dens / sum(dens)

s1 <- colSums(centered_kernel_mat_at_sampled)
summary(s1)

s2 <-  as.numeric(dim(centered_kernel_mat_at_sampled)[1] * (prob_sampled_base) %*% t(centered_kernel_mat_at_sampled))
summary(s2)

s12<- colSums(centered_kernel_mat_at_sampled) - dim(centered_kernel_mat_at_sampled)[1] * (prob_sampled_base) %*% t(centered_kernel_mat_at_sampled)

s12t <- as.numeric(s12)
summary(as.numeric(s12))

plot(x = sampled_x,s1)
plot(x = sampled_x,as.numeric(s2))
plot(x = sampled_x,s12t)

s3 <- as.vector(weights_hat_wo_grid) / prob_sampled
summary(s3)

plot(x = sampled_x, s3)
summary(s12t)

plot(x = sampled_x,s12t * lambda_hat - s3*tau_hat)

# for the case that lambda_hat = 1, tau_hat = 1/1350, n = 1000

library(quantreg)
library(ggplot2)
library(dplyr)
library(MASS)
# Example usage

set.seed(7)
# Define the weights for the mixture distribution
mixture_weights <- c(1/10,1/10,1/2,1/10,1/10,1/10)

# Define the parameters for the normal distributions
# First distribution: N(0, 1)
means <- c(0, -2,-1, 0,1,2)
sds <- c(1,0.1,0.1,0.1,0.1,0.1)

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



lambda_hat <- 4.5
tau_hat <- 0.015


weights_hat_wo_grid_new <- get_weights_wo_grid(lambda_hat =lambda_hat,
                                           tau_hat = tau_hat,
                                           centered_kernel_mat_at_sampled,
                                           sampled_x = sampled_x,
                                           min_x = min(x_grid),
                                           max_x = max(x_grid),
                                           max_iteration = 1000,
                                           NRstepsize = 0.1,
                                           print_trace = T
)




#~ Convert the data to a data frame for use with ggplot2
plot_data <- data.frame(sampled_x = sampled_x, weights_hat = as.vector(weights_hat_wo_grid_new))

# Create the ggplot
p <- ggplot(plot_data, aes(x = sampled_x, y = weights_hat)) +
  geom_point(color  = "black") +
  geom_line(color = "blue") +
  labs(x = "Sampled x",
       y = "Weights Hat")+
  ggtitle(paste('Weights Hat vs Sampled x, lambda_hat =',
                format(lambda_hat,digits = 3,scientific = T),'and tau_hat =',format(tau_hat,digits = 3,scientific = T))) +
  theme_bw()

print(p)


probs_new <- get_dens_or_prob(centered_kernel_mat_at_sampled,
                          centered_kernel_mat_at_grid,
                          centered_kernel_self_grid,
                          sampled_x,x_grid,
                          lambda_hat, as.vector(weights_hat_wo_grid_new),
                          type_of_p_is_prob = FALSE,
                          type_of_q_is_prob = FALSE,
                          method_of_p_calculation = "ordinary")





#probs <- get_dens_wo_grid(centered_kernel_mat_at_sampled,
#                          -3.1,3.1,
#                          sampled_x,
#                          lambda_hat, as.vector(weights_hat_wo_grid))


kef_df <- data.frame(grid = x_grid, kef_pdf = probs_new$grid_x)
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


# Precompute formatted values
lambda_hat_val <- format(lambda_hat, digits = 3, scientific = TRUE)
tau_hat_val <- format(tau_hat, digits = 3, scientific = TRUE)

# Plot
ggplot() +
  geom_histogram(aes(x = sampled_x, y = ..density..), bins = 60, fill = 'gray', alpha = 1, color = 'black') +
  geom_line(data = true_density_df, aes(x = grid, y = true_pdf, color = 'True Density'), linewidth = 1) +
  geom_line(data = kde_adaptive_df, aes(x = grid, y = kde_adaptive_pdf, color = 'KDE Adaptive'), linewidth = 1) +
  geom_line(data = kef_df, aes(x = grid, y = kef_pdf, color = 'KEF'), linewidth = 1) +
  scale_color_manual(name = "Type of density", values = c('True Density' = 'red', 'KDE Adaptive' = 'blue', 'KEF' = 'orange')) +
  ggtitle(
    bquote(
      "Histogram and Kernel Density Estimate for " ~
        lambda == .(lambda_hat_val) ~ "and" ~ tau == .(tau_hat_val)
    )
  ) +
  xlab('Value') +
  ylab('Density') +
  theme_bw() +
  theme(legend.position = "bottom")



dens_new <- get_dens_wo_grid(centered_kernel_mat_at_sampled,
                         min_x = -3.1,
                         max_x = 3.1,
                         sampled_x,
                         lambda_hat,
                         as.vector(weights_hat_wo_grid_new))

# Find the base measure of samples
sample_mid_points <- get_middle_points_grid(-3.1, sampled_x, 3.1)
base_measure_weights <- sample_mid_points[-1] - sample_mid_points[-length(sample_mid_points)]

dens_sampled_base_new <- dens_new * base_measure_weights

prob_sampled_base_new <- dens_sampled_base_new / sum(dens_sampled_base_new)
prob_sampled_new <- dens_new / sum(dens_new)

s1_new <- colSums(centered_kernel_mat_at_sampled)
summary(s1_new)

s2_new <-  as.numeric(dim(centered_kernel_mat_at_sampled)[1] * (prob_sampled_base_new) %*% t(centered_kernel_mat_at_sampled))
summary(s2_new)

s12_new<- colSums(centered_kernel_mat_at_sampled) - dim(centered_kernel_mat_at_sampled)[1] * (prob_sampled_base_new) %*% t(centered_kernel_mat_at_sampled)

s12t_new <- as.numeric(s12_new)
summary(as.numeric(s12_new))

plot(x = sampled_x,s1_new)
plot(x = sampled_x,as.numeric(s2_new))
plot(x = sampled_x,s12t_new)

s3_new <- as.vector(weights_hat_wo_grid_new) / prob_sampled_new
summary(s3_new)

plot(x = sampled_x, s3_new)
summary(s12t_new)

summary(s3_new)[1]/summary(s12t_new)[1] * lambda_hat

plot(x = sampled_x,s12t_new * lambda_hat - s3_new*tau_hat)


dev.off()
