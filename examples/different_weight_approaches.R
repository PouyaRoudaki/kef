library(quantreg)
library(ggplot2)
library(dplyr)
library(MASS)
# Example usage

set.seed(7)
# Define the weights for the mixture distribution
mixture_weights <- c(1/2, 1/6, 1/6, 1/6)

# Define the parameters for the normal distributions
# First distribution: N(0, 1)
means <- c(0, -1, 0, 1)
sds <- c(1, 0.1, 0.1, 0.1)

sampled_x <- sort(normal_mixture(400, means, sds, mixture_weights))
x_grid <-  seq(-3.1,3.1,length.out = 400)
# centering_grid <- sampled_x This doesn't work because using this centering grid the kernel mean embedding is zero.
centering_grid <- runif(min = -3.1,max = 3.1,n = 400)

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

lambda_hat <- 45
tau_hat <- 1.5



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
  labs(title = "Weights Hat vs Sampled x",
       x = "Sampled x",
       y = "Weights Hat")+
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

kef_df <- data.frame(grid = x_grid, kef_pdf = probs$grid_x)

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
  ggtitle(paste('Histogram and Kernel Density Estimate for lambda_hat',
                round(lambda_hat,2),'and tau_hat',round(tau_hat,2))) +
  xlab('Value') +
  ylab('Density') +
  theme_bw() +
  theme(legend.position = "bottom")



# Define the range and number of points
lambda_hat_min <- 1
lambda_hat_max <- 100
num_points <- 10

# Generate the log scale grid
lambda_hat_grid <- seq(lambda_hat_min, lambda_hat_max, length.out = num_points)

# Define the range and number of points
tau_hat_min <- 0.1
tau_hat_max <- 100
num_points <- 10

# Generate the log scale grid
tau_hat_grid <- seq(tau_hat_min, tau_hat_max, length.out = num_points)

p_vec <- get_dens_or_prob(centered_kernel_mat_at_sampled, centered_kernel_mat_at_grid,
                          centered_kernel_self_grid, sampled_x,
                          x_grid, lambda_hat, weights_hat,
                          type_of_p_is_prob = TRUE,
                          type_of_q_is_prob = TRUE,
                          method_of_p_calculation = "neighborhood_grid")$sampled_x

# Open a connection to the file
#sink("examples/max_marginal_likelihood_result_4th_approach_20_iter_first_try.txt")
for (j in 1:10) {

  # printing step of iteration
  print(paste("Step:",j))
  cat('------------------------------------------------------------\n')
  # maximum marginal likelihood
  print("p_vec summary:")
  print(summary(p_vec))
  cat('------------------------------------------------------------\n')
  marginal_likelihood_df <- compute_marginal_likelihood_grid(centered_kernel_mat_at_sampled,
                                                             centered_kernel_mat_at_grid,
                                                             centered_kernel_self_grid,
                                                             sampled_x,
                                                             x_grid,
                                                             p_vec = p_vec,
                                                             lambda_grid = lambda_hat_grid,
                                                             tau_grid = tau_hat_grid,
                                                             MC_iterations = 1000,
                                                             type_of_p_is_prob = FALSE,
                                                             type_of_q_is_prob = FALSE,
                                                             method_of_p_calculation = "ordinary")

  # Arrange dataframe by the 3rd column in decreasing order
  marginal_likelihood_df_sorted <- marginal_likelihood_df %>%
    arrange(desc(marginal_log_likelihood))

  # Print header
  cat('------------------------------------------------------------\n')
  cat(sprintf("lambda_hat.    tau_hat.      marginal_log_likelihood. \n"))
  cat('------------------------------------------------------------\n')

  # Print first 20 rows in the specified format
  for(i in 1:20) {
    cat(sprintf("%4.2e        %8.2e     %8.2e\n", marginal_likelihood_df_sorted$lambda_hat[i],
                marginal_likelihood_df_sorted$tau_hat[i],
                marginal_likelihood_df_sorted$marginal_log_likelihood[i]))
  }


  # Updating p by tau_hat and lambda_hat which give us the maximum of marginal likelihood

  loop_continue <- TRUE
  while (loop_continue == TRUE) {
    error_flag <- FALSE

    tau_hat <- marginal_likelihood_df$tau_hat[which.max(marginal_likelihood_df$marginal_log_likelihood)]
    lambda_hat <- marginal_likelihood_df$lambda_hat[which.max(marginal_likelihood_df$marginal_log_likelihood)]

    print(paste('Chosen_lambda_hat =',lambda_hat,',', 'Chosen_tau_hat =', tau_hat))

    tryCatch(get_weights(lambda_hat =lambda_hat, tau_hat = tau_hat,
                         centered_kernel_mat_at_sampled, centered_kernel_mat_at_grid,
                         centered_kernel_self_grid, sampled_x = sampled_x,
                         x_grid = x_grid,
                         type_of_p_is_prob = TRUE,
                         type_of_q_is_prob = TRUE,
                         method_of_p_calculation = "neighborhood_grid"),
             error = function(e) {
               # Code to handle error
               cat("An error occurred because of non-invertible Hessian: ", e$message, "\n")
               marginal_likelihood_df$marginal_log_likelihood[which.max(marginal_likelihood_df$marginal_log_likelihood)] <<- -10^6
               error_flag <<- TRUE
               return(NULL)
             })
    #print(marginal_likelihood_df$marginal_log_likelihood[which.max(marginal_likelihood_df$marginal_log_likelihood)])
    #print(error_flag)
    if (!error_flag){
      loop_continue <- FALSE
    }


  }


  weights_hat <- get_weights(lambda_hat =lambda_hat, tau_hat = tau_hat,
                             centered_kernel_mat_at_sampled, centered_kernel_mat_at_grid,
                             centered_kernel_self_grid, x_grid = x_grid)


  probs <- get_dens_or_prob(centered_kernel_mat_at_sampled, centered_kernel_mat_at_grid,
                            centered_kernel_self_grid,
                            sampled_x,
                            x_grid,
                            lambda_hat,
                            weights_hat,
                            type_of_p_is_prob = TRUE,
                            type_of_q_is_prob = TRUE,
                            method_of_p_calculation = "neighborhood_grid")

  p_vec <- probs$sampled_x
  cat('============================================================\n')
}
#sink()

library(plotly)
plot_ly(x=marginal_likelihood_df$lambda_hat,
        y=marginal_likelihood_df$tau_hat,
        z=marginal_likelihood_df$marginal_log_likelihood, type="surface")
plot_ly(x = marginal_likelihood_df$lambda_hat,
        y = marginal_likelihood_df$tau_hat,
        z = marginal_likelihood_df$marginal_log_likelihood,
        type = "scatter3d",
        mode = "markers") %>%
  layout(scene = list(
    xaxis = list(title = "Lambda Hat"),
    yaxis = list(title = "Tau Hat"),
    zaxis = list(title = "Marginal Log Likelihood")
  ))

# Convert the data to a data frame for use with ggplot2
plot_data <- data.frame(sampled_x = sampled_x, weights_hat = weights_hat[1,])

# Create the ggplot
p <- ggplot(plot_data, aes(x = sampled_x, y = weights_hat)) +
  geom_point(color  = "black") +
  geom_line(color = "blue") +
  labs(title = "Weights Hat vs Sampled x",
       x = "Sampled x",
       y = "Weights Hat")+
  theme_bw()

print(p)


probs <- get_dens_or_prob(centered_kernel_mat_at_sampled, centered_kernel_mat_at_grid,
                          centered_kernel_self_grid, sampled_x,
                          x_grid, lambda_hat, weights_hat,
                          type_of_p_is_prob = FALSE,
                          type_of_q_is_prob = FALSE,
                          method_of_p_calculation = "ordinary")

kef_df <- data.frame(grid = x_grid, kef_pdf = probs$grid_x)

true_density <- sapply(x_grid, function(x){
  result <- mixture_weights[1] * dnorm(x, mean = means[1], sd = sds[1]) +
    mixture_weights[2] * dnorm(x, mean = means[2], sd = sds[2]) +
    mixture_weights[3] * dnorm(x, mean = means[3], sd = sds[3]) +
    mixture_weights[4] * dnorm(x, mean = means[4], sd = sds[4])

  return(result)
})


true_density_df <- data.frame(grid = x_grid, true_pdf = true_density)


# Perform the adaptive KDE
kde_adaptive <- akj(sampled_x,x_grid,kappa = 0.35,alpha = 0.9)
kde_adaptive_df <- data.frame(grid = x_grid, kde_adaptive_pdf = kde_adaptive$dens)


ggplot() +
  geom_histogram(aes(x = sampled_x, y = ..density..), bins = 60, fill = 'gray', alpha = 1, color = 'black') +
  geom_line(data = true_density_df, aes(x = grid, y = true_pdf, color = 'True Density'), linewidth = 1) +
  geom_line(data = kde_adaptive_df, aes(x = grid, y = kde_adaptive_pdf, color = 'KDE Adaptive'), linewidth = 1) +
  geom_line(data = kef_df, aes(x = grid, y = kef_pdf, color = 'KEF'), linewidth = 1) +
  scale_color_manual(name = "Type of density", values = c('True Density' = 'red', 'KDE Adaptive' = 'blue', 'KEF' = 'orange')) +
  ggtitle(paste('Histogram and Kernel Density Estimate for lambda_hat',
                round(lambda_hat,2),'and tau_hat',round(tau_hat,2))) +
  xlab('Value') +
  ylab('Density') +
  theme_bw() +
  theme(legend.position = "bottom")

## Repeat


# Define the range and number of points
lambda_hat_min <- 0.1
lambda_hat_max <- 20
num_points <- 20


# Generate the log scale grid
lambda_hat_grid <- seq(lambda_hat_min, lambda_hat_max, length.out = num_points)

# Define the range and number of points
tau_hat_min <- 0.1
tau_hat_max <- 20
num_points <- 20


# Generate the log scale grid
tau_hat_grid <- seq(tau_hat_min, tau_hat_max, length.out = num_points)

p_vec <- get_dens_or_prob(centered_kernel_mat_at_sampled, centered_kernel_mat_at_grid,
                          centered_kernel_self_grid, sampled_x,
                          x_grid, lambda_hat, weights_hat,
                          type_of_p_is_prob = TRUE,
                          type_of_q_is_prob = TRUE,
                          method_of_p_calculation = "neighborhood_grid")$sampled_x

# Open a connection to the file
#sink("examples/max_marginal_likelihood_result_4th_approach_20_iter_second_try.txt")
for (j in 1:10) {

  # printing step of iteration
  print(paste("Step:",j))
  cat('------------------------------------------------------------\n')
  # maximum marginal likelihood
  print("p_vec summary:")
  print(summary(p_vec))
  cat('------------------------------------------------------------\n')
  marginal_likelihood_df <- compute_marginal_likelihood_grid(centered_kernel_mat_at_sampled,
                                                             centered_kernel_mat_at_grid,
                                                             centered_kernel_self_grid,
                                                             sampled_x,
                                                             x_grid,
                                                             p_vec = p_vec,
                                                             lambda_grid = lambda_hat_grid,
                                                             tau_grid = tau_hat_grid,
                                                             MC_iterations = 1000,
                                                             type_of_p_is_prob = FALSE,
                                                             type_of_q_is_prob = FALSE,
                                                             method_of_p_calculation = "ordinary")

  # Arrange dataframe by the 3rd column in decreasing order
  marginal_likelihood_df_sorted <- marginal_likelihood_df %>%
    arrange(desc(marginal_log_likelihood))

  # Print header
  cat('------------------------------------------------------------\n')
  cat(sprintf("lambda_hat.    tau_hat.      marginal_log_likelihood. \n"))
  cat('------------------------------------------------------------\n')

  # Print first 20 rows in the specified format
  for(i in 1:20) {
    cat(sprintf("%4.2e        %8.2e     %8.2e\n", marginal_likelihood_df_sorted$lambda_hat[i],
                marginal_likelihood_df_sorted$tau_hat[i],
                marginal_likelihood_df_sorted$marginal_log_likelihood[i]))
  }


  # Updating p by tau_hat and lambda_hat which give us the maximum of marginal likelihood

  loop_continue <- TRUE
  while (loop_continue == TRUE) {
    error_flag <- FALSE

    tau_hat <- marginal_likelihood_df$tau_hat[which.max(marginal_likelihood_df$marginal_log_likelihood)]
    lambda_hat <- marginal_likelihood_df$lambda_hat[which.max(marginal_likelihood_df$marginal_log_likelihood)]

    print(paste('Chosen_lambda_hat =',lambda_hat,',', 'Chosen_tau_hat =', tau_hat))

    tryCatch(get_weights(lambda_hat =lambda_hat, tau_hat = tau_hat,
                         centered_kernel_mat_at_sampled, centered_kernel_mat_at_grid,
                         centered_kernel_self_grid, sampled_x = sampled_x,
                         x_grid = x_grid,
                         type_of_p_is_prob = TRUE,
                         type_of_q_is_prob = TRUE,
                         method_of_p_calculation = "neighborhood_grid"),
             error = function(e) {
               # Code to handle error
               cat("An error occurred because of non-invertible Hessian: ", e$message, "\n")
               marginal_likelihood_df$marginal_log_likelihood[which.max(marginal_likelihood_df$marginal_log_likelihood)] <<- -10^6
               error_flag <<- TRUE
               return(NULL)
             })
    #print(marginal_likelihood_df$marginal_log_likelihood[which.max(marginal_likelihood_df$marginal_log_likelihood)])
    #print(error_flag)
    if (!error_flag){
      loop_continue <- FALSE
    }


  }


  #weights_hat <- get_weights(lambda_hat =lambda_hat, tau_hat = tau_hat,
  #                           centered_kernel_mat_at_sampled, centered_kernel_mat_at_grid,
  #                           centered_kernel_self_grid, x_grid = x_grid)


  probs <- get_dens_or_prob(centered_kernel_mat_at_sampled, centered_kernel_mat_at_grid,
                            centered_kernel_self_grid,
                            sampled_x,
                            x_grid,
                            lambda_hat,
                            weights_hat,
                            type_of_p_is_prob = TRUE,
                            type_of_q_is_prob = TRUE,
                            method_of_p_calculation = "neighborhood_grid")

  p_vec <- probs$sampled_x
  cat('============================================================\n')
}
#sink()
#lambda_hat
#tau_hat

#weights_hat <- get_weights(lambda_hat =11.42, tau_hat = 0.53,
#                                                      centered_kernel_mat_at_sampled, centered_kernel_mat_at_grid,
#                                                      centered_kernel_self_grid, x_grid = x_grid)

# Convert the data to a data frame for use with ggplot2
plot_data <- data.frame(sampled_x = sampled_x, weights_hat = weights_hat[1,])

# Create the ggplot
p <- ggplot(plot_data, aes(x = sampled_x, y = weights_hat)) +
  geom_point(color  = "black") +
  geom_line(color = "blue") +
  labs(title = "Weights Hat vs Sampled x",
       x = "Sampled x",
       y = "Weights Hat")+
  theme_bw()

print(p)


probs <- get_dens_or_prob(centered_kernel_mat_at_sampled, centered_kernel_mat_at_grid,
                          centered_kernel_self_grid, sampled_x,
                          x_grid, lambda_hat, weights_hat,
                          type_of_p_is_prob = FALSE,
                          type_of_q_is_prob = FALSE,
                          method_of_p_calculation = "ordinary")

kef_df <- data.frame(grid = x_grid, kef_pdf = probs$grid_x)

true_density <- sapply(x_grid, function(x){
  result <- mixture_weights[1] * dnorm(x, mean = means[1], sd = sds[1]) +
    mixture_weights[2] * dnorm(x, mean = means[2], sd = sds[2]) +
    mixture_weights[3] * dnorm(x, mean = means[3], sd = sds[3]) +
    mixture_weights[4] * dnorm(x, mean = means[4], sd = sds[4])

  return(result)
})


true_density_df <- data.frame(grid = x_grid, true_pdf = true_density)


# Perform the adaptive KDE
kde_adaptive <- akj(sampled_x,x_grid,kappa = 0.35,alpha = 0.9)
kde_adaptive_df <- data.frame(grid = x_grid, kde_adaptive_pdf = kde_adaptive$dens)


ggplot() +
  geom_histogram(aes(x = sampled_x, y = ..density..), bins = 60, fill = 'gray', alpha = 1, color = 'black') +
  geom_line(data = true_density_df, aes(x = grid, y = true_pdf, color = 'True Density'), linewidth = 1) +
  geom_line(data = kde_adaptive_df, aes(x = grid, y = kde_adaptive_pdf, color = 'KDE Adaptive'), linewidth = 1) +
  geom_line(data = kef_df, aes(x = grid, y = kef_pdf, color = 'KEF'), linewidth = 1) +
  scale_color_manual(name = "Type of density", values = c('True Density' = 'red', 'KDE Adaptive' = 'blue', 'KEF' = 'orange')) +
  ggtitle(paste('Histogram and Kernel Density Estimate for lambda_hat',
                round(lambda_hat,2),'and tau_hat',round(tau_hat,2))) +
  xlab('Value') +
  ylab('Density') +
  theme_bw() +
  theme(legend.position = "bottom")


################################################################################
tau_hat <- 0.8
lambda_hat <- 1
weights_hat <- get_weights(lambda_hat =lambda_hat, tau_hat = tau_hat,
                                                      centered_kernel_mat_at_sampled, centered_kernel_mat_at_grid,
                                                      centered_kernel_self_grid, x_grid = x_grid)

# Convert the data to a data frame for use with ggplot2
plot_data <- data.frame(sampled_x = sampled_x, weights_hat = weights_hat[1,])

# Create the ggplot
p <- ggplot(plot_data, aes(x = sampled_x, y = weights_hat)) +
  geom_point(color  = "black") +
  geom_line(color = "blue") +
  labs(title = "Weights Hat vs Sampled x",
       x = "Sampled x",
       y = "Weights Hat")+
  theme_bw()

print(p)


probs <- get_dens_or_prob(centered_kernel_mat_at_sampled, centered_kernel_mat_at_grid,
                          centered_kernel_self_grid, sampled_x,
                          x_grid, lambda_hat, weights_hat,
                          type_of_p_is_prob = FALSE,
                          type_of_q_is_prob = FALSE,
                          method_of_p_calculation = "ordinary")

kef_df <- data.frame(grid = x_grid, kef_pdf = probs$grid_x)

true_density <- sapply(x_grid, function(x){
  result <- mixture_weights[1] * dnorm(x, mean = means[1], sd = sds[1]) +
    mixture_weights[2] * dnorm(x, mean = means[2], sd = sds[2]) +
    mixture_weights[3] * dnorm(x, mean = means[3], sd = sds[3]) +
    mixture_weights[4] * dnorm(x, mean = means[4], sd = sds[4])

  return(result)
})


true_density_df <- data.frame(grid = x_grid, true_pdf = true_density)


# Perform the adaptive KDE
kde_adaptive <- akj(sampled_x,x_grid,kappa = 0.35,alpha = 0.9)
kde_adaptive_df <- data.frame(grid = x_grid, kde_adaptive_pdf = kde_adaptive$dens)


ggplot() +
  geom_histogram(aes(x = sampled_x, y = ..density..), bins = 60, fill = 'gray', alpha = 1, color = 'black') +
  geom_line(data = true_density_df, aes(x = grid, y = true_pdf, color = 'True Density'), linewidth = 1) +
  geom_line(data = kde_adaptive_df, aes(x = grid, y = kde_adaptive_pdf, color = 'KDE Adaptive'), linewidth = 1) +
  geom_line(data = kef_df, aes(x = grid, y = kef_pdf, color = 'KEF'), linewidth = 1) +
  scale_color_manual(name = "Type of density", values = c('True Density' = 'red', 'KDE Adaptive' = 'blue', 'KEF' = 'orange')) +
  ggtitle(paste('Histogram and Kernel Density Estimate for lambda_hat',
                round(lambda_hat,2),'and tau_hat',round(tau_hat,4))) +
  xlab('Value') +
  ylab('Density') +
  theme_bw() +
  theme(legend.position = "bottom")


