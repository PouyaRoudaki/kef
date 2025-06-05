library(spatstat)
library(ks)
library(ggplot2)

################################################################################
###########################     Our 3 modal Density 1000 #######################
################################################################################


set.seed(7)
# Define the weights for the mixture distribution
mixture_weights <- c(1/2,1/6,1/6,1/6)

# Define the parameters for the normal distributions
# First distribution: N(0, 1)
means <- c(0,-1, 0,1)
sds <- c(1,0.1,0.1,0.1)

sampled_x <- sort(normal_mixture(100, means, sds, mixture_weights))
n <- length(sampled_x)

x_grid <-  seq(-3.1,3.1,length.out = 4*n)

lambda <- 1
tau <- (lambda^2)/1350

kef_res <- kef(sampled_x,grid = x_grid,lambda = 1, tau = (lambda^2)/1350)

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
  geom_histogram(aes(x = sampled_x, y = ..density..), fill = 'gray', alpha = 1, color = 'black') +
  geom_line(data = true_density_df_sampled, aes(x = grid, y = true_pdf, color = 'True Density'), linewidth = 2) +
  #geom_point(data = true_density_df_sampled, aes(x = grid, y = weights_var, color = 'Weights Var'), size = 1) +
  geom_line(data = kde_adaptive_df, aes(x = grid, y = kde_adaptive_pdf, color = 'KDE Fixed'), linewidth = 2) +
  geom_line(data = kde_fixed_df, aes(x = grid, y = kde_fixed_pdf, color = 'KDE Adaptive'), linewidth = 2) +
  geom_line(data = kef_df, aes(x = grid, y = kef_pdf, color = 'KEF'), linewidth = 2) +
  scale_color_manual(name = "Type of density", values = c('True Density' = 'red', 'KDE Adaptive' = 'blue','KDE Fixed' = 'limegreen', 'KEF' = 'orange')) +
  ggtitle(paste('Kernel Density Estimate for lambda_hat =',
                format(lambda,digits = 3,scientific = T),'and tau_hat =',format(tau,digits = 3,scientific = T))) +
  xlab('Value') +
  ylab('Density') +
  theme_bw() +
  theme(legend.position = "bottom")



################################################################################
###########################   Without BASE MEASURE    ##########################
################################################################################


set.seed(7)
# Define the weights for the mixture distribution
mixture_weights <- c(1/2,1/6,1/6,1/6)

# Define the parameters for the normal distributions
# First distribution: N(0, 1)
means <- c(0,-1, 0,1)
sds <- c(1,0.1,0.1,0.1)

sampled_x <- sort(normal_mixture(100, means, sds, mixture_weights))
n <- length(sampled_x)

x_grid <-  seq(-3.1,3.1,length.out = 4*n)

lambda <- 1
tau <- (lambda^2)/1350

kef_res <- kef(sampled_x,grid = x_grid,lambda = 10, tau = 10,
                q_with_base = FALSE)

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
  geom_histogram(aes(x = sampled_x, y = ..density..), fill = 'gray', alpha = 1, color = 'black') +
  geom_line(data = true_density_df_sampled, aes(x = grid, y = true_pdf, color = 'True Density'), linewidth = 2) +
  #geom_point(data = true_density_df_sampled, aes(x = grid, y = weights_var, color = 'Weights Var'), size = 1) +
  geom_line(data = kde_adaptive_df, aes(x = grid, y = kde_adaptive_pdf, color = 'KDE Fixed'), linewidth = 2) +
  geom_line(data = kde_fixed_df, aes(x = grid, y = kde_fixed_pdf, color = 'KDE Adaptive'), linewidth = 2) +
  geom_line(data = kef_df, aes(x = grid, y = kef_pdf, color = 'KEF'), linewidth = 2) +
  scale_color_manual(name = "Type of density", values = c('True Density' = 'red', 'KDE Adaptive' = 'blue','KDE Fixed' = 'limegreen', 'KEF' = 'orange')) +
  ggtitle(paste('Kernel Density Estimate for lambda_hat =',
                format(lambda,digits = 3,scientific = T),'and tau_hat =',format(tau,digits = 3,scientific = T))) +
  xlab('Value') +
  ylab('Density') +
  theme_bw() +
  theme(legend.position = "bottom")

################################################################################
######################   With BASE MEASURE  both p and q  ######################
################################################################################


set.seed(7)
# Define the weights for the mixture distribution
mixture_weights <- c(1/2,1/6,1/6,1/6)

# Define the parameters for the normal distributions
# First distribution: N(0, 1)
means <- c(0,-1, 0,1)
sds <- c(1,0.1,0.1,0.1)

sampled_x <- sort(normal_mixture(100, means, sds, mixture_weights))
n <- length(sampled_x)

x_grid <-  seq(-3.1,3.1,length.out = 4*n)

lambda <- 1
tau <- (lambda^2)/1350

kef_res <- kef(sampled_x,grid = x_grid,lambda = 10, tau = 20,
               q_with_base = FALSE,normalised_q = TRUE,normalised_p = TRUE,
               p_with_base = FALSE)

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
  geom_histogram(aes(x = sampled_x, y = ..density..), fill = 'gray', alpha = 1, color = 'black') +
  geom_line(data = true_density_df_sampled, aes(x = grid, y = true_pdf, color = 'True Density'), linewidth = 2) +
  #geom_point(data = true_density_df_sampled, aes(x = grid, y = weights_var, color = 'Weights Var'), size = 1) +
  geom_line(data = kde_adaptive_df, aes(x = grid, y = kde_adaptive_pdf, color = 'KDE Fixed'), linewidth = 2) +
  geom_line(data = kde_fixed_df, aes(x = grid, y = kde_fixed_pdf, color = 'KDE Adaptive'), linewidth = 2) +
  geom_line(data = kef_df, aes(x = grid, y = kef_pdf, color = 'KEF'), linewidth = 2) +
  scale_color_manual(name = "Type of density", values = c('True Density' = 'red', 'KDE Adaptive' = 'blue','KDE Fixed' = 'limegreen', 'KEF' = 'orange')) +
  ggtitle(paste('Kernel Density Estimate for lambda_hat =',
                format(lambda,digits = 3,scientific = T),'and tau_hat =',format(tau,digits = 3,scientific = T))) +
  xlab('Value') +
  ylab('Density') +
  theme_bw() +
  theme(legend.position = "bottom")

################################################################################
######################   With BASE MEASURE  for q   ############################
################################################################################


set.seed(7)
# Define the weights for the mixture distribution
mixture_weights <- c(1/2,1/6,1/6,1/6)

# Define the parameters for the normal distributions
# First distribution: N(0, 1)
means <- c(0,-1, 0,1)
sds <- c(1,0.1,0.1,0.1)

sampled_x <- sort(normal_mixture(100, means, sds, mixture_weights))
n <- length(sampled_x)

x_grid <-  seq(-3.1,3.1,length.out = 4*n)

lambda <- 1
tau <- (lambda^2)/1350

kef_res <- kef(sampled_x,grid = x_grid,lambda = 1, tau = 1/1350,
               q_with_base = TRUE,normalised_q = FALSE,normalised_p = FALSE,
               p_with_base = FALSE)

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
  geom_histogram(aes(x = sampled_x, y = ..density..), fill = 'gray', alpha = 1, color = 'black') +
  geom_line(data = true_density_df_sampled, aes(x = grid, y = true_pdf, color = 'True Density'), linewidth = 2) +
  #geom_point(data = true_density_df_sampled, aes(x = grid, y = weights_var, color = 'Weights Var'), size = 1) +
  geom_line(data = kde_adaptive_df, aes(x = grid, y = kde_adaptive_pdf, color = 'KDE Fixed'), linewidth = 2) +
  geom_line(data = kde_fixed_df, aes(x = grid, y = kde_fixed_pdf, color = 'KDE Adaptive'), linewidth = 2) +
  geom_line(data = kef_df, aes(x = grid, y = kef_pdf, color = 'KEF'), linewidth = 2) +
  scale_color_manual(name = "Type of density", values = c('True Density' = 'red', 'KDE Adaptive' = 'blue','KDE Fixed' = 'limegreen', 'KEF' = 'orange')) +
  ggtitle(paste('Kernel Density Estimate for lambda_hat =',
                format(lambda,digits = 3,scientific = T),'and tau_hat =',format(tau,digits = 3,scientific = T))) +
  xlab('Value') +
  ylab('Density') +
  theme_bw() +
  theme(legend.position = "bottom")


################################################################################
######################   Without BASE MEASURE without normalising  ################################
################################################################################


set.seed(7)
# Define the weights for the mixture distribution
mixture_weights <- c(1/2,1/6,1/6,1/6)

# Define the parameters for the normal distributions
# First distribution: N(0, 1)
means <- c(0,-1, 0,1)
sds <- c(1,0.1,0.1,0.1)

sampled_x <- sort(normal_mixture(100, means, sds, mixture_weights))
n <- length(sampled_x)

x_grid <-  seq(-3.1,3.1,length.out = 4*n)

lambda <- 1
tau <- (lambda^2)/1350

kef_res <- kef(sampled_x,grid = x_grid,lambda = 1, tau = 1/1350,
               q_with_base = FALSE, normalised_q = FALSE, normalised_p = FALSE,
               p_with_base = FALSE)

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
  geom_histogram(aes(x = sampled_x, y = ..density..), fill = 'gray', alpha = 1, color = 'black') +
  geom_line(data = true_density_df_sampled, aes(x = grid, y = true_pdf, color = 'True Density'), linewidth = 2) +
  #geom_point(data = true_density_df_sampled, aes(x = grid, y = weights_var, color = 'Weights Var'), size = 1) +
  geom_line(data = kde_adaptive_df, aes(x = grid, y = kde_adaptive_pdf, color = 'KDE Fixed'), linewidth = 2) +
  geom_line(data = kde_fixed_df, aes(x = grid, y = kde_fixed_pdf, color = 'KDE Adaptive'), linewidth = 2) +
  geom_line(data = kef_df, aes(x = grid, y = kef_pdf, color = 'KEF'), linewidth = 2) +
  scale_color_manual(name = "Type of density", values = c('True Density' = 'red', 'KDE Adaptive' = 'blue','KDE Fixed' = 'limegreen', 'KEF' = 'orange')) +
  ggtitle(paste('Kernel Density Estimate for lambda_hat =',
                format(lambda,digits = 3,scientific = T),'and tau_hat =',format(tau,digits = 3,scientific = T))) +
  xlab('Value') +
  ylab('Density') +
  theme_bw() +
  theme(legend.position = "bottom")

