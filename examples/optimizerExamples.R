set.seed(7)

mixture_weights <- c(1/2, 1/6, 1/6, 1/6)

# Define the parameters for the normal distributions
# First distribution: N(0, 1)
means <- c(0, -1, 0, 1)
sds <- c(1, 0.1, 0.1, 0.1)

sampled_x <- sort(normal_mixture(100, means, sds, mixture_weights))
x_grid <-  seq(-3.1,3.1,length.out = 400)
# centering_grid <- sampled_x This doesn't work because using this centering grid the kernel mean embedding is zero.
#centering_grid <- runif(min = -3.1,max = 3.1,n = 1000)


centered_kernel_mat_at_sampled <- centered_kernel_matrix(first_vec_kernel = sampled_x,
                                                         second_vec_kernel = sampled_x,
                                                         centering_grid = x_grid,
                                                         hurst_coef = 0.5)

# Save the entire global environment to a file
#save.image(file = "my_environment.RData")



lambda_hat_grid <- c(1,2,80,3)
tau_hat_grid <- c(1)

library(pracma)

jackknife_normalised_weight_error_grid_inner_parallelized(centered_kernel_mat_at_sampled,
                                               sampled_x,
                                              min_x=-3.1,
                                              max_x=+3.1,
                                              lambda_hat_grid,
                                              tau_hat_grid,
                                              cloud_computing = F
                                              )

loocv_error_grid_inner_parallelized(centered_kernel_mat_at_sampled,
                                               sampled_x,
                                               min_x=-3.1,
                                               max_x=+3.1,
                                               lambda_hat_grid,
                                               tau_hat_grid,
                                               cloud_computing = F
)

# Load libraries
library(ggplot2)
library(dplyr)


# Replace 'data' with your dataset name
data <- result

# Calculate threshold (min mean + corresponding se)
min_mean <- min(data$jackknife_err_mean)
corresponding_se <- data$jackknife_err_se[which.min(data$jackknife_err_mean)]
threshold <- min_mean + corresponding_se



# Add a column to flag cells for highlighting
data <- data %>%
  mutate(
    highlight = jackknife_err_mean > threshold,
    is_min = jackknife_err_mean == min_mean,  # Identify the minimum value
  )



selected_tau <- data[which.min(data$jackknife_err_mean[which(data$highlight==T)]),"tau_hat"]
selected_lambda <- data[which.min(data$jackknife_err_mean[which(data$highlight==T)]),"lambda_hat"]

data <- data %>%
  mutate( selected = (tau_hat == selected_tau) & (lambda_hat == selected_lambda))

# Heatmap with highlights
ggplot(data, aes(x = tau_hat, y = lambda_hat)) +
  # Heatmap of mean values
  geom_tile(aes(fill = jackknife_err_mean)) +
  scale_fill_gradient(low = "blue", high = "red", name = "Mean Value") +

  # Highlight cells above threshold with yellow border
  geom_tile(data = subset(data, highlight == TRUE),
            aes(x = tau_hat, y = lambda_hat),
            fill = NA,color = "yellow", linewidth = 0.1) +

  # Highlight the tile with the minimum value using geom_point
  geom_point(data = subset(data, is_min == TRUE),
             aes(x = tau_hat, y = lambda_hat),
             fill= NA,color = "green", shape = 15, size = 2) +

  # Highlight the tile with the minimum value using geom_point
  geom_point(data = subset(data, selected == TRUE),
             aes(x = tau_hat, y = lambda_hat),
             fill= NA,color = "green", shape = 15, size = 2) +

  # Add labels and theme
  labs(
    title = "Heatmap of mean jackknife error",
    x = "tau_hat",
    y = "lambda_hat"
  ) +
  theme_minimal()


min_1se <- min(data$jackknife_err_mean[which(data$highlight==T)])
data[which.min(data$jackknife_err_mean[which(data$highlight==T)]),]
