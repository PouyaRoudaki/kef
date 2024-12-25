lambda_grid <- as.vector(outer(1:9, 10^c(-2:2), "*"))
lambda_grid_trimmed <- lambda_grid[lambda_grid<60]


lambda_grid_trimmed <- seq(30,60,by =1)
library(ggplot2)
compute_marginal_likelihood_grid(centered_kernel_mat_at_sampled,
                                 min_x = -3.1,
                                 max_x = 3.1,
                                 sampled_x,
                                 lambda_grid = lambda_grid_trimmed,
                                 tau_hat = 1,
                                 initial_lambda_hat = 35,
                                 MC_iterations = 1000)



sample_mid_points <- get_middle_points_grid(-3.1, sampled_x, 3.1)
base_measure_weights <- sample_mid_points[-1] - sample_mid_points[-length(sample_mid_points)]

plot(base_measure_weights^(-1),x = sampled_x)
