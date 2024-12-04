

mixture_weights <- c(1/2, 1/6, 1/6, 1/6)

# Define the parameters for the normal distributions
# First distribution: N(0, 1)
means <- c(0, -1, 0, 1)
sds <- c(1, 0.1, 0.1, 0.1)

sampled_x <- sort(normal_mixture(100, means, sds, mixture_weights))
x_grid <-  seq(-3.1,3.1,length.out = 400)
# centering_grid <- sampled_x This doesn't work because using this centering grid the kernel mean embedding is zero.
centering_grid <- runif(min = -3.1,max = 3.1,n = 1000)


centered_kernel_mat_at_sampled <- centered_kernel_matrix_parallel(first_vec_kernel = sampled_x,
                                                         second_vec_kernel = sampled_x,
                                                         centering_grid = centering_grid,
                                                         hurst_coef = 0.5)
centered_kernel_mat_at_grid <- centered_kernel_matrix_parallel(first_vec_kernel = sampled_x,
                                                      second_vec_kernel = x_grid,
                                                      centering_grid = centering_grid,
                                                      hurst_coef = 0.5)
centered_kernel_self_grid <- diag(centered_kernel_matrix_parallel(first_vec_kernel = x_grid,
                                                         second_vec_kernel = x_grid,
                                                         centering_grid = centering_grid,
                                                         hurst_coef = 0.5))

# Save the entire global environment to a file
#save.image(file = "my_environment.RData")



lambda_hat_grid <- c(1)
tau_hat_grid <- c(0.1,1)

library(pracma)

jackknife_weight_error_grid(centered_kernel_mat_at_sampled,
                                        centered_kernel_mat_at_grid,
                                        centered_kernel_self_grid,
                                        sampled_x,
                                        x_grid,
                                        lambda_hat_grid,
                                        tau_hat_grid,
                                        type_of_p_is_prob = TRUE,
                                        type_of_q_is_prob = TRUE,
                                        method_of_p_calculation = "ordinary")



