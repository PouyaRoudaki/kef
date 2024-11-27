# Start the timer
start_time <- proc.time()


mixture_weights <- c(1/2, 1/6, 1/6, 1/6)

# Define the parameters for the normal distributions
# First distribution: N(0, 1)
means <- c(0, -1, 0, 1)
sds <- c(1, 0.1, 0.1, 0.1)

sampled_x <- sort(normal_mixture(100, means, sds, mixture_weights))
x_grid <-  seq(-3.1,3.1,length.out = 400)
# centering_grid <- sampled_x This doesn't work because using this centering grid the kernel mean embedding is zero.
centering_grid <- runif(min = -3.1,max = 3.1,n = 100)


centered_kernel_mat_at_sampled_parallel_p <- centered_kernel_matrix_parallel(first_vec_kernel = sampled_x,
                                                                           second_vec_kernel = sampled_x,
                                                                           centering_grid = centering_grid,
                                                                           hurst_coef = 0.5)

centered_kernel_mat_at_sampled_parallel <- centered_kernel_matrix_parallel(first_vec_kernel,
                                                                           second_vec_kernel,
                                                                           centering_grid,
                                                                           hurst_coef)

centered_kernel_mat_at_sampled_p <- centered_kernel_matrix(first_vec_kernel = sampled_x,
                                                         second_vec_kernel = sampled_x,
                                                         centering_grid = centering_grid,
                                                         hurst_coef = 0.5)
# Check if matrices have the same elements
isTRUE(all.equal(centered_kernel_mat_at_sampled_parallel_p, centered_kernel_mat_at_sampled_p))




