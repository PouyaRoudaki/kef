library(ggplot2)
n_fixed_u <- 1000
n_grid <- 4000
n_fixed_x <- 500
lambda_true <- 1

set.seed(7)

fixed_u_vec <- seq(-3.1, 3.1,length.out = n_fixed_u)
#x_grid <-  seq(-4,4,length.out =40)
x_grid_cv <- seq(-3.1, 3.1,length.out = n_grid)

w_u_vec <- sapply(fixed_u_vec,function(x){-0.2*exp(-(10^(-6))*abs(x)^40)*(3.1-abs(x)^1)^(0.01)*(sin((-abs(x)+pi/4)*(2*pi)))})
w_u_df <- data.frame(x = fixed_u_vec, weights = w_u_vec)

sigma_w <- sd(w_u_vec)
tau_true <- var(w_u_vec)^-1

# Create the plot
ggplot(w_u_df, aes(x = x, y = weights)) +
  geom_hline(yintercept = 0, color = "darkgray", linetype = "solid", linewidth = 0.5) +
  geom_vline(xintercept = 0, color = "darkgray", linetype = "solid", linewidth = 0.5) +
  geom_line(color = "lightblue", linewidth = .1) +
  geom_point(color = "darkblue", size = 0.5) +
  ggtitle(bquote("Weights:" ~
                   m ==.(n_fixed_u) * "," ~
                   lambda == .(lambda_true) * "," ~
                   tau == .(tau_true))) +
  xlab("x") +
  ylab("Weights") +
  theme_bw()



centering_grid <- x_grid_cv



centered_kernel_mat_at_sampled_cv <- centered_kernel_matrix(first_vec_kernel = fixed_u_vec,
                                                            second_vec_kernel = fixed_u_vec,
                                                            centering_grid = centering_grid,
                                                            hurst_coef = 0.5)
centered_kernel_mat_at_grid_cv <- centered_kernel_matrix(first_vec_kernel = fixed_u_vec,
                                                         second_vec_kernel = x_grid_cv,
                                                         centering_grid = centering_grid,
                                                         hurst_coef = 0.5)
centerd_kernel_self_grid_cv <- diag(centered_kernel_matrix(first_vec_kernel = x_grid_cv,
                                                           second_vec_kernel = x_grid_cv,
                                                           centering_grid = centering_grid,
                                                           hurst_coef = 0.5))

probs_cv <- get_dens_or_prob(centered_kernel_mat_at_sampled_cv,
                             centered_kernel_mat_at_grid_cv,
                             centerd_kernel_self_grid_cv,
                             fixed_u_vec,x_grid_cv,
                             lambda_true, w_u_vec,
                             type_of_p_is_prob = FALSE,
                             type_of_q_is_prob = FALSE,
                             method_of_p_calculation = "ordinary")

true_density_df <- data.frame(grid = x_grid_cv, true_pdf = probs_cv$grid_x)

ggplot() +
  geom_line(data = true_density_df, aes(x = grid, y = true_pdf, color = 'True Density'), linewidth = 1) +
  #geom_line(data = kde_adaptive_df, aes(x = grid, y = kde_adaptive_pdf, color = 'KDE Adaptive'), linewidth = 1) +
  #geom_line(data = kef_df, aes(x = grid, y = kef_pdf, color = 'KEF'), linewidth = 1) +
  scale_color_manual(name = "Type of density", values = c('True Density' = 'blue')) +
  ggtitle(paste('Histogram and Kernel Density Estimate for lambda =', lambda_true,'and tau =',round(tau_true,2))) +
  xlab('u') +
  ylab('Density') +
  theme_bw() +
  theme(legend.position = "bottom")


## Take samples X1, . . . , Xn ∼ p(x) Using invert cdf transform, take samples X1, . . . , Xn ∼ p(x).
# Normalize the density values to obtain probabilities
p_x_vec <- probs_cv$grid_x

probabilities <- p_x_vec / sum(p_x_vec)
# Create the cumulative density function (CDF)
cdf <- cumsum(probabilities)
# plot(cdf)
# I took n number of samples from x which n << m
sampled_x_vec <- numeric(n_fixed_x)
sampled_indeces <- integer()
# Invert cdf transform
j = 0
while (j <= n_fixed_x) {
  u <- runif(1)
  sample_index <- which(cdf >= u)[1]
  if(!(sample_index %in% sampled_indeces)){
    sampled_indeces <- c(sampled_indeces,sample_index)
    sampled_x_vec[j] <- x_grid_cv[sample_index]
    j = j+1
  }
}

# Sort sampled x
sampled_x_vec <- sort(sampled_x_vec)

hist(sampled_x_vec,breaks = 50)

centering_grid <- x_grid_cv
centered_kernel_mat_x_u_cv <- centered_kernel_matrix(first_vec_kernel = fixed_u_vec,
                                                     second_vec_kernel = sampled_x_vec,
                                                     centering_grid = centering_grid,
                                                     hurst_coef = 0.5)
theta_n_evaluated_at_sampled_x <- (w_u_vec %*% centered_kernel_mat_x_u_cv)

centered_kernel_mat_x_x_cv <- centered_kernel_matrix(first_vec_kernel = sampled_x_vec,
                                                     second_vec_kernel = sampled_x_vec,
                                                     centering_grid = centering_grid,
                                                     hurst_coef = 0.5)




w_x_vec <- as.numeric(solve(centered_kernel_mat_x_x_cv) %*% t(theta_n_evaluated_at_sampled_x))
w_x_df <- data.frame( x = sampled_x_vec, weights = w_x_vec)

ggplot(w_x_df, aes(x = x, y = weights)) +
  geom_hline(yintercept = 0, color = "darkgray", linetype = "solid", linewidth = 0.5) +
  geom_vline(xintercept = 0, color = "darkgray", linetype = "solid", linewidth = 0.5) +
  geom_line(aes(color = "w_x_sampled"), linewidth = .5) +
  geom_point(color = "darkblue", size = 1) +
  geom_point(data = w_u_df, aes(x = x, y = weights), color = "darkgreen", size = 1) +
  geom_line(data = w_u_df, aes(x = x, y = weights, color = "w_u_grid"), linewidth = 0.3) +
  scale_color_manual(name = "Type of weights", values = c('w_x_sampled' = 'lightblue', 'w_u_grid' = 'chartreuse4')) +
  ggtitle(bquote( w[x]* ":" ~
                    m ==.(n_fixed_u) * "," ~
                    lambda == .(lambda_true) * "," ~
                    tau[w] == .(tau_true))) +
  xlab("x") +
  ylab(bquote( w[x] )) +
  theme_bw()


centering_grid_cv <- x_grid_cv



centered_kernel_mat_x_x_cv <- centered_kernel_matrix(first_vec_kernel = sampled_x_vec,
                                                     second_vec_kernel = sampled_x_vec,
                                                     centering_grid = centering_grid,
                                                     hurst_coef = 0.5)
centered_kernel_mat_x_grid_cv <- centered_kernel_matrix(first_vec_kernel = sampled_x_vec,
                                                        second_vec_kernel = x_grid_cv,
                                                        centering_grid = centering_grid,
                                                        hurst_coef = 0.5)
centerd_kernel_self_grid_cv <- diag(centered_kernel_matrix(first_vec_kernel = x_grid_cv,
                                                           second_vec_kernel = x_grid_cv,
                                                           centering_grid = centering_grid,
                                                           hurst_coef = 0.5))

probs_n_cv <- get_dens_or_prob(centered_kernel_mat_x_x_cv,
                               centered_kernel_mat_x_grid_cv,
                               centerd_kernel_self_grid_cv,
                               sampled_x_vec,x_grid_cv,
                               lambda_true, w_x_vec,
                               type_of_p_is_prob = FALSE,
                               type_of_q_is_prob = FALSE,
                               method_of_p_calculation = "ordinary")

true_n_density_df <- data.frame(grid = x_grid_cv, n_pdf = probs_n_cv$grid_x)

ggplot() +
  geom_line(data = true_density_df, aes(x = grid, y = true_pdf, color = 'p(x)'), linewidth = 1, alpha =1) +
  #geom_line(data = kde_adaptive_df, aes(x = grid, y = kde_adaptive_pdf, color = 'KDE Adaptive'), linewidth = 1) +
  geom_line(data = true_n_density_df, aes(x = grid, y = n_pdf, color = 'p_n(x)'), linewidth = 1, alpha = 1) +
  scale_color_manual(name = "Type of density", values = c('p(x)' = 'blue', 'p_n(x)' = 'chartreuse3')) +
  ggtitle(paste('Histogram and Kernel Density Estimate for lambda =', lambda_true, 'and tau =', round(1/var(w_x_vec),2))) +
  xlab('u') +
  ylab('Density') +
  theme_bw() +
  theme(legend.position = "bottom")

################################################################################

lambda_hat <- 0.055
tau_hat <- 0.055^2/1350

x_grid <- seq(-3.1, 3.1,length.out = 400)


centering_grid_cv <- x_grid



centered_kernel_mat_at_sampled <- centered_kernel_matrix(first_vec_kernel = sampled_x_vec,
                                                         second_vec_kernel = sampled_x_vec,
                                                         centering_grid = centering_grid_cv,
                                                         hurst_coef = 0.5)
centered_kernel_mat_at_grid <- centered_kernel_matrix(first_vec_kernel = sampled_x_vec,
                                                      second_vec_kernel = x_grid,
                                                      centering_grid = centering_grid_cv,
                                                      hurst_coef = 0.5)
centerd_kernel_self_grid <- diag(centered_kernel_matrix(first_vec_kernel = x_grid,
                                                        second_vec_kernel = x_grid,
                                                        centering_grid = centering_grid_cv,
                                                        hurst_coef = 0.5))



weights_hat <- get_weights_wo_grid(lambda_hat =lambda_hat,
                                   tau_hat = tau_hat,
                                   centered_kernel_mat_at_sampled,
                                   sampled_x = sampled_x_vec,
                                   min_x = min(x_grid),
                                   max_x = max(x_grid),
                                   print_trace = T
)



# Convert the data to a data frame for use with ggplot2
plot_data <- data.frame(sampled_x = sampled_x_vec, weights_hat = weights_hat[1,])

# Create the ggplot
p <- ggplot(plot_data, aes(x = sampled_x, y = weights_hat)) +
  geom_point(color = "darkred", size = 1) +
  geom_line(aes(color = "weights_estimation"), linewidth = 0.6) +
  geom_point(data = w_x_df, aes(x = x, y = weights), color = "black", size = 1) +
  geom_line(data = w_x_df, aes(x = x, y = weights, color = "true_weights"), linewidth = 0.6) +
  scale_color_manual(name = "Type of weights", values = c('weights_estimation' = 'red', 'true_weights' = 'blue')) +
  labs(title = "True Weights and Weights Estimations",
       x = "Sampled x",
       y = "Weights Hat") +
  theme_bw()+
  theme(legend.position = "bottom")

print(p)

estimated_probs_NR_cv <- get_dens_or_prob(centered_kernel_mat_at_sampled,
                                          centered_kernel_mat_at_grid,
                                          centerd_kernel_self_grid,
                                          sampled_x_vec,x_grid,
                                          lambda_hat, as.vector(weights_hat),
                                          type_of_p_is_prob = FALSE,
                                          type_of_q_is_prob = FALSE,
                                          method_of_p_calculation = "ordinary")



estimated_density_df <- data.frame(grid = x_grid, estimated_pdf = estimated_probs_NR_cv$grid_x)

ggplot() +
  geom_line(data = true_density_df, aes(x = grid, y = true_pdf, color = 'p(x)'), linewidth = 1, alpha =1) +
  geom_line(data = true_n_density_df, aes(x = grid, y = n_pdf, color = 'p_n(x)'), linewidth = 1, alpha = 1) +
  geom_line(data = estimated_density_df, aes(x = grid, y = estimated_pdf, color = 'KEF'), linewidth = 1) +
  scale_color_manual(name = "Type of density", values = c('p(x)' = 'blue', 'p_n(x)' = 'chartreuse3', 'KEF' = 'orange')) +
  ggtitle(paste('Histogram and Kernel Density Estimate for lambda =', lambda_hat, 'and tau =', format(tau_hat,digits = 3,scientific = T))) +
  xlab('u') +
  ylab('Density') +
  theme_bw() +
  theme(legend.position = "bottom")

################################################################################

s1_new <- colSums(centered_kernel_mat_at_sampled)
summary(s1_new)

# Find the base measure of samples
sample_mid_points <- get_middle_points_grid(-3.1, sampled_x_vec, 3.1)
base_measure_weights <- sample_mid_points[-1] - sample_mid_points[-length(sample_mid_points)]

dens_n <- probs_n_cv$sampled_x
dens_sampled_base_n <- dens_n * base_measure_weights
prob_sampled_base_n <- dens_sampled_base_n / sum(dens_sampled_base_n)
prob_sampled_n <- dens_n / sum(dens_n)

s2_new <-  as.numeric(dim(centered_kernel_mat_at_sampled)[1] * (prob_sampled_base_n) %*% t(centered_kernel_mat_at_sampled))
summary(s2_new)

s12_new<- colSums(centered_kernel_mat_at_sampled) - dim(centered_kernel_mat_at_sampled)[1] * (prob_sampled_base_n) %*% t(centered_kernel_mat_at_sampled)

s12t_new <- as.numeric(s12_new)
summary(as.numeric(s12_new))

plot(x = sampled_x_vec,s1_new)
plot(x = sampled_x_vec,as.numeric(s2_new))
plot(x = sampled_x_vec,s12t_new)

s3_new <- as.vector(w_x_vec) / prob_sampled_n
summary(s3_new)

plot(x = sampled_x_vec, s3_new)
summary(s12t_new)

summary(s3_new)/summary(s12t_new)

plot(x = sampled_x_vec,s12t_new * lambda_hat - s3_new*tau_hat)

################################################################################
output <- vector("numeric",length = length(sampled_x_vec))
# Find the base measure of samples
sample_mid_points <- get_middle_points_grid(-3.1, sampled_x_vec, 3.1)
base_measure_weights <- sample_mid_points[-1] - sample_mid_points[-length(sample_mid_points)]

n_x <- length(sampled_x_vec)
dens <- probs_n_cv$sampled_x
dens_sampled_base <- dens * base_measure_weights
prob_sampled_base <- dens_sampled_base / sum(dens_sampled_base)
prob_sampled <- dens / sum(dens)


for( i in 1: length(sampled_x_vec)){
  output[i] <- w_x_vec[i]/
    (probs_n_cv$sampled_x[i] *
       (colSums(centered_kernel_mat_at_sampled)[i] - n_x * (prob_sampled_base %*% centered_kernel_mat_at_sampled)[1,i] ))
}
which(abs(output)>10)
