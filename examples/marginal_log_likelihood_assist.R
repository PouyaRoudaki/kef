library(quantreg)
library(ggplot2)
library(plotly)
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
x_grid <-  rep(0,400)
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

  tau = 1/135
  # Get the number of sampled points
  n <- nrow(centered_kernel_mat_at_sampled)


  plot(sampled_x,p_vec)

  # Sample weights w_i from a normal distribution N(0, p(x_i)/tau)
  # Multiply each column by corresponding vector value
  w_sampled <- sweep(std_rnorm_matrix, 2,  p_vec/tau, `*`)
  #w_sampled <- adjusted_mc_rnorm(std_rnorm,n = n,new_var_vec = p_vec/tau)

  hist(w_sampled[,1])
  data <- as.data.frame(cbind(rep(1:length(sampled_x), each = MC_iterations),
                              as.numeric(w_sampled)))

  data <- data %>% mutate(sampled_x = sampled_x[V1])

  data1 <- data %>%
    group_by(V1) %>%
    summarise(var(V2))

  print(plot(sampled_x, sqrt(data1$`var(V2)`)))

  print(plot(data$sampled_x,data$V2))
  # Calculate the probability for each set of sampled weights using the custom 'get_dens_wo_grid' function
  densities_for_given_weights <- t(apply(w_sampled, 1, function(w_vector) {

    #plot(sampled_x,unique(sign(w_vector))*w_vector)
    get_dens_wo_grid(centered_kernel_mat_at_sampled,
                     min_x,
                     max_x,
                     sampled_x,
                     lambda,
                     w_vector)
  }))

  data <- as.data.frame(cbind(rep(1:length(sampled_x), each = MC_iterations),
                              as.numeric(densities_for_given_weights)))

  data <- data %>% mutate(sampled_x = sampled_x[V1])

  data1 <- data %>%
    group_by(V1) %>%
    summarise(mean(V2))

  print(plot(sampled_x, sqrt(data1$`mean(V2)`)))

  #print(plot(data$sampled_x,data$V2))

  #print(dim(probabilities_for_given_weights))

  #sample_mid_points <- get_middle_points_grid(min_x, sampled_x, max_x)
  #base_measure_weights <- sample_mid_points[-1] - sample_mid_points[-length(sample_mid_points)]

  #prob_inverse_base <- probabilities_for_given_weights * exp(-base_measure_weights)

  dim(densities_for_given_weights)
  # Compute the likelihood vector by taking the mean across Monte Carlo iterations
  likelihood_vector <- colMeans(densities_for_given_weights, na.rm = TRUE)

  #print(length(likelihood_vector))
  print(plot(sampled_x,likelihood_vector))

  # Create a data frame with your data
  plot_data <- data.frame(x = sampled_x, likelihood = likelihood_vector)

  # Plot using ggplot
  plot <- ggplot(plot_data, aes(x = x, y = likelihood)) +
    geom_line() +  # Line plot
    geom_point() +  # Add points to emphasize sampled values
    labs(title = paste("Densities vs. Sampled X,", expression(lambda),"=",lambda
                       ,",", expression(tau),"=",tau),
         x = "Sampled X",
         y = "Densities") +
    theme_bw()  # Use a clean theme

  print(plot)

  # Compute the log of the marginal likelihood by summing the log of the likelihood vector
  marginal_log_likelihood <- sum(log(likelihood_vector))

  # Print a separator (seems to be for debugging purposes)
  #print("-")

  # Return the computed marginal log likelihood
  return(marginal_log_likelihood)
}

#' Title: Compute Marginal Likelihood over a Grid of Parameters
#'
#' This function computes the marginal likelihood for each combination of lambda and tau
#' in the provided grids.
#'
#' @param centered_kernel_mat_at_sampled A matrix of centered kernel values at sampled points.
#' @param centered_kernel_mat_at_grid A matrix of centered kernel values at grid points.
#' @param centered_kernel_self_grid A matrix of centered kernel self-values at grid points.
#' @param x_grid A grid of x values where the function is evaluated.
#' @param p_vec A vector of initial probabilities (default is a vector of 1's).
#' @param lambda_grid A grid of lambda values to be evaluated.
#' @param tau_grid A grid of tau values to be evaluated.
#' @param MC_iterations Number of Monte Carlo iterations to perform.
#'
#' @return A data frame containing lambda, tau, and their corresponding marginal log likelihoods.
#' @export
#'
compute_marginal_likelihood_grid <- function(centered_kernel_mat_at_sampled,
                                             min_x,
                                             max_x,
                                             sampled_x,
                                             lambda_grid,
                                             tau_grid,
                                             initial_lambda = 1,
                                             initial_w = rep(0,length(sampled_x)),
                                             MC_iterations) {

  min_x=-3.1
  max_x=3.1

  # Compute the marginal likelihood for each combination in the grid using mapply
  t <- 1

  n <- length(sampled_x)

  MC_iterations<-1000
  # Generate matrix with each row independently sampled from Uniform(0,1)
  std_rnorm_matrix <- matrix(rnorm(MC_iterations * n, mean = 0, sd = 1),
                           nrow = MC_iterations,
                           ncol = n)

  hist(std_rnorm_matrix[,1])

  lambda <- 1
  w_vec <- weights_hat_wo_grid
  dens_vec <- get_dens_wo_grid(centered_kernel_mat_at_sampled,
                               min_x,
                               max_x,
                               sampled_x,
                               lambda,
                               w_vec)

  p_vec <- dens_vec/sum(dens_vec)
  plot(sampled_x,p_vec)

  while (t <= 10) {

    plot_ml <- data.frame(
      lambda = numeric(),
      tau = numeric(),
      mll = numeric()
    )


    max_marginal_log_likelihood <- -Inf

    for (lambda in lambda_grid) {
      for (tau in tau_grid) {
        cat(paste0("-"))


        #print(length(p_vec))
        marginal_log_likelihood <- marginal_log_likelihood(
          centered_kernel_mat_at_sampled,
          sampled_x,
          min_x,
          max_x,
          p_vec,
          lambda,
          tau,
          std_rnorm,
          MC_iterations)

        #print(lambda_hat)
        #print(marginal_log_likelihood)
        #print(max_marginal_log_likelihood)

        #print(paste(lambda_hat,",",marginal_log_likelihood))

        plot_ml <- rbind(plot_ml, data.frame(
          lambda = lambda,
          tau = tau,
          mll = marginal_log_likelihood
        ))

        #print(plot_ml)

        if(marginal_log_likelihood > max_marginal_log_likelihood){
          max_marginal_log_likelihood <- marginal_log_likelihood
          max_likelihood_lambda <- lambda
          max_likelihood_tau <- tau
        }
      }
    }

    lambda <- max_likelihood_lambda
    tau <- max_likelihood_tau

    w_vec <- get_weights_wo_grid_mll(lambda_t = lambda,
                                     tau_t = tau,
                                     centered_kernel_mat_at_sampled = centered_kernel_mat_at_sampled,
                                     sampled_x = sampled_x,
                                     min_x = min_x,
                                     max_x = max_x,
                                     p_vec_t_1 = p_vec,
                                     print_trace =  F
    )

    dens_vec <- get_dens_wo_grid(centered_kernel_mat_at_sampled,
                                 min_x,
                                 max_x,
                                 sampled_x,
                                 lambda,
                                 w_vec)

    p_vec <- dens_vec/sum(dens_vec)

    # Add an extra newline after the loop
    cat("\n")
    print(paste0("Iteration = ", t,
                 ", lambda_hat = ", lambda,
                 ", tau_hat = ", tau,
                 ", max_marginal_log_likelihood = ",max_marginal_log_likelihood))

    library(ggplot2)
    library(viridis)  # For scale_fill_viridis_c()
    # Continuous heatmap plot using geom_tile
    p_ml <- ggplot(plot_ml, aes(x = log10(lambda), y = log10(tau), fill = mll)) +
      geom_tile() +
      scale_fill_viridis_c(option = "D", direction = -1) +  # Continuous color scale
      theme_minimal() +
      labs(
        title = paste0("Heatmap of Marginal Log Likelihood, Iteration: ", t),
        x = expression(log[10](lambda)),
        y = expression(log[10](tau)),
        fill = "Marginal Log-Likelihood"
      ) +
      theme(
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, size = 14)
      )

    # Print the plot
    print(p_ml)

    t <- t +1
  }
}



for (j in c(1)) {
  worst <- t(rbind(rank(-densities_for_given_weights[,j])[(rank(-densities_for_given_weights[,j])<6)],
                   which(rank(-densities_for_given_weights[,j])<6)))

  worst_df <- as.data.frame(worst)
  worst_df <- worst_df %>% arrange(V1)
  worst_df

  for (i in 1:5) {
    # Prepare data frame for plotting
    df <- data.frame(
      sampled_x = sampled_x,
      w_values = w_sampled[worst_df[i, 2], ],
      cen_brownian_kernel = centered_kernel_mat_at_sampled[j, ]
    )

    # Compute scaling factor for dual axis transformation
    scale_factor <- max(df$w_values) / max(df$cen_brownian_kernel)

    # Plot with dual y-axis for better visualization
    plot <- ggplot(df, aes(x = sampled_x)) +
      # First y-axis: w_values
      geom_point(aes(y = w_values, color = "w_values"), size = 2) +
      geom_smooth(aes(y = w_values), method = "loess", color = "red", se = FALSE) +

      # Second y-axis: cen_brownian_kernel transformed for scaling
      geom_point(aes(y = cen_brownian_kernel * scale_factor, color = "cen_brownian_kernel"),
                 size = 2, alpha = 0.5) +

      # Manual color mapping
      scale_color_manual(values = c("w_values" = "blue", "cen_brownian_kernel" = "springgreen3")) +

      # Configure dual y-axis
      scale_y_continuous(
        name = "w_sampled",
        sec.axis = sec_axis(~ . / scale_factor, name = "centered brownian kernel")
      ) +

      # Add plot labels and theme
      labs(
        title = paste0("Sample ",j ,", Worst densities: number ",i , ", ", worst_df[i, 2],
                       "th iteration of MC: w.h = ",
                       round(w_sampled[worst_df[i, 2], ] %*% centered_kernel_mat_at_sampled[,j],2)),
        x = "sampled_x",
        color = "Legend"
      ) +
      theme_bw()

    # Print plot
    print(plot)
  }

  best <- t(rbind(rank(densities_for_given_weights[,j])[(rank(densities_for_given_weights[,j])<6)],
                  which(rank(densities_for_given_weights[,j])<6)))

  best_df <- as.data.frame(best)
  best_df <- best_df %>% arrange(V1)
  best_df

  for (i in 1:5) {
    # Prepare data frame for plotting
    df <- data.frame(
      sampled_x = sampled_x,
      w_values = w_sampled[best_df[i, 2], ],
      cen_brownian_kernel = centered_kernel_mat_at_sampled[j, ]
    )

    # Compute scaling factor for dual axis transformation
    scale_factor <- max(df$w_values) / max(df$cen_brownian_kernel)

    # Plot with dual y-axis for better visualization
    plot <- ggplot(df, aes(x = sampled_x)) +
      # First y-axis: w_values
      geom_point(aes(y = w_values, color = "w_values"), size = 2) +
      geom_smooth(aes(y = w_values), method = "loess", color = "red", se = FALSE) +

      # Second y-axis: cen_brownian_kernel transformed for scaling
      geom_point(aes(y = cen_brownian_kernel * scale_factor, color = "cen_brownian_kernel"),
                 size = 2, alpha = 0.5) +

      # Manual color mapping
      scale_color_manual(values = c("w_values" = "blue", "cen_brownian_kernel" = "springgreen3")) +

      # Configure dual y-axis
      scale_y_continuous(
        name = "w_sampled",
        sec.axis = sec_axis(~ . / scale_factor, name = "centered brownian kernel")
      ) +

      # Add plot labels and theme
      labs(
        title = paste0("Sample ",j,", Best densities: number ",i , ", ", best_df[i, 2],
                       "th iteration of MC: w.h = ",
                       round(w_sampled[best_df[i, 2], ] %*% centered_kernel_mat_at_sampled[,j],2)),
        x = "sampled_x",
        color = "Legend"
      ) +
      theme_bw()

    # Print plot
    print(plot)
  }
}


