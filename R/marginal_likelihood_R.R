#' Title: Compute Marginal Log Likelihood
#'
#' This function calculates the marginal likelihood for a given set of parameters.
#'
#' @param centered_kernel_mat_at_sampled A matrix of centered kernel values at sampled points.
#' @param centered_kernel_mat_at_grid A matrix of centered kernel values at grid points.
#' @param centered_kernel_self_grid A matrix of centered kernel self-values at grid points.
#' @param x_grid A grid of x values where the function is evaluated.
#' @param dens_vec A vector of initial probabilities (default is a vector of 1's).
#' @param lambda A scalar parameter, lambda, to be used in the likelihood computation.
#' @param tau A scalar parameter, tau, to be used in the likelihood computation.
#' @param MC_iterations Number of Monte Carlo iterations to perform.
#' @param censoring True if the margnials is censored.
#'
#' @return The log of the marginal likelihood.
#' @export
#'
marginal_log_likelihood_R <- function(centered_kernel_mat_at_sampled,
                                sampled_x,
                                min_x,
                                max_x,
                                p_vec = rep(1,nrow(centered_kernel_mat_at_sampled)),
                                lambda,
                                tau,
                                std_rnorm_matrix,
                                MC_iterations,
                                censoring = FALSE){

  # Get the number of sampled points
  n <- nrow(centered_kernel_mat_at_sampled)



  # Sample weights w_i from a normal distribution N(0, p(x_i)/tau)
  w_sampled <- sweep(std_rnorm_matrix, 2,  sqrt(p_vec/tau), `*`)

  #data <- as.data.frame(cbind(rep(1:length(sampled_x), MC_iterations),
  #                            as.numeric(w_sampled)))

  #data <- data %>% mutate(sampled_x = sampled_x[V1])

  #data1 <- data %>%
  #  group_by(V1) %>%
  #  summarise(var(V2))

  #print(plot(sampled_x, sqrt(data1$`var(V2)`)))

  #print(plot(data$sampled_x,data$V2))
  # Calculate the probability for each set of sampled weights using the custom 'get_dens_wo_grid' function
  densities_for_given_weights <- apply(w_sampled, 1, function(w_vector) {
    as.vector(get_dens_wo_grid(centered_kernel_mat_at_sampled,
                     min_x,
                     max_x,
                     sampled_x,
                     lambda,
                     w_vector))
  })

  #print(dim(probabilities_for_given_weights))

  # Assuming your matrix is called densities_for_given_weights this will calculate
  # the marginal joint conditional density for a simulation of Monte Carlo
  monte_carlo_likelihood <- apply(densities_for_given_weights, 2, prod)

  #print(dim(monte_carlo_likelihood))
  #print(length(monte_carlo_likelihood))
  # Find the marginal likelihood: average of all the marginal joint conditional
  # densities for all Monte Carlo simulations.
  marginal_likelihood <- mean(monte_carlo_likelihood)

  # Compute the likelihood vector by taking the mean across Monte Carlo iterations
  #likelihood_vector <- rowMeans(densities_for_given_weights, na.rm = F)

  #sample_mid_points <- get_middle_points_grid(min_x, sampled_x, max_x)
  #base_measure_weights <- sample_mid_points[-1] - sample_mid_points[-length(sample_mid_points)]

  # Define number of neighbors
  #k <- sqrt(100)

  # Find nearest neighbors
  #neighbors <- get.knn(as.matrix(sampled_x), k = k)

  # Compute local standard deviation
  #local_sd <- sapply(1:length(sampled_x), function(i) {
  #  neighbor_indices <- neighbors$nn.index[i, ]
  #  neighbor_values <- sampled_x[neighbor_indices]
  #  sd(neighbor_values)
  #})

  #likelihood_vector <- likelihood_vector/exp(((max(sampled_x) - min(sampled_x))*local_sd))

  #weight_average <- colMeans(w_sampled, na.rm = TRUE)

  #likelihood_vector <- as.vector(get_dens_wo_grid(centered_kernel_mat_at_sampled,
  #                             min_x,
  #                             max_x,
  #                             sampled_x,
  #                             lambda,
  #                             weight_average))

  # Compute the likelihood vector by taking the median across Monte Carlo iterations
  #likelihood_vector <- apply(densities_for_given_weights, 1, median, na.rm = TRUE)

  #if(censoring){
  #  likelihood_vector_trimmed <- likelihood_vector[(5):(length(likelihood_vector)-4)]
  #  sampled_x_trimmed <- sampled_x[(5):(length(sampled_x)-4)]
  #  # Normalize the density by the integral over the grid
  #  normalizing_cte <- pracma::trapz( sampled_x_trimmed, likelihood_vector_trimmed)  # trapz is from pracma package

  #  norm_likelihood_vector <- likelihood_vector_trimmed/normalizing_cte

  #  sampled_x <- sampled_x_trimmed
  #}else{
  #  normalizing_cte <- pracma::trapz( sampled_x, likelihood_vector)  # trapz is from pracma package

  #  norm_likelihood_vector <- likelihood_vector/normalizing_cte
  #}


  # Initial guess (avoid zeros)
  #z0 <- rep(1, length(sampled_x))

  # Define the function in log-space
  #f <- function(z) {
  #  x <- exp(z)  # Ensure positivity
  #  S <- sum(base_measure_weights * x)
  #  return(base_measure_weights * x / S - norm_likelihood_vector)
  #}

  # Solve using Broyden's method
  #result <- nleqslv(z0, f, method = "Broyden")

  #norm_likelihood_vector_base <- exp(result$x)
  #norm_likelihood_vector_base <- norm_likelihood_vector_base/sum(norm_likelihood_vector_base)

  #norm_likelihood_vector <- norm_likelihood_vector_new
  #norm_likelihood_vector <- likelihood_vector_me
  #normalizing_cte <- pracma::trapz( sampled_x, likelihood_vector)  # trapz is from pracma package

  #norm_likelihood_vector <- likelihood_vector/normalizing_cte

  #likelihood_vector_std <- norm_likelihood_vector/sum(norm_likelihood_vector)
  #sample_mid_points <- get_middle_points_grid(min_x, sampled_x, max_x)
  #base_measure_weights <- sample_mid_points[-1] - sample_mid_points[-length(sample_mid_points)]
  #likelihood_with_base <- likelihood_vector * base_measure_weights

  #print(length(likelihood_vector))
  #print(plot(sampled_x,likelihood_vector))


  #norm_likelihood_vector <- abs(norm_likelihood_vector - base_measure_weights)

  # Create a data frame with your data
  #plot_data <- data.frame(x = sampled_x, likelihood = norm_likelihood_vector)

  #log_lambda <- log10(lambda)
  #log_tau <- log10(tau)
  # Plot using ggplot
  #plot <- ggplot(plot_data, aes(x = x, y = likelihood)) +
  #  geom_line(color = "black", size = 1) +  # Line plot with color and size
  #  geom_point(color = "black", size = 2) +  # Add points with color and size to emphasize sampled values
  #  labs(title = bquote("Likelihood vs. Sampled X, "
  #                      #~ log[10](lambda) == .(format(log_lambda, digits = 3, scientific = TRUE)) ~ ","
  #                         ~ log[10](tau) == .(format(log_tau, digits = 3, scientific = TRUE)) ~
  #                        "," ~ log-likelihood == .(format(sum(log(norm_likelihood_vector)), digits = 3, scientific = TRUE))),
  #       x = "Sampled X",
  #       y = "Likelihood") +
  #  theme_bw() +  # Use a clean theme
  #  theme(plot.title = element_text(hjust = 0.5))  # Center the title

  #print(plot)

  # Compute the log of the marginal likelihood by summing the log of the likelihood vector
  marginal_log_likelihood <- log(marginal_likelihood)

  # Print a separator (seems to be for debugging purposes)
  #print("-")

  # Return the computed marginal log likelihood
  return(marginal_log_likelihood)
}

#' Title: Compute Marginal Likelihood over a Grid of Parameters Parallelly.
#'
#' This function computes the marginal likelihood for each combination of lambda and tau
#' in the provided grids.
#'
#' @param centered_kernel_mat_at_sampled A matrix of centered kernel values at sampled points.
#' @param centered_kernel_mat_at_grid A matrix of centered kernel values at grid points.
#' @param centered_kernel_self_grid A matrix of centered kernel self-values at grid points.
#' @param x_grid A grid of x values where the function is evaluated.
#' @param p_vec A vector of initial probabilities (default is a vector of 1's).
#' @param hyperparam_grid A dataframe including pairs of lambda and tau values.
#' @param MC_iterations Number of Monte Carlo iterations to perform.
#' @param censoring True if the margnials is censored.
#'
#' @return A data frame containing lambda, tau, and their corresponding marginal log likelihoods.
#' @export
#'

compute_marginal_likelihood_grid_parallel_R <- function(centered_kernel_mat_at_sampled,
                                             min_x,
                                             max_x,
                                             sampled_x,
                                             hyperparam_grid,
                                             initial_lambda = 1,
                                             initial_w = rep(0, length(sampled_x)),
                                             MC_iterations,
                                             max.iterations = 1,
                                             censoring = F) {

  t <- 1
  n <- length(sampled_x)

  std_rnorm_matrix <- matrix(stats::rnorm(MC_iterations * n, mean = 0, sd = 1),
                             nrow = MC_iterations,
                             ncol = n)

  lambda <- initial_lambda
  w_vec <- initial_w
  dens_vec <- as.numeric(get_dens_wo_grid(centered_kernel_mat_at_sampled,
                               min_x,
                               max_x,
                               sampled_x,
                               lambda,
                               w_vec))

  p_vec <- dens_vec / sum(dens_vec)

  doParallel::registerDoParallel(cores = parallel::detectCores() - 1)

  while (t <= max.iterations) {
    plot_ml <- data.frame(lambda = numeric(), tau = numeric(), mll = numeric())
    max_marginal_log_likelihood <- -Inf

    results <- foreach::foreach(i = 1:nrow(hyperparam_grid), .combine = 'rbind',
                                .packages = c("dplyr","ggplot2","plotly")) %dopar% {

      lambda <- hyperparam_grid[i,1]
      tau <- hyperparam_grid[i,2]

      marginal_log_likelihood <- marginal_log_likelihood_R(
        centered_kernel_mat_at_sampled,
        sampled_x,
        min_x,
        max_x,
        p_vec,
        lambda,
        tau,
        std_rnorm_matrix,
        MC_iterations,
        censoring
      )

        data.frame(lambda = lambda, tau = tau, mll = marginal_log_likelihood)
    }

    plot_ml <- results

    best <- plot_ml[which.max(plot_ml$mll), ]
    lambda <- best$lambda
    tau <- best$tau
    max_marginal_log_likelihood <- best$mll



    cat("\nIteration:", t, ",lambda_hat:", lambda, ",tau_hat:", tau, ",max_mll:", round(max_marginal_log_likelihood, 2), "\n")

    #cat("dim of p_vec:", dim(p_vec), "length_p_vec", length(p_vec), "dim of kernel_matrix", dim(centered_kernel_mat_at_sampled), "length samples", length(sampled_x))

    w_vec <- get_weights_wo_grid_mll(lambda, tau, centered_kernel_mat_at_sampled, sampled_x, min_x, max_x, p_vec, FALSE)
    dens_vec <- as.numeric(get_dens_wo_grid(centered_kernel_mat_at_sampled, min_x, max_x, sampled_x, lambda, w_vec))
    p_vec <- dens_vec / sum(dens_vec)



    plot_ml <- plot_ml %>% filter( mll > -200)
    #cat(min(log10(plot_ml$lambda)))
    #cat(max(log10(plot_ml$lambda)))
    #x_grid <- seq(min(log10(plot_ml$lambda)), max(log10(plot_ml$lambda)), length.out = 100)
    #y_grid <- seq(min(log10(plot_ml$tau)), max(log10(plot_ml$tau)), length.out = 100)

    #interp_data <- with(plot_ml, akima::interp(log10(lambda), log10(tau), mll, xo = x_grid, yo = y_grid, duplicate = "mean"))
    #interp_df <- expand.grid(x = interp_data$x, y = interp_data$y)
    #interp_df$z <- as.vector(interp_data$z)

    #p_ml <- plotly::plot_ly() %>%
    #  plotly::add_surface(x = interp_data$x, y = interp_data$y, z = matrix(interp_df$z, nrow = length(interp_data$x), ncol = length(interp_data$y), byrow = TRUE), colorscale = 'Viridis') %>%
    #  plotly::layout(title = paste0("Marginal Log Likelihood-Iteration:", t, ", Censored:",censoring),
    #                 scene = list(xaxis = list(title = "log(λ)"), yaxis = list(title = "log(τ)"), zaxis = list(title = "MLL")))

    #print(p_ml)

    #scatter <- plotly::plot_ly(plot_ml, x = ~log10(lambda), y = ~log10(tau), z = ~mll, type = 'scatter3d', mode = 'markers',
    #                           marker = list(size = 4, color = ~mll, colorscale = 'Viridis')) %>%
    #  plotly::layout(title = paste0("Scatter plot Marginal Log Likelihood-Iteration:", t),
    #                 scene = list(xaxis = list(title = "log(λ)"), yaxis = list(title = "log(τ)"), zaxis = list(title = "MLL")))

    #print(scatter)

    t <- t + 1
  }

  doParallel::stopImplicitCluster()

  return(plot_ml)
}

#' Title: Compute Marginal Likelihood over a Grid of Parameters

#' This function computes the marginal likelihood for each combination of lambda and tau
#' in the provided grids.

#' @param centered_kernel_mat_at_sampled A matrix of centered kernel values at sampled points.
#' @param centered_kernel_mat_at_grid A matrix of centered kernel values at grid points.
#' @param centered_kernel_self_grid A matrix of centered kernel self-values at grid points.
#' @param x_grid A grid of x values where the function is evaluated.
#' @param p_vec A vector of initial probabilities (default is a vector of 1's).
#' @param hyperparam_grid A dataframe including pairs of lambda and tau values.
#' @param MC_iterations Number of Monte Carlo iterations to perform.
#' @param max.iterations Maximum number of iterations (default: 5).
#' @param seed An integer that controls the randomness.

#' @return A data frame containing lambda, tau, and their corresponding marginal log likelihoods.
#' @export

compute_marginal_likelihood_grid_R <- function(centered_kernel_mat_at_sampled,
                                               min_x,
                                               max_x,
                                               sampled_x,
                                               hyperparam_grid,
                                               initial_lambda = 1,
                                               initial_w = rep(0, length(sampled_x)),
                                               MC_iterations,
                                               max.iterations = 5,
                                               seed = 1) {

  t <- 1
  n <- length(sampled_x)

  # Generate matrix with each row independently sampled from Uniform(0,1)
  set.seed(seed)
  std_rnorm_matrix <- matrix(rnorm(MC_iterations * n, mean = 0, sd = 1),
                             nrow = MC_iterations,
                             ncol = n)

  lambda <- initial_lambda
  w_vec <- initial_w
  dens_vec <- as.numeric(get_dens_wo_grid(centered_kernel_mat_at_sampled,
                               min_x,
                               max_x,
                               sampled_x,
                               lambda,
                               w_vec))

  p_vec <- dens_vec / sum(dens_vec)

  p_vec_init <- p_vec

  while (t <= max.iterations) {

    plot_ml <- data.frame(
      lambda = numeric(),
      tau = numeric(),
      mll = numeric()
    )

    max_marginal_log_likelihood <- -Inf
    #max_likelihood_lambda <- -Inf
    #max_likelihood_tau <- -Inf

    for (i in (1:dim(hyperparam_grid)[1])) {
      lambda <- hyperparam_grid[i, 1]
      tau <- hyperparam_grid[i, 2]

      cat(paste0("-"))

      marginal_log_likelihood <- marginal_log_likelihood(centered_kernel_mat_at_sampled,
                                                         sampled_x,
                                                         min_x,
                                                         max_x,
                                                         p_vec = p_vec,
                                                         lambda,
                                                         tau,
                                                         std_rnorm_matrix,
                                                         MC_iterations,
                                                         parallel_computing = TRUE)


      plot_ml <- rbind(plot_ml, data.frame(
        lambda = lambda,
        tau = tau,
        mll = marginal_log_likelihood
      ))

      if (marginal_log_likelihood > max_marginal_log_likelihood & (!is.nan(marginal_log_likelihood))) {
        max_marginal_log_likelihood <- marginal_log_likelihood
        max_likelihood_lambda <- lambda
        max_likelihood_tau <- tau
      }
    }

    lambda <- max_likelihood_lambda
    tau <- max_likelihood_tau

    # Update weights
    w_vec <- get_weights_wo_grid_BBsolve(lambda_hat = lambda,
                                           tau_hat = tau,
                                           centered_kernel_mat_at_sampled = centered_kernel_mat_at_sampled,
                                           sampled_x = sampled_x,
                                           min_x = min_x,
                                           max_x = max_x,
                                           prior_variance_p_vector = p_vec,
                                           print_trace = FALSE)

    dens_vec <- as.numeric(get_dens_wo_grid(centered_kernel_mat_at_sampled,
                                 min_x,
                                 max_x,
                                 sampled_x,
                                 lambda,
                                 w_vec))

    p_vec <- dens_vec / sum(dens_vec)

    #output_file <- "marginal_likelihood_results.txt"

    cat(paste0("Iteration = ", t,
               ", lambda_hat = ", lambda,
               ", tau_hat = ", tau,
               ", max_marginal_log_likelihood = ", round(max_marginal_log_likelihood, 2),
               ", The ratio: ", lambda^2/tau ,"\n")
    )
        #file = output_file, append = TRUE)

    plot_ml <- plot_ml %>% filter(mll > -140)

    #plot <- ggplot(plot_ml, aes(x = log10(tau), y = mll)) +
    #  geom_point(color = "red") +
    #  geom_line(color = "red") +
    #  geom_vline(xintercept = log10(1/1350), linetype = "dotted", color = "blue", size = 1) + # Dotted vertical line
    #  labs(
    #    title = bquote("Likelihood vs." ~ log[10](tau) ~
    #                     "," ~ log[10](lambda) == .(0)),
    #    x = bquote(~ log[10](tau)),
    #    y = "Likelihood"
    #  ) +
    #  theme_bw()

    #print(plot)
    x_grid <- seq(min(log10(plot_ml$lambda)), max(log10(plot_ml$lambda)), length.out = 100)
    y_grid <- seq(min(log10(plot_ml$tau)), max(log10(plot_ml$tau)), length.out = 100)

    interp_data <- with(plot_ml, akima::interp(log10(lambda), log10(tau), mll, xo = x_grid, yo = y_grid, duplicate = "mean"))
    interp_df <- expand.grid(x = interp_data$x, y = interp_data$y)
    interp_df$z <- as.vector(interp_data$z)

    p_ml <- plotly::plot_ly() %>%
      plotly::add_surface(x = interp_data$x, y = interp_data$y, z = matrix(interp_df$z, nrow = length(interp_data$x), ncol = length(interp_data$y), byrow = TRUE), colorscale = 'Viridis') %>%
      plotly::layout(title = paste0("Marginal Log Likelihood-Iteration:", t, ", Seed:",seed),
                     scene = list(xaxis = list(title = "log10(λ)"), yaxis = list(title = "log10(τ)"), zaxis = list(title = "MLL")))
    #print(p_ml)

    scatter <- plot_ly(
      plot_ml,
      x = ~log10(lambda), y = ~log10(tau), z = ~mll,
      type = 'scatter3d', mode = 'markers',
      marker = list(size = 4, color = ~mll, colorscale = 'Viridis')
    ) %>%
      layout(
        title = paste0("Scatter plot Marginal Log Likelihood-Iteration:", t, ", Seed:", seed),
        scene = list(
          xaxis = list(title = "log10(\u03bb)"),
          yaxis = list(title = "log10(\u03c4)"),
          zaxis = list(title = "MLL")
        )
      )

    print(scatter)

    t <- t + 1
  }

  lst <- list()
  lst[[1]] <- plot_ml
  lst[[2]] <- interp_df
  lst[[3]] <- std_rnorm_matrix
  lst[[4]] <- p_vec_init
  return(lst)
}


#' Title: Compute Marginal Likelihood over a Grid of Parameters

#' This function computes the marginal likelihood for each combination of lambda and tau
#' in the provided grids.

#' @param centered_kernel_mat_at_sampled A matrix of centered kernel values at sampled points.
#' @param centered_kernel_mat_at_grid A matrix of centered kernel values at grid points.
#' @param centered_kernel_self_grid A matrix of centered kernel self-values at grid points.
#' @param x_grid A grid of x values where the function is evaluated.
#' @param p_vec A vector of initial probabilities (default is a vector of 1's).
#' @param hyperparam_grid A dataframe including pairs of lambda and tau values.
#' @param MC_iterations Number of Monte Carlo iterations to perform.
#' @param max.iterations Maximum number of iterations (default: 5).
#' @param seed An integer that controls the randomness.

#' @return A data frame containing lambda, tau, and their corresponding marginal log likelihoods.
#' @export

compute_marginal_likelihood_grid_R_pint <- function(centered_kernel_mat_at_sampled,
                                               min_x,
                                               max_x,
                                               sampled_x,
                                               hyperparam_grid,
                                               initial_lambda = 1,
                                               initial_tau = 1/1350,
                                               MC_iterations,
                                               max.iterations = 5,
                                               seed = 1) {

  t <- 1
  n <- length(sampled_x)

  # Generate matrix with each row independently sampled from Uniform(0,1)
  set.seed(seed)
  std_rnorm_matrix <- matrix(rnorm(MC_iterations * n, mean = 0, sd = 1),
                             nrow = MC_iterations,
                             ncol = n)

  lambda <- initial_lambda
  tau <- initial_tau

  grid <- seq(from = min_x, to = max_x,length.out = 4*n)

  dens_vec <- kef(sampled_x,grid,lambda = lambda,tau = tau)$probs_sample

  p_vec <- dens_vec / sum(dens_vec)

  p_vec_init <- p_vec

  while (t <= max.iterations) {

    plot_ml <- data.frame(
      lambda = numeric(),
      tau = numeric(),
      mll = numeric()
    )

    max_marginal_log_likelihood <- -Inf
    #max_likelihood_lambda <- -Inf
    #max_likelihood_tau <- -Inf

    for (i in (1:dim(hyperparam_grid)[1])) {
      lambda <- hyperparam_grid[i, 1]
      tau <- hyperparam_grid[i, 2]

      cat(paste0("-"))

      marginal_log_likelihood <- marginal_log_likelihood(centered_kernel_mat_at_sampled,
                                                         sampled_x,
                                                         min_x,
                                                         max_x,
                                                         p_vec = p_vec,
                                                         lambda,
                                                         tau,
                                                         std_rnorm_matrix,
                                                         MC_iterations,
                                                         parallel_computing = TRUE)


      plot_ml <- rbind(plot_ml, data.frame(
        lambda = lambda,
        tau = tau,
        mll = marginal_log_likelihood
      ))

      if (marginal_log_likelihood > max_marginal_log_likelihood & (!is.nan(marginal_log_likelihood))) {
        max_marginal_log_likelihood <- marginal_log_likelihood
        max_likelihood_lambda <- lambda
        max_likelihood_tau <- tau
      }
    }

    lambda <- max_likelihood_lambda
    tau <- max_likelihood_tau

    # Update weights
    w_vec <- get_weights_wo_grid_BBsolve(lambda_hat = lambda,
                                         tau_hat = tau,
                                         centered_kernel_mat_at_sampled = centered_kernel_mat_at_sampled,
                                         sampled_x = sampled_x,
                                         min_x = min_x,
                                         max_x = max_x,
                                         prior_variance_p_vector = p_vec,
                                         print_trace = FALSE)

    dens_vec <- as.numeric(get_dens_wo_grid(centered_kernel_mat_at_sampled,
                                            min_x,
                                            max_x,
                                            sampled_x,
                                            lambda,
                                            w_vec))

    p_vec <- dens_vec / sum(dens_vec)

    #output_file <- "marginal_likelihood_results.txt"

    cat(paste0("Iteration = ", t,
               ", lambda_hat = ", lambda,
               ", tau_hat = ", tau,
               ", max_marginal_log_likelihood = ", round(max_marginal_log_likelihood, 2),
               ", The ratio: ", lambda^2/tau ,"\n")
    )
    #file = output_file, append = TRUE)

    plot_ml <- plot_ml %>% filter(mll > -140)

    #plot <- ggplot(plot_ml, aes(x = log10(tau), y = mll)) +
    #  geom_point(color = "red") +
    #  geom_line(color = "red") +
    #  geom_vline(xintercept = log10(1/1350), linetype = "dotted", color = "blue", size = 1) + # Dotted vertical line
    #  labs(
    #    title = bquote("Likelihood vs." ~ log[10](tau) ~
    #                     "," ~ log[10](lambda) == .(0)),
    #    x = bquote(~ log[10](tau)),
    #    y = "Likelihood"
    #  ) +
    #  theme_bw()

    #print(plot)
    x_grid <- seq(min(log10(plot_ml$lambda)), max(log10(plot_ml$lambda)), length.out = 100)
    y_grid <- seq(min(log10(plot_ml$tau)), max(log10(plot_ml$tau)), length.out = 100)

    interp_data <- with(plot_ml, akima::interp(log10(lambda), log10(tau), mll, xo = x_grid, yo = y_grid, duplicate = "mean"))
    interp_df <- expand.grid(x = interp_data$x, y = interp_data$y)
    interp_df$z <- as.vector(interp_data$z)

    p_ml <- plotly::plot_ly() %>%
      plotly::add_surface(x = interp_data$x, y = interp_data$y, z = matrix(interp_df$z, nrow = length(interp_data$x), ncol = length(interp_data$y), byrow = TRUE), colorscale = 'Viridis') %>%
      plotly::layout(title = paste0("Marginal Log Likelihood-Iteration p_int:", t, ", Seed:",seed),
                     scene = list(xaxis = list(title = "log10(λ)"), yaxis = list(title = "log10(τ)"), zaxis = list(title = "MLL")))
    print(p_ml)

    scatter <- plot_ly(
      plot_ml,
      x = ~log10(lambda), y = ~log10(tau), z = ~mll,
      type = 'scatter3d', mode = 'markers',
      marker = list(size = 4, color = ~mll, colorscale = 'Viridis')
    ) %>%
      layout(
        title = paste0("Scatter plot Marginal Log Likelihood-Iteration p_int:", t, ", Seed:", seed),
        scene = list(
          xaxis = list(title = "log10(\u03bb)"),
          yaxis = list(title = "log10(\u03c4)"),
          zaxis = list(title = "MLL")
        )
      )

    print(scatter)

    t <- t + 1
  }

  lst <- list()
  lst[[1]] <- plot_ml
  lst[[2]] <- interp_df
  lst[[3]] <- std_rnorm_matrix
  lst[[4]] <- p_vec_init
  return(lst)
}


#' Optimize Marginal Log-Likelihood using L-BFGS-B Optimization
#'
#' This function finds the optimal `lambda` and `tau` by maximizing the marginal log-likelihood
#' using an iterative approach. Instead of a grid search, it employs L-BFGS-B optimization
#' to efficiently search for the best parameters.
#'
#' @param centered_kernel_mat_at_sampled The kernel matrix centered at sampled points.
#' @param min_x Minimum value of the sampled domain.
#' @param max_x Maximum value of the sampled domain.
#' @param sampled_x Vector of sampled points.
#' @param initial_lambda Initial value for lambda (default: 1).
#' @param initial_tau Initial value for tau (default: 1).
#' @param initial_w Initial weights vector (default: zeros of sampled_x length).
#' @param MC_iterations Number of Monte Carlo iterations.
#' @param max.iterations Maximum number of iterations (default: 1).
#' @param censoring Boolean indicating whether censoring is applied (default: FALSE).
#'
#' @return A list containing:
#' 	\item{lambda}{Optimized value of lambda.}
#' 	\item{tau}{Optimized value of tau.}
#' 	\item{max_marginal_log_likelihood}{Maximum marginal log-likelihood value.}
#'
#' @examples
#' \dontrun{
#' result <- optimize_marginal_log_likelihood(
#'   centered_kernel_mat_at_sampled, min_x, max_x, sampled_x,
#'   MC_iterations = 1000, max.iterations = 5
#' )
#' print(result)
#' }
#'
#' @export
optimize_marginal_log_likelihood_R <- function(centered_kernel_mat_at_sampled,
                                             min_x,
                                             max_x,
                                             sampled_x,
                                             initial_lambda = 1,
                                             initial_w = rep(0, length(sampled_x)),
                                             MC_iterations,
                                             max.iterations = 1,
                                             censoring = FALSE) {

  t <- 1
  n <- length(sampled_x)

  # Generate matrix with each row independently sampled from Normal(0,1)
  std_rnorm_matrix <- matrix(rnorm(MC_iterations * n, mean = 0, sd = 1),
                             nrow = MC_iterations,
                             ncol = n)

  lambda <- initial_lambda
  w_vec <- initial_w

  dens_vec <- as.numeric(get_dens_wo_grid(centered_kernel_mat_at_sampled,
                                          min_x,
                                          max_x,
                                          sampled_x,
                                          lambda,
                                          w_vec))

  p_vec <- dens_vec / sum(dens_vec)

  while (t <= max.iterations) {

    cat(paste0("Iteration: ", t, "\n"))

    # Define objective function for optimization
    objective_function <- function(params) {
      log_lambda <- params[1]  # Optimizing log(lambda)
      theta <- params[2]  # Free parameter for tau reparameterization
      log_tau <- log_lambda - 9.671 + exp(theta)  # Enforcing constraint
      lambda <- exp(log_lambda)
      tau <- exp(log_tau)



      -marginal_log_likelihood_R(
        centered_kernel_mat_at_sampled,
        sampled_x,
        min_x,
        max_x,
        p_vec,
        lambda,
        tau,
        std_rnorm_matrix,
        MC_iterations,
        censoring)
    }

    cat(paste0("lambda: ", lambda,",tau: ", exp(log(lambda) - 9.671 + (.Machine$double.xmin)), "\n"))

    # Optimization using L-BFGS-B (bounded optimization)
    opt_result <- optim(
      par = c(log(lambda), log(.Machine$double.xmin) ),  # Initial values for log(lambda) and theta
      fn = objective_function,
      method = "L-BFGS-B",
      lower = c(log(1e-1), log(.Machine$double.xmin)),
      upper = c(log(1e1), 5)
    )

    # Retrieve optimal lambda and tau
    log_lambda <- opt_result$par[1]
    theta <- opt_result$par[2]
    log_tau <- log_lambda - 9.678 + exp(theta)
    lambda <- exp(log_lambda)
    tau <- exp(log_tau)

    cat(paste0("Optimized lambda: ", lambda, ", tau: ", tau, ", MLL: ", -opt_result$value,
               ", The ratio: ", lambda^2/tau ,"\n"))

    # Update weights
    w_vec <- get_weights_wo_grid_mll(lambda_t = lambda,
                                     tau_t = tau,
                                     centered_kernel_mat_at_sampled = centered_kernel_mat_at_sampled,
                                     sampled_x = sampled_x,
                                     min_x = min_x,
                                     max_x = max_x,
                                     p_vec_t_1 = p_vec,
                                     print_trace = FALSE)

    # Update density and p_vec
    dens_vec <- as.numeric(get_dens_wo_grid(centered_kernel_mat_at_sampled,
                                            min_x,
                                            max_x,
                                            sampled_x,
                                            lambda,
                                            w_vec))

    p_vec <- dens_vec / sum(dens_vec)

    t <- t + 1
  }

  return(list(lambda = lambda, tau = tau, max_marginal_log_likelihood = -opt_result$value))
}

