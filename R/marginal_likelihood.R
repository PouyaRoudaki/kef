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
marginal_log_likelihood <- function(centered_kernel_mat_at_sampled,
                                sampled_x,
                                min_x,
                                max_x,
                                p_vec = rep(1,nrow(centered_kernel_mat_at_sampled)),
                                lambda,
                                tau,
                                std_rnorm_matrix,
                                MC_iterations,
                                censoring = T){

  # Get the number of sampled points
  n <- nrow(centered_kernel_mat_at_sampled)



  # Sample weights w_i from a normal distribution N(0, p(x_i)/tau)
  w_sampled <- sweep(std_rnorm_matrix, 2,  p_vec/tau, `*`)

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
    get_dens_wo_grid(centered_kernel_mat_at_sampled,
                     min_x,
                     max_x,
                     sampled_x,
                     lambda,
                     w_vector)
  })

  #print(dim(probabilities_for_given_weights))



  # Compute the likelihood vector by taking the mean across Monte Carlo iterations
  likelihood_vector <- rowMeans(densities_for_given_weights, na.rm = TRUE)

  if(censoring){
    likelihood_vector_trimmed <- likelihood_vector[(5):(length(likelihood_vector)-4)]
    sampled_x_trimmed <- sampled_x[(5):(length(sampled_x)-4)]
    # Normalize the density by the integral over the grid
    normalizing_cte <- pracma::trapz( sampled_x_trimmed, likelihood_vector_trimmed)  # trapz is from pracma package

    norm_likelihood_vector <- likelihood_vector_trimmed/normalizing_cte

    sampled_x <- sampled_x_trimmed
  }else{
    normalizing_cte <- pracma::trapz( sampled_x, likelihood_vector)  # trapz is from pracma package

    norm_likelihood_vector <- likelihood_vector/normalizing_cte
  }

  #normalizing_cte <- pracma::trapz( sampled_x, likelihood_vector)  # trapz is from pracma package

  #norm_likelihood_vector <- likelihood_vector/normalizing_cte

  #likelihood_vector_std <- norm_likelihood_vector/sum(norm_likelihood_vector)
  #sample_mid_points <- get_middle_points_grid(min_x, sampled_x, max_x)
  #base_measure_weights <- sample_mid_points[-1] - sample_mid_points[-length(sample_mid_points)]
  #likelihood_with_base <- likelihood_vector * base_measure_weights

  #print(length(likelihood_vector))
  #print(plot(sampled_x,likelihood_vector))

  # Create a data frame with your data
  plot_data <- data.frame(x = sampled_x, likelihood = norm_likelihood_vector)

  log_lambda <- log10(lambda)
  log_tau <- log10(tau)
  # Plot using ggplot
  #plot <- ggplot(plot_data, aes(x = x, y = likelihood)) +
  #  geom_line() +  # Line plot
  #  geom_point() +  # Add points to emphasize sampled values
  #  labs(title = bquote("Densities vs. Sampled X, " ~ log[10](lambda) == .(format(log_lambda, digits = 3, scientific = TRUE)) ~
  #                        "," ~ log[10](tau) == .(format(log_tau, digits = 3, scientific = TRUE))),
  #       x = "Sampled X",
  #       y = "Densities") +
  #  theme_bw()  # Use a clean theme

  #print(plot)

  # Compute the log of the marginal likelihood by summing the log of the likelihood vector
  marginal_log_likelihood <- sum(log(norm_likelihood_vector))

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
#' @param censoring True if the margnials is censored.
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
                                             MC_iterations,
                                             max.iterations = 1,
                                             censoring = T) {


  # Compute the marginal likelihood for each combination in the grid using mapply
  t <- 1

  n <- length(sampled_x)

  # Generate matrix with each row independently sampled from Uniform(0,1)
  std_rnorm_matrix <- matrix(rnorm(MC_iterations * n, mean = 0, sd = 1),
                             nrow = MC_iterations,
                             ncol = n)

  lambda <- initial_lambda
  w_vec <- initial_w
  dens_vec <- get_dens_wo_grid(centered_kernel_mat_at_sampled,
                            min_x,
                            max_x,
                            sampled_x,
                            lambda_hat,
                            w_vec)

  p_vec <- dens_vec/sum(dens_vec)

  while (t <= max.iterations) {

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
          std_rnorm_matrix,
          MC_iterations,
          censoring)

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

        if(marginal_log_likelihood > max_marginal_log_likelihood & (!is.nan(marginal_log_likelihood))) {
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
                   ", max_marginal_log_likelihood = ",round(max_marginal_log_likelihood,2)))


    plot_ml <- plot_ml %>% filter(mll > -200)

    # Define a finer grid for smoother interpolation
    x_grid <- seq(min(log10(plot_ml$lambda)), max(log10(plot_ml$lambda)), length.out = 100)
    y_grid <- seq(min(log10(plot_ml$tau)), max(log10(plot_ml$tau)), length.out = 100)

    # Interpolate the scattered data onto the fine grid
    interp_data <- with(plot_ml, interp(
      x = log10(lambda),
      y = log10(tau),
      z = mll,
      xo = x_grid,
      yo = y_grid,
      duplicate = "mean"
    ))

    # Convert to structured format
    interp_df <- expand.grid(x = interp_data$x, y = interp_data$y)
    interp_df$z <- as.vector(interp_data$z)

    # Plot the smooth surface
    p_ml <- plot_ly() %>%
      add_surface(
        x = interp_data$x,
        y = interp_data$y,
        z = matrix(interp_df$z, nrow = length(interp_data$x), ncol = length(interp_data$y),byrow = T),
        colorscale = 'Viridis'
      ) %>%
      layout(
        title = paste0("Marginal Log Likelihood-Iteration:", t),
        scene = list(
          xaxis = list(title = "log(λ)"),
          yaxis = list(title = "log(τ)"),
          zaxis = list(title = "MLL")
        )
      )



    # Print the plot
    print(p_ml)

    scatter <- plot_ly(
      plot_ml,
      x = ~log10(lambda), y = ~log10(tau), z = ~mll,
      type = 'scatter3d', mode = 'markers',
      marker = list(size = 4, color = ~mll, colorscale = 'Viridis')
    ) %>%
      layout(
        title =  paste0("Scatter plot Marginal Log Likelihood-Iteration:", t),
        scene = list(
          xaxis = list(title = "log(λ)"),
          yaxis = list(title = "log(τ)"),
          zaxis = list(title = "MLL")
        )
      )
    print(scatter)

    t <- t +1
  }

  lst <- list()
  lst[[1]] <- plot_ml
  lst[[2]] <- interp_df
  return(lst)
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
#' @param lambda_grid A grid of lambda values to be evaluated.
#' @param tau_grid A grid of tau values to be evaluated.
#' @param MC_iterations Number of Monte Carlo iterations to perform.
#' @param censoring True if the margnials is censored.
#'
#' @return A data frame containing lambda, tau, and their corresponding marginal log likelihoods.
#' @export
#'

compute_marginal_likelihood_grid_parallel <- function(centered_kernel_mat_at_sampled,
                                             min_x,
                                             max_x,
                                             sampled_x,
                                             lambda_grid,
                                             tau_grid,
                                             initial_lambda = 1,
                                             initial_w = rep(0, length(sampled_x)),
                                             MC_iterations,
                                             max.iterations = 1,
                                             censoring = T) {

  t <- 1
  n <- length(sampled_x)

  std_rnorm_matrix <- matrix(stats::rnorm(MC_iterations * n, mean = 0, sd = 1),
                             nrow = MC_iterations,
                             ncol = n)

  lambda <- initial_lambda
  w_vec <- initial_w
  dens_vec <- get_dens_wo_grid(centered_kernel_mat_at_sampled,
                               min_x,
                               max_x,
                               sampled_x,
                               lambda,
                               w_vec)

  p_vec <- dens_vec / sum(dens_vec)

  doParallel::registerDoParallel(cores = parallel::detectCores() - 1)

  while (t <= max.iterations) {
    plot_ml <- data.frame(lambda = numeric(), tau = numeric(), mll = numeric())
    max_marginal_log_likelihood <- -Inf

    results <- foreach::foreach(lambda = lambda_grid, .combine = 'rbind', .packages = c("dplyr","ggplot2","plotly")) %:%
      foreach::foreach(tau = tau_grid, .combine = 'rbind') %dopar% {

        marginal_log_likelihood <- marginal_log_likelihood(
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

    w_vec <- get_weights_wo_grid_mll(lambda, tau, centered_kernel_mat_at_sampled, sampled_x, min_x, max_x, p_vec, FALSE)
    dens_vec <- get_dens_wo_grid(centered_kernel_mat_at_sampled, min_x, max_x, sampled_x, lambda, w_vec)
    p_vec <- dens_vec / sum(dens_vec)

    cat("\nIteration:", t, ",lambda_hat:", lambda, ",tau_hat:", tau, ",max_mll:", round(max_marginal_log_likelihood, 2), "\n")

    plot_ml <- plot_ml %>% filter( mll > -200)
    #cat(min(log10(plot_ml$lambda)))
    #cat(max(log10(plot_ml$lambda)))
    x_grid <- seq(min(log10(plot_ml$lambda)), max(log10(plot_ml$lambda)), length.out = 100)
    y_grid <- seq(min(log10(plot_ml$tau)), max(log10(plot_ml$tau)), length.out = 100)

    interp_data <- with(plot_ml, akima::interp(log10(lambda), log10(tau), mll, xo = x_grid, yo = y_grid, duplicate = "mean"))
    interp_df <- expand.grid(x = interp_data$x, y = interp_data$y)
    interp_df$z <- as.vector(interp_data$z)

    p_ml <- plotly::plot_ly() %>%
      plotly::add_surface(x = interp_data$x, y = interp_data$y, z = matrix(interp_df$z, nrow = length(interp_data$x), ncol = length(interp_data$y), byrow = TRUE), colorscale = 'Viridis') %>%
      plotly::layout(title = paste0("Marginal Log Likelihood-Iteration:", t),
                     scene = list(xaxis = list(title = "log(λ)"), yaxis = list(title = "log(τ)"), zaxis = list(title = "MLL")))

    print(p_ml)

    #scatter <- plotly::plot_ly(plot_ml, x = ~log10(lambda), y = ~log10(tau), z = ~mll, type = 'scatter3d', mode = 'markers',
    #                           marker = list(size = 4, color = ~mll, colorscale = 'Viridis')) %>%
    #  plotly::layout(title = paste0("Scatter plot Marginal Log Likelihood-Iteration:", t),
    #                 scene = list(xaxis = list(title = "log(λ)"), yaxis = list(title = "log(τ)"), zaxis = list(title = "MLL")))

    #print(scatter)

    t <- t + 1
  }

  doParallel::stopImplicitCluster()

  return(list(plot_ml, interp_df))
}
