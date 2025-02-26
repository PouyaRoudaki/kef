#' Estimate Weights Using the Newton-Raphson Method
#'
#' This function estimates the weight vector using an iterative Newton-Raphson method.
#' The method updates weights based on the provided kernel matrices, regularization
#' parameters, sampled points, and grid points.
#'
#' @param lambda_hat A scalar representing the estimated lambda parameter, which
#'        controls the contribution of the kernel matrices to the weight estimation.
#' @param tau_hat A scalar representing the estimated tau parameter, which
#'        determines the regularization strength applied to the weights.
#' @param centered_kernel_mat_at_sampled A square matrix (n x n) representing the centered
#'        kernel matrix evaluated at the sampled points, where n is the number of sampled points.
#' @param centered_kernel_mat_at_grid A matrix (n x m) representing the centered kernel
#'        matrix evaluated at the grid points, where n is the number of sampled points
#'        and m is the number of grid points.
#' @param centered_kernel_self_grid A vector of length m representing the diagonal of
#'        the centered kernel matrix evaluated at the grid points.
#' @param sampled_x A vector of sampled points for which the weights are to be estimated.
#' @param x_grid A vector of grid points that define the evaluation domain for the kernel.
#' @param type_of_p_is_prob Logical; if TRUE, treats the "p" function as a probability density.
#' @param type_of_q_is_prob Logical; if TRUE, treats the "q" function as a probability density.
#' @param method_of_p_calculation A string indicating the method used to calculate probabilities
#'        for the "p" function. Default is "ordinary".
#' @param print_trace Logical; if TRUE, prints progress updates during the Newton-Raphson iterations.
#'
#' @return A numeric vector of length n representing the estimated weight vector for the sampled points.
#'
#' @export
get_weights <- function(lambda_hat,
                        tau_hat,
                        centered_kernel_mat_at_sampled,
                        centered_kernel_mat_at_grid,
                        centered_kernel_self_grid,
                        sampled_x,
                        x_grid,
                        type_of_p_is_prob=TRUE,
                        type_of_q_is_prob=TRUE,
                        method_of_p_calculation="ordinary",
                        print_trace = FALSE) {

  max_iteration <- 2000  # Maximum number of iterations for the Newton-Raphson method
  NRstepsize <- 0.1  # Step size for the Newton-Raphson update
  n <- nrow(centered_kernel_mat_at_sampled)  # Number of sampled points
  weight_hat_vec <- rep(0, n)  # Initialize the weight vector with zeros
  #s <- rep(1000, n)
  #weight_hat_change <- rep(1000, n)
  #counter <- 1
  #while ((norm(s, p = 2) > 10^(-10)) & (norm(weight_hat_change, p = 2) > 10^(-10))) {


  for (i in 1:max_iteration) {
    # Calculate probabilities for sampled and grid points
    probs <- get_dens_or_prob(centered_kernel_mat_at_sampled,
                              centered_kernel_mat_at_grid,
                       centered_kernel_self_grid,
                       sampled_x,
                       x_grid,
                       lambda_hat,
                       weight_hat_vec,
                       type_of_p_is_prob,
                       type_of_q_is_prob,
                       method_of_p_calculation)

    prob_sampled_x <- probs$sampled_x
    prob_grid_x <- probs$grid_x

    #print(probs)
    #print(length(probs$grid_x))
    #print(length(prob_grid_x))
    #print(dim(centered_kernel_mat_at_grid))
    # Compute the gradient (s)
    s <- lambda_hat * (colSums(centered_kernel_mat_at_sampled) -
                         n * prob_grid_x %*% t(centered_kernel_mat_at_grid)) -
      tau_hat * weight_hat_vec / prob_sampled_x

    # Compute the inverse of the Hessian matrix
    if(type_of_q_is_prob == TRUE){
      Hessian <- lambda_hat^2 * n * (centered_kernel_mat_at_grid %*%
                             diag(prob_grid_x) %*% t(centered_kernel_mat_at_grid) -
                             (centered_kernel_mat_at_grid %*% prob_grid_x) %*%
                             t(centered_kernel_mat_at_grid %*% prob_grid_x)  ) +
                             diag(tau_hat / prob_sampled_x)

      Hessian_inv <- solve(Hessian)
    }else{
      Hessian_inv <- solve(lambda_hat^2 * n * centered_kernel_mat_at_grid %*%
                             diag(prob_grid_x) %*% t(centered_kernel_mat_at_grid) +
                             diag(tau_hat / prob_sampled_x))
    }

    #print(s)
    #print(Hessian_inv)
    # Update the weight vector using the Newton-Raphson method
    #weight_hat_change <- NRstepsize * s %*% Hessian_inv
    weight_hat_vec <- weight_hat_vec + NRstepsize * s %*% Hessian_inv

    # Print progress every 10% of the iterations or at the first iteration
    if ((i %% round(max_iteration / 10) == 0 || i == 1) & print_trace == TRUE) {
      print(paste("Iteration", i, ": ||s||_2 =", pracma::Norm(s)))
      #print(summary(as.vector(Hessian)))
    }
    #counter = counter + 1
  }

  return(weight_hat_vec)
}

#' Estimate Weights Using the Newton-Raphson Method
#'
#' This function estimates the weight vector using an iterative Newton-Raphson method.
#' The method updates weights based on the provided kernel matrices, regularization
#' parameters, sampled points, and grid points.
#'
#' @param lambda_hat A scalar representing the estimated lambda parameter, which
#'        controls the contribution of the kernel matrices to the weight estimation.
#' @param tau_hat A scalar representing the estimated tau parameter, which
#'        determines the regularization strength applied to the weights.
#' @param centered_kernel_mat_at_sampled A square matrix (n x n) representing the centered
#'        kernel matrix evaluated at the sampled points, where n is the number of sampled points.
#' @param sampled_x A vector of sampled points for which the weights are to be estimated.
#' @param min_x A scalar representing the minimum value of the domain.
#' @param max_x A scalar representing the maximum value of the domain.
#' @param print_trace Logical; if TRUE, prints progress updates during the Newton-Raphson iterations.
#'
#' @return A numeric vector of length n representing the estimated weight vector for the sampled points.
#'
#' @export
get_weights_wo_grid <- function(lambda_hat,
                        tau_hat,
                        centered_kernel_mat_at_sampled,
                        sampled_x,
                        min_x,
                        max_x,
                        print_trace = FALSE) {

  max_iteration <- 1000  # Maximum number of iterations for the Newton-Raphson method
  NRstepsize <- 0.5  # Step size for the Newton-Raphson update
  n <- nrow(centered_kernel_mat_at_sampled)  # Number of sampled points
  weight_hat_vec <- rep(0,n)  # Initialize the weight vector with zeros

  #print(summary(base_measure_weights))
  for (i in 1:max_iteration) {
    # Calculate probabilities for sampled and grid points
    dens <- get_dens_wo_grid(centered_kernel_mat_at_sampled,
                              min_x,
                              max_x,
                              sampled_x,
                              lambda_hat,
                              weight_hat_vec)

    # Find the base measure of samples
    sample_mid_points <- get_middle_points_grid(min_x, sampled_x, max_x)
    base_measure_weights <- sample_mid_points[-1] - sample_mid_points[-length(sample_mid_points)]

    dens_sampled_base <- dens * base_measure_weights

    prob_sampled_base <- dens_sampled_base / sum(dens_sampled_base)
    prob_sampled <- dens / sum(dens)

    s <- lambda_hat * (colSums(centered_kernel_mat_at_sampled) -
                         n * (prob_sampled_base) %*% t(centered_kernel_mat_at_sampled)) -
      tau_hat * weight_hat_vec / prob_sampled

    # Compute the inverse of the Hessian matrix
    Hessian <- lambda_hat^2 * n * (centered_kernel_mat_at_sampled %*%
                                                 diag(prob_sampled_base) %*%
                                                 t(centered_kernel_mat_at_sampled) -
                                                 (centered_kernel_mat_at_sampled %*% prob_sampled_base) %*%
                                                 t(centered_kernel_mat_at_sampled %*% prob_sampled_base)) +
                                                 diag(tau_hat / prob_sampled)



    Hessian_inv <- solve(Hessian)
    #print(summary(as.vector(Hessian_inv)))

    #print(s)
    #print(Hessian_inv)
    # Update the weight vector using the Newton-Raphson method
    #weight_hat_change <- NRstepsize * s %*% Hessian_inv
    weight_hat_vec <- weight_hat_vec + NRstepsize * s %*% Hessian_inv

    # Print progress every 10% of the iterations or at the first iteration
    if ((i %% round(max_iteration / 10) == 0 || i == 1) & print_trace == TRUE) {

      print(paste("Iteration", i, ": ||s||_2 =", format(pracma::Norm(s), digits = 3, scientific = T) ))
      #print(summary(as.vector(Hessian)))
      #print(summary(as.vector(solve(Hessian))))
    }
    #counter = counter + 1
  }

  return(weight_hat_vec)
}


#' Estimate Weights Using the Newton-Raphson Method for Monte Carlo Approach.
#'
#' This function estimates the weight vector using an iterative Newton-Raphson method.
#' The method updates weights based on the provided kernel matrices, regularization
#' parameters, sampled points, and grid points.
#'
#' @param lambda_t A scalar representing the estimated lambda parameter, which
#'        controls the contribution of the kernel matrices to the weight estimation.
#' @param tau_t A scalar representing the estimated tau parameter, which
#'        determines the regularization strength applied to the weights.
#' @param centered_kernel_mat_at_sampled A square matrix (n x n) representing the centered
#'        kernel matrix evaluated at the sampled points, where n is the number of sampled points.
#' @param sampled_x A vector of sampled points for which the weights are to be estimated.
#' @param min_x A scalar representing the minimum value of the domain.
#' @param max_x A scalar representing the maximum value of the domain.
#' @param p_vec_t_1 A vector representing the probabilities of step t-1
#' @param print_trace Logical; if TRUE, prints progress updates during the Newton-Raphson iterations.
#'
#' @return A numeric vector of length n representing the estimated weight vector for the sampled points.
#'
#' @export
get_weights_wo_grid_mll <- function(lambda_t,
                                tau_t,
                                centered_kernel_mat_at_sampled,
                                sampled_x,
                                min_x,
                                max_x,
                                p_vec_t_1,
                                print_trace = FALSE) {

  max_iteration <- 1000  # Maximum number of iterations for the Newton-Raphson method
  NRstepsize <- 0.5  # Step size for the Newton-Raphson update
  n <- nrow(centered_kernel_mat_at_sampled)  # Number of sampled points
  weight_t_vec <- rep(0,n)  # Initialize the weight vector with zeros

  #print(summary(base_measure_weights))
  for (i in 1:max_iteration) {
    # Calculate probabilities for sampled and grid points
    dens <- as.numeric(get_dens_wo_grid(centered_kernel_mat_at_sampled,
                             min_x,
                             max_x,
                             sampled_x,
                             lambda_t,
                             weight_t_vec))

    # Find the base measure of samples
    sample_mid_points <- as.numeric(get_middle_points_grid(min_x, sampled_x, max_x))
    base_measure_weights <- sample_mid_points[-1] - sample_mid_points[-length(sample_mid_points)]

    dens_sampled_base <- dens * base_measure_weights

    prob_sampled_base <- dens_sampled_base / sum(dens_sampled_base)
    prob_sampled <- p_vec_t_1

    s <- lambda_t * (colSums(centered_kernel_mat_at_sampled) -
                         n * (prob_sampled_base) %*% t(centered_kernel_mat_at_sampled)) -
      tau_t * weight_t_vec / prob_sampled

    # Compute the inverse of the Hessian matrix
    Hessian <- lambda_t^2 * n * (centered_kernel_mat_at_sampled %*%
                                     diag(prob_sampled_base) %*%
                                     t(centered_kernel_mat_at_sampled) -
                                     (centered_kernel_mat_at_sampled %*% prob_sampled_base) %*%
                                     t(centered_kernel_mat_at_sampled %*% prob_sampled_base)) +
      diag(tau_t / prob_sampled)



    Hessian_inv <- solve(Hessian)
    #print(summary(as.vector(Hessian_inv)))

    #print(s)
    #print(Hessian_inv)
    # Update the weight vector using the Newton-Raphson method
    #weight_hat_change <- NRstepsize * s %*% Hessian_inv
    weight_t_vec <- weight_t_vec + NRstepsize * s %*% Hessian_inv

    # Print progress every 10% of the iterations or at the first iteration
    if ((i %% round(max_iteration / 10) == 0 || i == 1) & print_trace == TRUE) {

      print(paste("Iteration", i, ": ||s||_2 =", pracma::Norm(s)))
      #print(summary(as.vector(Hessian)))
      #print(summary(as.vector(solve(Hessian))))
    }
    #counter = counter + 1
  }

  return(weight_t_vec)
}


#' Estimate Weights Using the Fixed-point Method
#'
#' This function estimates the weight vector using an iterative fixed_point method.
#' The method updates weights based on the provided kernel matrices, regularization
#' parameters, sampled points, and grid points.
#'
#' @param lambda_hat A scalar representing the estimated lambda parameter, which
#'        controls the contribution of the kernel matrices to the weight estimation.
#' @param tau_hat A scalar representing the estimated tau parameter, which
#'        determines the regularization strength applied to the weights.
#' @param centered_kernel_mat_at_sampled A square matrix (n x n) representing the centered
#'        kernel matrix evaluated at the sampled points, where n is the number of sampled points.
#' @param centered_kernel_mat_at_grid A matrix (n x m) representing the centered kernel
#'        matrix evaluated at the grid points, where n is the number of sampled points
#'        and m is the number of grid points.
#' @param centered_kernel_self_grid A vector of length m representing the diagonal of
#'        the centered kernel matrix evaluated at the grid points.
#' @param sampled_x A vector of sampled points for which the weights are to be estimated.
#' @param x_grid A vector of grid points that define the evaluation domain for the kernel.
#' @param type_of_p_is_prob Logical; if TRUE, treats the "p" function as a probability density.
#' @param type_of_q_is_prob Logical; if TRUE, treats the "q" function as a probability density.
#' @param method_of_p_calculation A string indicating the method used to calculate probabilities
#'        for the "p" function. Default is "ordinary".
#' @param print_trace Logical; if TRUE, prints progress updates during the Newton-Raphson iterations.
#'
#' @return A numeric vector of length n representing the estimated weight vector for the sampled points.
#'
#' @export

get_weights_fixed_point <- function(lambda_hat,
                                    tau_hat,
                                    centered_kernel_mat_at_sampled,
                                    centered_kernel_mat_at_grid,
                                    centered_kernel_self_grid,
                                    sampled_x,
                                    x_grid,
                                    type_of_p_is_prob=TRUE,
                                    type_of_q_is_prob=TRUE,
                                    method_of_p_calculation="ordinary",
                                    print_trace = FALSE){

  n <- nrow(centered_kernel_mat_at_sampled)  # Number of sampled points
  weight_hat_vec <- rep(1,n)  # Initialize the weight vector with zeros

  one_vec <- rep(1,n)# Initialize the weight vector with 1

  # Find the base measure of samples
  sample_mid_points <- get_middle_points_grid(x_grid[1], sampled_x, x_grid[length(x_grid)])
  base_measure_weights <- sample_mid_points[-1] - sample_mid_points[-length(sample_mid_points)]


  weight_hat_change <- rep(1000, n)
  counter <- 1

  while(pracma::Norm(weight_hat_change) > 10^(-14)){
    # Calculate probabilities for sampled and grid points
    probs <- get_dens_or_prob(centered_kernel_mat_at_sampled,
                              centered_kernel_mat_at_grid,
                              centered_kernel_self_grid,
                              sampled_x,
                              x_grid,
                              lambda_hat,
                              weight_hat_vec,
                              type_of_p_is_prob,
                              type_of_q_is_prob,
                              method_of_p_calculation)

    prob_sampled_x <- probs$sampled_x

    old_weight_hat_vec <- weight_hat_vec

    weight_hat_vec <- lambda_hat/tau_hat * prob_sampled_x *
      (centered_kernel_mat_at_sampled %*% (one_vec -
                                             n * prob_sampled_x * base_measure_weights))

    weight_hat_vec <- as.vector(weight_hat_vec)

    weight_hat_change <- weight_hat_vec - old_weight_hat_vec

    if ((counter %% 10 == 0 || counter == 1) & print_trace == TRUE) {
      #print(length(old_weight_hat_vec))
      #print(length(weight_hat_vec))
      print(paste("Iteration", counter, ": weights change norm =", pracma::Norm(weight_hat_change)))
    }

    counter = counter + 1
  }

  return(weight_hat_vec)
}

#' Estimate Weights Using the Gradient Descent Method
#'
#' This function estimates the weight vector using an iterative gradient descent method.
#' The method updates weights based on the provided kernel matrices, regularization
#' parameters, sampled points, and grid points.
#'
#' @param lambda_hat A scalar representing the estimated lambda parameter, which
#'        controls the contribution of the kernel matrices to the weight estimation.
#' @param tau_hat A scalar representing the estimated tau parameter, which
#'        determines the regularization strength applied to the weights.
#' @param centered_kernel_mat_at_sampled A square matrix (n x n) representing the centered
#'        kernel matrix evaluated at the sampled points, where n is the number of sampled points.
#' @param centered_kernel_mat_at_grid A matrix (n x m) representing the centered kernel
#'        matrix evaluated at the grid points, where n is the number of sampled points
#'        and m is the number of grid points.
#' @param centered_kernel_self_grid A vector of length m representing the diagonal of
#'        the centered kernel matrix evaluated at the grid points.
#' @param sampled_x A vector of sampled points for which the weights are to be estimated.
#' @param x_grid A vector of grid points that define the evaluation domain for the kernel.
#' @param type_of_p_is_prob Logical; if TRUE, treats the "p" function as a probability density.
#' @param type_of_q_is_prob Logical; if TRUE, treats the "q" function as a probability density.
#' @param method_of_p_calculation A string indicating the method used to calculate probabilities
#'        for the "p" function. Default is "ordinary".
#' @param print_trace Logical; if TRUE, prints progress updates during the Newton-Raphson iterations.
#'
#' @return A numeric vector of length n representing the estimated weight vector for the sampled points.
#'
#' @export
get_weights_gd <- function(lambda_hat,
                                    tau_hat,
                                    centered_kernel_mat_at_sampled,
                                    centered_kernel_mat_at_grid,
                                    centered_kernel_self_grid,
                                    sampled_x,
                                    x_grid,
                                    type_of_p_is_prob=TRUE,
                                    type_of_q_is_prob=TRUE,
                                    method_of_p_calculation="ordinary",
                                    print_trace = FALSE){

  n <- nrow(centered_kernel_mat_at_sampled)  # Number of sampled points
  weight_hat_vec <- rep(0,n)  # Initialize the weight vector with zeros

  one_vec <- rep(1,n)# Initialize the weight vector with 1

  # Find the base measure of samples
  sample_mid_points <- get_middle_points_grid(x_grid[1], sampled_x, x_grid[length(x_grid)])
  #base_measure_weights <- sample_mid_points[-1] - sample_mid_points[-length(sample_mid_points)]


  weight_hat_change <- rep(1000, n)
  counter <- 1

  GD_step_size <- 0.1^counter

  while(pracma::Norm(weight_hat_change) > 10^(-14)){
    # Calculate probabilities for sampled and grid points
    probs <- get_dens_or_prob(centered_kernel_mat_at_sampled,
                              centered_kernel_mat_at_grid,
                              centered_kernel_self_grid,
                              sampled_x,
                              x_grid,
                              lambda_hat,
                              weight_hat_vec,
                              type_of_p_is_prob,
                              type_of_q_is_prob,
                              method_of_p_calculation)

    prob_sampled_x <- probs$sampled_x

    print("weight")
    print(summary(weight_hat_vec))

    print("prob")
    print(summary(prob_sampled_x))
    old_weight_hat_vec <- weight_hat_vec

    s <- lambda_hat * (colSums(centered_kernel_mat_at_sampled) -
                         n * (base_measure_weights * prob_sampled_x) %*% t(centered_kernel_mat_at_sampled)) -
      tau_hat * weight_hat_vec / prob_sampled_x

    print(dim(s))

    weight_hat_vec <- weight_hat_vec + 0.00001^counter * s

    weight_hat_vec <- as.vector(weight_hat_vec)

    weight_hat_change <- weight_hat_vec - old_weight_hat_vec

    if ((counter %% 1 == 0 || counter == 1) & print_trace == TRUE) {
      print(summary(weight_hat_change))
      #print(length(weight_hat_vec))
      print(paste("Iteration", counter, ": weights change norm =", pracma::Norm(weight_hat_change)))
    }

    counter = counter + 1
  }

  return(weight_hat_vec)
}

#' Estimate Weights Using the Gradient Descent Method
#'
#' This function estimates the weight vector using an iterative gradient descent method.
#' The method updates weights based on the provided kernel matrices, regularization
#' parameters, sampled points, and grid points.
#'
#' @param lambda_hat A scalar representing the estimated lambda parameter, which
#'        controls the contribution of the kernel matrices to the weight estimation.
#' @param tau_hat A scalar representing the estimated tau parameter, which
#'        determines the regularization strength applied to the weights.
#' @param centered_kernel_mat_at_sampled A square matrix (n x n) representing the centered
#'        kernel matrix evaluated at the sampled points, where n is the number of sampled points.
#' @param centered_kernel_mat_at_grid A matrix (n x m) representing the centered kernel
#'        matrix evaluated at the grid points, where n is the number of sampled points
#'        and m is the number of grid points.
#' @param centered_kernel_self_grid A vector of length m representing the diagonal of
#'        the centered kernel matrix evaluated at the grid points.
#' @param sampled_x A vector of sampled points for which the weights are to be estimated.
#' @param x_grid A vector of grid points that define the evaluation domain for the kernel.
#' @param type_of_p_is_prob Logical; if TRUE, treats the "p" function as a probability density.
#' @param type_of_q_is_prob Logical; if TRUE, treats the "q" function as a probability density.
#' @param method_of_p_calculation A string indicating the method used to calculate probabilities
#'        for the "p" function. Default is "ordinary".
#' @param print_trace Logical; if TRUE, prints progress updates during the Newton-Raphson iterations.
#'
#' @return A numeric vector of length n representing the estimated weight vector for the sampled points.
#'
#' @export
get_weights_gd <- function(lambda_hat,
                           tau_hat,
                           centered_kernel_mat_at_sampled,
                           centered_kernel_mat_at_grid,
                           centered_kernel_self_grid,
                           sampled_x,
                           x_grid,
                           type_of_p_is_prob=TRUE,
                           type_of_q_is_prob=TRUE,
                           method_of_p_calculation="ordinary",
                           print_trace = FALSE){

  n <- nrow(centered_kernel_mat_at_sampled)  # Number of sampled points
  weight_hat_vec <- rep(0,n)  # Initialize the weight vector with zeros

  one_vec <- rep(1,n)# Initialize the weight vector with 1

  # Find the base measure of samples
  sample_mid_points <- get_middle_points_grid(x_grid[1], sampled_x, x_grid[length(x_grid)])
  #base_measure_weights <- sample_mid_points[-1] - sample_mid_points[-length(sample_mid_points)]


  weight_hat_change <- rep(1000, n)
  counter <- 1

  GD_step_size <- 0.1^counter

  while(pracma::Norm(weight_hat_change) > 10^(-14)){
    # Calculate probabilities for sampled and grid points
    probs <- get_dens_or_prob(centered_kernel_mat_at_sampled,
                              centered_kernel_mat_at_grid,
                              centered_kernel_self_grid,
                              sampled_x,
                              x_grid,
                              lambda_hat,
                              weight_hat_vec,
                              type_of_p_is_prob,
                              type_of_q_is_prob,
                              method_of_p_calculation)

    prob_sampled_x <- probs$sampled_x

    print("weight")
    print(summary(weight_hat_vec))

    print("prob")
    print(summary(prob_sampled_x))
    old_weight_hat_vec <- weight_hat_vec

    s <- lambda_hat * (colSums(centered_kernel_mat_at_sampled) -
                         n * (base_measure_weights * prob_sampled_x) %*% t(centered_kernel_mat_at_sampled)) -
      tau_hat * weight_hat_vec / prob_sampled_x

    print(dim(s))

    weight_hat_vec <- weight_hat_vec + 0.00001^counter * s

    weight_hat_vec <- as.vector(weight_hat_vec)

    weight_hat_change <- weight_hat_vec - old_weight_hat_vec

    if ((counter %% 1 == 0 || counter == 1) & print_trace == TRUE) {
      print(summary(weight_hat_change))
      #print(length(weight_hat_vec))
      print(paste("Iteration", counter, ": weights change norm =", pracma::Norm(weight_hat_change)))
    }

    counter = counter + 1
  }

  return(weight_hat_vec)
}


