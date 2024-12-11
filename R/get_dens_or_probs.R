#' Density at Sampled Points
#'
#' This function calculates the density values at sampled data points using
#' a given centered kernel matrix, a scaling parameter (lambda_hat), and a vector of weights (weight_hat_vec).
#'
#' @param centered_kernel_mat_at_sampled A matrix where each column represents the
#' centered kernel values at different sampled data points.
#' @param lambda_hat A scalar parameter used to scale the contribution of the kernel
#' values and weights in the exponential function.
#' @param weight_hat_vec A vector of weights that are applied to the centered kernel
#' matrix in the computation.
#'
#' @return A vector containing the computed density values at each sampled data point.
#' @export
#'
#' @examples
#' # Example usage of the function
#' density_at_sampled_x(centered_kernel_mat_at_sampled, lambda_hat, weight_hat_vec)
density_at_sampled_x <- function(centered_kernel_mat_at_sampled, lambda_hat, weight_hat_vec) {

  # Extract the diagonal elements of the kernel matrix, which represent the
  # self-kernel values (i.e., the kernel values of each point with itself).
  diag_vals <- diag(centered_kernel_mat_at_sampled)

  # Compute the density values at each sampled data point
  p_vec <- sapply(1:nrow(centered_kernel_mat_at_sampled), function(i) {
    # For each sampled point, calculate the exponential of the linear combination of weights
    # and centered kernel values, adjusted by lambda_hat, and subtract 0.5 times the
    # self-kernel value for that point.
    exp(lambda_hat * (weight_hat_vec %*% centered_kernel_mat_at_sampled[,i] - 0.5 * diag_vals[i]))
  })

  # Return the vector of density values
  return(p_vec)
}



#' Density at Grid Points
#'
#' This function calculates the density values at specific grid points using
#' a given centered kernel matrix, self-grid kernel values, a scaling parameter (lambda_hat),
#' and a vector of weights (weight_hat_vec).
#'
#' @param centered_kernel_mat_at_grid A matrix where each column represents the
#' centered kernel values at different grid points.
#' @param centerd_kernel_self_grid A vector containing the centered kernel values
#' for the grid points when compared to themselves.
#' @param lambda_hat A scalar parameter used to scale the contribution of the kernel
#' values and weights in the exponential function.
#' @param weight_hat_vec A vector of weights that are applied to the centered kernel
#' matrix in the computation.
#'
#' @return A vector containing the computed density values at each grid point.
#' @export
#'
#' @examples
#' # Example usage of the function
#' density_at_grid(centered_kernel_mat_at_grid, centerd_kernel_self_grid, lambda_hat, weight_hat_vec)
density_at_grid <- function(centered_kernel_mat_at_grid, centerd_kernel_self_grid, lambda_hat, weight_hat_vec) {

  # Compute the density values at each grid point
  p_vec <- sapply(1:ncol(centered_kernel_mat_at_grid), function(i) {
    # For each grid point, calculate the exponential of the linear combination of weights
    # and centered kernel values, adjusted by lambda_hat, and subtract 0.5 times the
    # self-grid kernel value.
    exp(lambda_hat * (weight_hat_vec %*% centered_kernel_mat_at_grid[,i] - 0.5 * centerd_kernel_self_grid[i]))
  })

  # Return the vector of density values
  return(p_vec)
}



#' Compute Densities for Sampled and Grid Points
#'
#' This function computes densities for both sampled points and grid points.
#' The normalization is done by dividing the densities by the integral over the grid points.
#'
#' @param centered_kernel_mat_at_sampled A square matrix (n x n) where n is the number of sampled points.
#'        The matrix should be a centered kernel matrix evaluated at the sampled points.
#' @param centered_kernel_mat_at_grid A matrix (n x m) where n is the number of sampled points
#'        and m is the number of grid points. The matrix should be a centered kernel matrix
#'        evaluated at the grid points.
#' @param centered_kernel_self_grid A vector (length m) representing the diagonal of the centered kernel
#'        matrix evaluated at the grid points.
#' @param sampled_x A vector of sampled points where the probabilities are evaluated.
#' @param x_grid A vector of grid points where the probabilities are evaluated.
#' @param lambda_hat A scalar value representing the estimated lambda.
#' @param weight_hat_vec A vector of weights (length n) corresponding to the sampled points.
#' @param type_of_p_is_prob A Boolean that specifies if the type of p is probability (TRUE) or it is density (FALSE).
#' @param type_of_q_is_prob A Boolean that specifies if the type of q is probability (TRUE) or it is density (FALSE).
#' @param method_of_p_calculation A string that shows by which approach we calculate the p. If it is ''ordinary'' it means that we used the ordinary approach of density at sampled points and if it is ''neighborhood_grid'' approach it means that we approximate sample densities or probabilities using grid densities.
#'
#' @return A list containing two elements:
#'         \item{sampled_x}{A vector of densities or normalized probabilities at the sampled points.}
#'         \item{grid_x}{A vector of densities or normalized probabilities at the grid points.}
#' @export
#'
#' @examples
#' # Example usage (assuming inputs are defined):
#' probs <- get_dens_or_prob(centered_kernel_mat_at_sampled, centered_kernel_mat_at_grid,
#'                    centered_kernel_self_grid, sampled_x,
#'                    x_grid, lambda_hat, weight_hat_vec)
get_dens_or_prob <- function(centered_kernel_mat_at_sampled,
                             centered_kernel_mat_at_grid,
                          centered_kernel_self_grid,
                          sampled_x,
                          x_grid,
                          lambda_hat,
                          weight_hat_vec,
                          type_of_p_is_prob = TRUE,
                          type_of_q_is_prob = TRUE,
                          method_of_p_calculation = "ordinary"){

  # Compute the probabilities at the sampled points
  dens_sampled_x <- density_at_sampled_x(centered_kernel_mat_at_sampled, lambda_hat, weight_hat_vec)

  # Compute the densities at the grid points
  dens_grid <- density_at_grid(centered_kernel_mat_at_grid, centered_kernel_self_grid, lambda_hat, weight_hat_vec)

  # Normalize the density by the integral over the grid
  normalizing_cte <- pracma::trapz( x_grid, dens_grid)  # trapz is from pracma package

  # Prepare the output as a list of normalized probabilities
  dens_list <- list()
  dens_list$sampled_x <- dens_sampled_x / normalizing_cte
  dens_list$grid_x <- dens_grid / normalizing_cte

  prob_list <- list()
  prob_list$sampled_x <- dens_list$sampled_x / sum(dens_list$sampled_x)
  prob_list$grid_x <- dens_list$grid_x / sum(dens_list$grid_x)

  #print(paste("inside get_dens_or_probs",length(prob_list$grid_x)))
  if(method_of_p_calculation == "neighborhood_grid"){

    approx_dens_or_prob <- get_grid_approx_dens_or_probs_vectorized(sampled_x, x_grid, dens_list)

    dens_list$sampled_x <- approx_dens_or_prob$dens

    prob_list$sampled_x <- approx_dens_or_prob$prob

  }

  result_list <- list()

  if(type_of_p_is_prob){
    result_list$sampled_x <- prob_list$sampled_x
  }else{
    result_list$sampled_x <- dens_list$sampled_x
  }

  if(type_of_q_is_prob){
    result_list$grid_x <- prob_list$grid_x
  }else{
    result_list$grid_x <- dens_list$grid_x
  }

  #print(paste("inside get_dens_or_probs",length(result_list$grid_x)))

  return(result_list)
}

#' Compute Densities for Sampled points without a grid
#'
#' This function computes densities for both sampled points and grid points.
#' The normalization is done by dividing the densities by the integral over the grid points.
#'
#' @param centered_kernel_mat_at_sampled A square matrix (n x n) where n is the number of sampled points.
#'        The matrix should be a centered kernel matrix evaluated at the sampled points.
#' @param sampled_x A vector of sampled points where the probabilities are evaluated.
#' @param min_x A scalar value representing the minimum of x domain.
#' @param max_x A scalar value representing the maximum of x domain.
#' @param lambda_hat A scalar value representing the estimated lambda.
#' @param weight_hat_vec A vector of weights (length n) corresponding to the sampled points.
#' @param type_of_p_is_prob A Boolean that specifies if the type of p is probability (TRUE) or it is density (FALSE).
#'
#' @return A list containing an element:
#'         \item{sampled_x}{A vector of densities or normalized probabilities at the sampled points.}
#' @export
#'
#' @examples
#' # Example usage (assuming inputs are defined):
#' probs <- get_dens_or_prob_wo_grid(centered_kernel_mat_at_sampled,
#'                                   sampled_x,
#'                    x_grid, lambda_hat, weight_hat_vec)
get_dens_or_prob_wo_grid <- function(centered_kernel_mat_at_sampled,
                             min_x,
                             max_x,
                             sampled_x,
                             lambda_hat,
                             weight_hat_vec,
                             type_of_p_is_prob = TRUE){

  # Compute the probabilities at the sampled points
  dens_sampled_x <- density_at_sampled_x(centered_kernel_mat_at_sampled, lambda_hat, weight_hat_vec)

  # Compute the densities at the grid points
  #dens_grid <- density_at_grid(centered_kernel_mat_at_grid, centered_kernel_self_grid, lambda_hat, weight_hat_vec)
  # Find the base measure of samples
  sample_mid_points <- get_middle_points_grid(x_grid[1], sampled_x, x_grid[length(x_grid)])
  base_measure_weights <- sample_mid_points[-1] - sample_mid_points[-length(sample_mid_points)]

  dens_base_sampled_x <- dens_sampled_x * base_measure_weights

  # Normalize the density by the integral over the grid
  normalizing_cte <- pracma::trapz( sampled_x, dens_base_sampled_x)  # trapz is from pracma package

  # Prepare the output as a list of normalized probabilities
  dens_list <- list()
  dens_list$sampled_x <- dens_sampled_x / normalizing_cte
  #dens_list$grid_x <- dens_grid / normalizing_cte

  prob_list <- list()
  prob_list$sampled_x <- dens_list$sampled_x / sum(dens_list$sampled_x)
  prob_list$grid_x <- dens_list$grid_x / sum(dens_list$grid_x)

  result_list <- list()

  if(type_of_p_is_prob){
    result_list$sampled_x <- prob_list$sampled_x
  }else{
    result_list$sampled_x <- dens_list$sampled_x
  }

  #print(paste("inside get_dens_or_probs",length(result_list$grid_x)))

  return(result_list)
}
