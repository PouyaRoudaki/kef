
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
#' @param method_of_p_calculation A string that shows by which approach we calculate the p. If it is ''ordinary'' it means that we used the ordinary approach of density at sampled points and if it is ''grid_based'' approach it means that we approximate sample densities or probabilities using grid densities.
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
                          method_of_p_calculation = 1){

  # Compute the probabilities at the sampled points
  dens_sampled_x <- density_at_sampled_x(centered_kernel_mat_at_sampled, lambda_hat, weight_hat_vec)

  # Compute the densities at the grid points
  dens_grid <- density_at_grid(centered_kernel_mat_at_grid, centered_kernel_self_grid, lambda_hat, weight_hat_vec)

  # Normalize the density by the integral over the grid
  normalizing_cte <- trapz(x_grid, dens_grid)  # trapz is from pracma package

  # Prepare the output as a list of normalized probabilities
  dens_list <- list()
  dens_list$sampled_x <- dens_sampled_x / normalizing_cte
  dens_list$grid_x <- dens_grid / normalizing_cte

  prob_list <- list()
  prob_list$sampled_x <- den_list$sampled_x / sum(den_list$sampled_x)
  prob_list$grid_x <- den_list$grid_x / sum(den_list$grid_x)


  if(method_of_p_calculation == "grid_based"){

    approx_dens_or_prob <- get_grid_approx_dens_or_probs(sampled_x, x_grid, dens_list)

    dens_list$sampled_x <- approx_dens_or_prob$dens

    prob_list$sampled_x <- approx_dens_or_prob$prob

  }

  result_list <- list()
  if(type_of_p_is_prob == TRUE & type_of_q_is_prob = TRUE){

    result_list$sampled_x <- prob_list$sampled_x
    result_list$grid_x <- prob_list$grid_x

  } else if(type_of_p_is_prob == TRUE & type_of_q_is_prob = FALSE){

    result_list$sampled_x <- prob_list$sampled_x
    result_list$grid_x <- dens_list$grid_x

  } else if(type_of_p_is_prob == FALSE & type_of_q_is_prob = TRUE){

    result_list$sampled_x <- dens_list$sampled_x
    result_list$grid_x <- prob_list$grid_x

  } else if(type_of_p_is_prob == FALSE & type_of_q_is_prob = FALSE){

    result_list$sampled_x <- dens_list$sampled_x
    result_list$grid_x <- dens_list$grid_x

  }


  return(result_list)
}
