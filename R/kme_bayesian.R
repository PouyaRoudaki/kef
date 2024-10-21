#' Regularized Estimator of Kernel Mean Embedding(RKMSE)
#'
#' Compute the regularized estimator of Kernel Mean Embedding. Same shrinkage in all directions
#'
#' @param evaluate_at The point at which to evaluate the kernel mean embedding estimator function.
#' @param list_fixed A list of fixed points for the kernel mean embedding estimator function.
#' @param kernel_type Type of kernel function (default is "rbf" for Radial Basis Function).
#' @param kernel_params A list of kernel parameters, including length_scale, degree, free_add, free_mult, and nu_matern.
#'
#' @return The regularized estimator of Kernel Mean Embedding.
#' @export
#'
#' @author Pouya Roudaki \email{roudaki.pouya@@gmail.com}.
#'
#' @references Bayesian Learning of Kernel Embeddings,
#' Association for Uncertainty in Artificial Intelligence
#' Equations (7), (8) and Appendix (A.3)
#'
#' @examples
#' bayesian_est_KME(evaluate_at = c(5, 3, 8), list_fixed = list(c(0, 1, 2), c(1, 2, 3)),
#' kernel_type = "rbf",
#' kernel_params = list(length_scale = 1, degree = 2, free_add = 0, free_mult = 1, nu_matern = 1),
#' dom_measure_type = "gaussian",
#' dom_measure_parmas = list(eta = 1),
#' tau = 1)




bayesian_est_KME <- function(evaluate_at,
                             list_fixed,
                             kernel_type = "rbf",
                             kernel_params = list(length_scale = 1, degree = 2,
                                                  free_add = 0, free_mult = 1, nu_matern = 1),
                             dom_measure_type = "gaussian",
                             dom_measure_parmas = list(eta = 1),
                             tau = 1,
                             precomputed = list(R = NULL, inverse_regularizer = NULL)) {

  # Calculate the required material for calculation of mean and variance of Gaussian Process posterior Equation (8)


  if (is.null(precomputed$inverse_regularizer)){

    # Find convolutional gram matrix

    if (is.null(precomputed$R)){
      R <- conv_gram_matrix(vec_list = list_fixed,
                            kernel_type = kernel_type,
                            kernel_params = kernel_params,
                            dom_measure_type = dom_measure_type,
                            dom_measure_parmas = dom_measure_parmas)
    } else {
      R <- precomputed$R
    }

    ## Find number of fixed points
    n <- length(list_fixed)

    ## inverse matrix of regularizer
    inverse_regularizer <- solve(R + tau^2/n * diag(n))
  } else {

    ## inverse matrix of regularizer
    inverse_regularizer <- precomputed$inverse_regularizer
  }


  r_double_star <- r_conv_kernel(x = evaluate_at,
                                 y = evaluate_at,
                                 kernel_type = kernel_type,
                                 length_scale = kernel_params$length_scale,
                                 degree = kernel_params$degree,
                                 free_add = kernel_params$free_add,
                                 free_mult = kernel_params$free_mult,
                                 nu_matern = kernel_params$nu_matern,
                                 dom_measure_type = dom_measure_type,
                                 dom_measure_eta = dom_measure_parmas$eta)

  R_star <-  sapply(list_fixed,
                    r_conv_kernel,
                    x = evaluate_at,
                    kernel_type = kernel_type,
                    length_scale = kernel_params$length_scale,
                    degree = kernel_params$degree,
                    free_add = kernel_params$free_add,
                    free_mult = kernel_params$free_mult,
                    nu_matern = kernel_params$nu_matern,
                    dom_measure_type = dom_measure_type,
                    dom_measure_eta = dom_measure_parmas$eta)

  # For optimizing the code auxiliary calculation
  aux_matrix <- t(R_star) %*% inverse_regularizer

  # Calculate the standard estimator of Kernel Mean Embedding
  mu_hat_vec <- sapply(list_fixed,
                       std_est_KME,
                       list_fixed,
                       kernel_type = kernel_type,
                       kernel_params = kernel_params)



  # Calculate the mean of GP posterior at evaluation point Equation (7)
  result_mu_gp <-  aux_matrix %*% mu_hat_vec

  # Calculate the variance of GP posterior at evaluation point Equation (7)
  result_var_gp <-  r_double_star - aux_matrix %*% R_star

  result_bayesian <- c(post_mean = result_mu_gp, post_var = result_var_gp)

  return(result_bayesian)
}
