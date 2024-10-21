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
#' @references Kernel Mean Shrinkage Estimators, Muandet et al.,
#' Journal of Machine Learning Research 17 (2016) 1-41 Submitted 5/14; Revised 9/15; Published 4/16
#' Equations (22) and (25).
#'
#' @examples
#' reg_est_KME(c(5, 3, 8), list(c(0, 1, 2), c(1, 2, 3)), kernel_type = "rbf", kernel_params = list(length_scale = 1, degree = 2, free_add = 0, free_mult = 1, nu_matern = 1))
reg_est_KME <- function(evaluate_at,
                        list_fixed,
                        kernel_type = "rbf",
                        kernel_params = list(length_scale = 1, degree = 2,
                                             free_add = 0, free_mult = 1, nu_matern = 1),
                        precomputed = list(gram = NULL, lambda = NULL)) {

  # In order to find the shrinkage parameter lambda we need to find "rho" and
  # "rho with stroke" see Equation (25) in the reference.

  if (is.null(precomputed$lambda)) {

    if (is.null(precomputed$gram)) {
      gram <- gram_matrix(vec_list = list_fixed,
                          kernel_type = kernel_type,
                          kernel_params = kernel_params)
    } else {
      gram <- precomputed$gram
    }

    n = length(list_fixed)

    rho = 1/(n^2) * sum(gram)

    rho_with_stroke = 1/n * sum(diag(gram))

    lambda = n*(rho_with_stroke-rho) / ((n-1)*(n*rho-rho_with_stroke))
  } else {
    lambda <- precomputed$lambda
  }

  # Calculate the standard estimator of Kernel Mean Embedding
  result_std_KME <- std_est_KME(evaluate_at,
                                list_fixed,
                                kernel_type = "rbf",
                                kernel_params = kernel_params)


  # Calculate the regularized estimate
  result_reg_est <- 1/(1+lambda) * result_std_KME

  return(result_reg_est)
}
