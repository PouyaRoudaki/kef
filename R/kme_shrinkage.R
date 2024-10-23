#' Shrinkage Estimator of Kernel Mean Embedding(SKMSE)
#'
#' Compute the shrinkage estimator of Kernel Mean Embedding. More shrinkage in the direction of less variation
#' (More shrinkage in the direction of eigen vectors with lower eigen values)
#'
#' @param evaluate_at The point at which to evaluate the kernel mean embedding estimator function.
#' @param list_fixed A list of fixed points for the kernel mean embedding estimator function.
#' @param kernel_type Type of kernel function (default is "rbf" for Radial Basis Function).
#' @param kernel_params A list of kernel parameters, including length_scale, degree, free_add, free_mult, and nu_matern.
#'
#' @return The shrinkage estimator of Kernel Mean Embedding.
#' @export
#'
#' @author Pouya Roudaki \email{roudaki.pouya@@gmail.com}.
#'
#' @references Kernel Mean Shrinkage Estimators, Muandet et al.,
#' Journal of Machine Learning Research 17 (2016) 1-41 Submitted 5/14; Revised 9/15; Published 4/16,
#' Equation (32)
#'
#'
#' @examples
#' shr_est_KME(c(5, 3, 8), list(c(0, 1, 2), c(1, 2, 3)), kernel_type = "rbf", kernel_params = list(length_scale = 1, degree = 2, free_add = 0, free_mult = 1, nu_matern = 1))
shr_est_KME <- function(evaluate_at,
                        list_fixed,
                        kernel_type = "rbf",
                        kernel_params = list(length_scale = 1, degree = 2,
                                             free_add = 0, free_mult = 1, nu_matern = 1, centering_param = 7),
                        lambda_tunner = 1,
                        precomputed = list(gram = NULL, lambda = NULL, inverse_regularizer = NULL, beta_s = NULL)) {

  #Journal of Machine Learning Research 17 (2016) 1-41 Submitted 5/14; Revised 9/15; Published 4/16,
  # Equation (32)


  # Find beta_s in "1/sqrt(n) * t(kernel_vec)  %*% beta_s"
  if (is.null(precomputed$beta_s)) {

    # Find number of fixed points
    n = length(list_fixed)

    # Find solve(gram + n*lambda*diag(n))
    if (is.null(precomputed$inverse_regularizer)) {

      # Find the gram matrix
      if (is.null(precomputed$gram)) {
        gram <- gram_matrix(vec_list = list_fixed,
                            kernel_type = kernel_type,
                            kernel_params = kernel_params,
                            centering_param = centering_param)
      } else {
        gram <- precomputed$gram
      }

      # Find lambda
      if (is.null(precomputed$lambda)) {

        # Set lambda to be equal to largest non-negative eigen-value of gram matrix
        gram_eigens <- eigen(gram,symmetric = T,only.values = T)$values

        # Set lambda to be equal to smallest non-negative eigen-value of gram matrix
        lambda <- lambda_tunner * min(gram_eigens[gram_eigens>=0])
      } else {

        # LOOCV function of Shrinkage estimator and find the optimal lambda
        lambda <- precomputed$lambda
      }

      inverse_regularizer <- solve(gram + n*lambda*diag(n))

    } else {
      inverse_regularizer <- precomputed$inverse_regularizer
    }

    # 1_n vector
    one_n <- rep(1,n)/n

    # Find beta
    beta_s <- sqrt(n) * inverse_regularizer %*% gram %*%  one_n

  } else {
    beta_s <- precomputed$beta_s
  }


  # Kernel evaluation: Output is a vector. k(x,x_i)
  kernel_vec <- sapply(list_fixed,
                       kernel,
                       x = evaluate_at,
                       type = kernel_type,
                       length_scale = kernel_params$length_scale,
                       degree = kernel_params$degree,
                       free_add = kernel_params$free_add,
                       free_mult = kernel_params$free_mult,
                       nu_matern = kernel_params$nu_matern)

  # Calculate the shrinkage estimate
  result_shr_est <- 1/sqrt(n) * t(kernel_vec)  %*% beta_s

  return(result_shr_est)
}
