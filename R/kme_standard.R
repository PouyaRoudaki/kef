#' Standard Estimator of Kernel Mean Embedding(KMSE)
#'
#' Compute the standard estimator of Kernel Mean Embedding.
#'
#' @param evaluate_at The point at which to evaluate the kernel mean embedding estimator function.
#' @param list_fixed A list of fixed points for the kernel mean embedding estimator function.
#' @param kernel_type Type of kernel function (default is "rbf" for Radial Basis Function).
#' @param kernel_params A list of kernel parameters, including length_scale, degree, free_add, free_mult, and nu_matern.
#'
#' @return The standardized estimate using Kernel Mean Embedding.
#' @export
#'
#' @author Pouya Roudaki \email{roudaki.pouya@@gmail.com}.
#'
#' @examples
#' std_est_KME(c(5, 3, 8), list(c(0, 1, 2), c(1, 2, 3)), kernel_type = "rbf", kernel_params = list(length_scale = 1, degree = 2, free_add = 0, free_mult = 1, nu_matern = 1))
std_est_KME <- function(evaluate_at,
                        list_fixed,
                        kernel_type = "rbf",
                        kernel_params = list(length_scale = 1, degree = 2,
                                             free_add = 0, free_mult = 1, nu_matern = 1)) {

  # Calculate the sum of kernel evaluations for each fixed point
  kernel_sums <- sum(sapply(list_fixed,
                            kernel,
                            x = evaluate_at,
                            type = kernel_type,
                            length_scale = kernel_params$length_scale,
                            degree = kernel_params$degree,
                            free_add = kernel_params$free_add,
                            free_mult = kernel_params$free_mult,
                            nu_matern = kernel_params$nu_matern))

  # Calculate the standardized estimate
  result_std_est <- 1/length(list_fixed) * kernel_sums

  return(result_std_est)
}
