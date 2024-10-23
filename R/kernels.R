#' @name kernels
#'
#' @param x            Kernel free argument, the location to evaluate KDE (single scalar or vector)
#' @param fixed        Kernel fixed argument (typically sample data vector or scalar)
#' @param type       Kernel name (\code{default = "rbf"})
#' @param length_scale       Kernel parameters, length scale
#' @param degree  Kernel parameters, polynomial degree
#' @param free_add  Kernel parameters, free additive parameter
#' @param free_mult   Kernel parameters, free multiplicative parameter
#' @param nu_matern  Kernel parameters, nu_matern for matern kernel
#'
#' , polynomial degree, free parameters, etc. for the kernel (as half-width of the kernel) or \code{NULL}
#'
#' @title Kernel functions
#'
#' @description Functions for commonly used kernels for kernel density estimation. The density and cumulative distribution functions are provided.
#'
#'
#'
#'
#'
#'
#' @seealso
#' \code{\link{norm}}, \code{\link{exp}}, \code{\link{sin}}, \code{\link{besselK}}
#'
#'
#'
#'
#' @author Pouya Roudaki \email{roudaki.pouya@@gmail.com}.
#'
#' @export
#'
#' @examples
#' kernel(c(1,2,3),c(2,1,3), type = "linear")
#'
#' @keywords kernel
#'
#'
#'

# Single Argument Approach
kernel <- function(x, fixed, type = "rbf", length_scale = 1, degree = 2,
                   free_add = 0, free_mult = 1,
                   nu_matern = 1, centering_param = 7) {
  type <- tolower(type)  # Convert type to lowercase for case-insensitivity

  switch(type,
         uniform = {
           # Implementation for the uniform (UNI) kernel
           # Required parameter is length_scale: the entire length of the uniform kernel
           # x vector dim
           d =  length(x)

           z = 2 * (norm(x - fixed, p = Inf)) / length_scale

           result_uni <- (abs(z) <= 1)

           return(result_uni)
         },
         linear = {
           # Implementation for the linear (LIN) kernel
           # Required parameter is free_param: constant that is added to the linear kernel
           # Inner product in R: as.numeric(x %*% y) as.numeric since we want a scalar, not an array

           result_lin <- free_mult * as.numeric(x %*% fixed) + free_add

           return(result_lin)
         },
         polynomial = {
           # Implementation for the polynomial (POL) kernel
           # Required parameters are the degree of the polynomial and the free_param
           # which is a constant that is added to the linear kernel before powering to the degree

           result_pol <- (free_mult * as.numeric(x %*% fixed) + free_add)^degree

           return(result_pol)
         },
         rq = {
           # Implementation for the rational quadratic (RQ) kernel
           # Required parameters are the degree of the polynomial and the length scale

           result_rq <- (1 + norm(x - fixed, p = 2)^2 / (2 * degree * (length_scale^2)))^(-degree)

           return(result_rq)
         },
         rbf = {
           # Implementation for the radial basis function (RBF) kernel
           # Required parameter is the length scale

           result_rbf <- free_mult * exp(-norm(x - fixed, p = 2)^2 / (2 * length_scale^2))

           return(result_rbf)
         },
         gibbs = {
           # Implementation for the Gibbs (GIB) kernel
           # Use different functions for the length scale in each direction
           #  k_{\text {Gibbs }}\left(\mathbf{x}, \mathbf{x}^{\prime}\right)=
           #  \prod_{p=1}^P\left(\frac{2 l_p(\mathbf{x}) l_p
           #  \left(\mathbf{x}^{\prime}\right)}{l_p^2(\mathbf{x})+
           #  l_p^2\left(\mathbf{x}^{\prime}\right)}\right)^{\frac{1}{2}}
           #  \exp \left(-\sum_{p=1}^P \frac{\left(x_p-x_p^{\prime}\right)^2}
           #  {l_p^2(\mathbf{x})+l_p^2\left(\mathbf{x}^{\prime}\right)}\right)
           #  where $x_p$ is the $p^{\text {th }}$ component of $\mathbf{x}$
           #
           #For example here we consider l_p = l*x_p
           #
           exp_term <- exp(-sum(mapply(aux_gibbs_exp_term, x, fixed, length_scale)))

           result_gibbs <- prod(mapply(aux_gibbs_kernel, x, fixed, length_scale, exp_term))

           return(result_gibbs)
         },
         periodic = {
           # Implementation for the periodic (PER) kernel
           # Required parameter is the length scale

           result_per <- free_mult * exp(-2 * (sin(norm(x - fixed, p = 1)))^2 / (length_scale^2))

           return(result_per)
         },
         ou = {
           # Implementation for the Ornstein-Uhlenbeck kernel (OU) kernel
           # Required parameter is the length scale
           result_ou <- exp(-norm(x - fixed, p = 2) / length_scale)

           return(result_ou)
         },
         matern = {
           # Implementation for the Matern kernel
           # Required parameters are the length scale and nu_matern
           # k_{\text {Matern }}\left(\mathbf{x}, \mathbf{x}^{\prime}\right)=
           # \frac{2^{1-\nu}}
           # {\Gamma(\nu)}\left(\frac{\sqrt{2 \nu}\left|\mathbf{x}-\mathbf{x}^{\prime}\right|}{l}\right)
           # K_\nu\left(\frac{\sqrt{2 \nu}\left|\mathbf{x}-\mathbf{x}^{\prime}\right|}{l}\right),
           #
           # where $K_\nu$ is a modified Bessel function.
           #

           aux <- sqrt(2 * nu_matern) * norm(x - fixed, p = 1) / length_scale

           result_matern <- (2^(1 - nu_matern)) / gamma(nu_matern) *
             (aux)^nu_matern * besselK(aux, nu = nu_matern, expon.scaled = FALSE)

           return(result_matern)
         },
         brownian = {
           # Implementation for the Brownian kernel
           # No parameter is required!
           result_brownian <- length_scale * 1/2 * (norm(x,1) + norm(fixed,1) - norm(x - fixed,1))
           #result_brownian <- length_scale^2 * min(x, fixed)

           return(result_brownian)
         },
         gaussian_centered_brownian = {
           # Implementation for the Brownian kernel
           # No parameter is required!
           result_centered_brownian <- length_scale * 1/2 *
             (-norm(x - fixed)
              -sqrt(2/pi) * (1- exp(-x^2/2)) + x*(1-2*pnorm(x))
              -sqrt(2/pi) * (1- exp(-fixed^2/2)) + fixed*(1-2*pnorm(fixed)))
           + length_scale * (-1 + sqrt(2))/sqrt(pi)
           #result_brownian <- length_scale^2 * min(x, fixed)

           return(result_centered_brownian)
         },
         uniform_centered_brownian = {
           # Implementation for the Brownian kernel
           # No parameter is required!

           u <- centering_param

           if (abs(u) >= abs(x) & abs(u) >= abs(fixed)) {
             result_centered_brownian <- length_scale / 2 * (-abs(x - fixed) +
                                                               (x^2 + fixed^2) / (2 * u) +
                                                               u / 3)
           } else if (abs(u) >= abs(x) & abs(u) < abs(fixed)) {
             result_centered_brownian <- length_scale / 2 * (abs(fixed) -
                                                               abs(fixed - x) +
                                                               x^2 / (2 * u) - u / 6)
           } else if (abs(u) < abs(x) & abs(u) >= abs(fixed)) {
             result_centered_brownian <- length_scale / 2 * (abs(x) -
                                                               abs(x - fixed) +
                                                               fixed^2 / (2 * u) - u / 6)
           } else if (abs(u) < abs(x) & abs(u) < abs(fixed)) {
             result_centered_brownian <- length_scale / 2 * (abs(x) +
                                                               abs(fixed) -
                                                               abs(x - fixed) -
                                                               2 * u / 3)
           }

           return(result_centered_brownian)
         },
         vectorized_uniform_centered_brownian = {
           # Implementation for the Brownian kernel
           # No parameter is required!

           u <- centering_param

           result_centered_brownian <- if_else(
             abs(u) >= abs(x) & abs(u) >= abs(fixed),
             length_scale / 2 * (-abs(x - fixed) + (x^2 + fixed^2) / (2 * u) + u / 3),
             if_else(
               abs(u) >= abs(x) & abs(u) < abs(fixed),
               length_scale / 2 * (abs(fixed) - abs(fixed - x) + x^2 / (2 * u) - u / 6),
               if_else(
                 abs(u) < abs(x) & abs(u) >= abs(fixed),
                 length_scale / 2 * (abs(x) - abs(x - fixed) + fixed^2 / (2 * u) - u / 6),
                 length_scale / 2 * (abs(x) + abs(fixed) - abs(x - fixed) - 2 * u / 3)
               )
             )
           )

           return(result_centered_brownian)
         },
         {
           stop("Invalid kernel type. Supported types: uniform, linear, polynomial,
                rq, rbf, gibbs, periodic, ou, matern, brownian,
                uniform_centered_brownian, gaussian_centered_brownian")
         }
  )
}
