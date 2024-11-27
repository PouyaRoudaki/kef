#' Convolution of a Kernel with Itself over a Dominated Finite Measure
#'
#' This function computes the convolution of a kernel with itself over a specified dominated finite measure.
#' It supports various kernel functions and dominated measures, defaulting to "rbf" for the kernel and "gaussian"
#' for the measure. The function allows for customization of kernel and measure parameters.
#'
#' @param x Numeric vector, the first input to the kernel function.
#' @param y Numeric vector, the second input to the kernel function.
#' @param kernel_type Character string specifying the kernel function type (default: "rbf").
#' @param length_scale       Kernel parameters, length scale
#' @param degree  Kernel parameters, polynomial degree
#' @param free_add  Kernel parameters, free additive parameter
#' @param free_mult   Kernel parameters, free multiplicative parameter
#' @param nu_matern  Kernel parameters, nu_matern for matern kernel
#' @param dom_measure_type Character stype = "2"tring indicating the dominated finite measure for convolution (default: "gaussian").
#' @param dom_measure_eta List of parameters for the dominated finite measure, with 'eta' as an example.
#'
#' @return The convolution of the kernel with itself under the specified dominated finite measure.
#' @export
#'
#' @author Pouya Roudaki \email{roudaki.pouya@@gmail.com}.
#' @examples
#' r_conv_kernel(x, y, "rbf", list(length_scale = 1), "gaussian", list(eta = 1))
#'


r_conv_kernel <- function(x,
                          y,
                          kernel_type = "rbf",
                          length_scale = 1,
                          degree = 2,
                          free_add = 0,
                          free_mult = 1,
                          nu_matern = 1,
                          dom_measure_type = "gaussian",
                          dom_measure_eta = 1 ) {
  # Combine kernel and measure types for case selection.
  combined_key <- paste(tolower(kernel_type), tolower(dom_measure_type), sep = "_")

  # Calculate dimensionality from the input vector.
  D <- length(x)
  # Square the length_scale parameter for kernel calculations.
  theta <- (length_scale)^2

  # Perform computation based on kernel and measure types.
  switch(combined_key,
         "rbf_lebesgue" = {
           # Convolution for RBF kernel with Lebesgue measure.
           result <- (pi)^(D/2) * theta^(D/2) * exp(-1/(4*theta) * (norm(x-y)^2))
           return(result)
         },
         "rbf_gaussian" = {
           # Access eta parameter from dominated measure parameters.
           eta <- dom_measure_eta
           # Convolution for RBF kernel with Gaussian measure.
           result <- (2*pi*eta^2)^(-0.5) * (2*pi)^(D/2) * (2/theta + 1/(eta^2))^(-D/2) *
             exp(-1/(4*theta) * (norm(x-y)^2)) *
             exp(-1/2 * (1/2 * theta + eta^2)^(-1) * (norm((x+y)/2)^2))
           return(result)
         },
         {
           # Error handling for unknown kernel-measure combinations.
           stop("Unknown combination of kernel_type and dom_measure_type")
         })
}
