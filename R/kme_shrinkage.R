#' Calculate the length of a vector after scaling.
#'
#' @param a Numeric element of a vector
#' @param length_scale General numeric scaling factor (default is 1)
#'
#' @return Numeric scaling factor
#' @export
#'
#'
#' @author Pouya Roudaki \email{roudaki.pouya@@gmail.com}.

aux_gibbs_length <- function(a, length_scale = 1) {
  return(a * length_scale)
}

#' Calculate the exponential term in the Gibbs kernel function.
#'
#' @param a Numeric element of first vector
#' @param b Numeric element of second vector
#' @param length_scale Numeric scaling factor (default is 1)
#'
#' @return Numeric value representing the exponential term
#' @export
#'
#' @author Pouya Roudaki \email{roudaki.pouya@@gmail.com}.
#'
#'
aux_gibbs_exp_term <- function(a, b, length_scale = 1) {
  l_a <- aux_gibbs_length(a, length_scale)
  l_b <- aux_gibbs_length(b, length_scale)

  denominator <- l_a^2 + l_b^2

  result_aux_exp_term <- (a - b)^2 / denominator
  return(result_aux_exp_term)
}

#' Calculate the Gibbs kernel between two vectors.
#'
#' @param a Numeric element of first vector
#' @param b Numeric element of second vector
#' @param length_scale Numeric scaling factor (default is 1)
#' @param exp_term Numeric value representing the exponential term
#'
#' @return Numeric value representing the Gibbs kernel
#' @export
#'
#' @author Pouya Roudaki \email{roudaki.pouya@@gmail.com}.
#'
#'
aux_gibbs_kernel <- function(a, b, length_scale = 1, exp_term) {
  l_a <- aux_gibbs_length(a, length_scale)
  l_b <- aux_gibbs_length(b, length_scale)

  denominator <- l_a^2 + l_b^2

  result_aux_gibbs_kernel <- sqrt(2 * l_a * l_b / denominator) * exp_term
  return(result_aux_gibbs_kernel)
}

#' Convert upper triangle vector to matrix
#'
#' @param upper_triangle_vector a vector including the upper triangle elements of
#' a matrix
#'
#' @return symmetric_matrix which is a symmetric matrix made by upper triangle vector
#' @export
#'
#' @author Pouya Roudaki \email{roudaki.pouya@@gmail.com}.
#'
#'
upper_tri_to_mat <- function(upper_triangle_vector){

  # Calculate the size of the square matrix
  n <- floor((sqrt(8 * length(upper_triangle_vector)+1)+1)/2)

  # Create a symmetric matrix
  symmetric_matrix <- matrix(0, nrow = n, ncol = n)

  # Fill the lower triangle with the values from the vector
  symmetric_matrix[lower.tri(symmetric_matrix, diag = F)] <- upper_triangle_vector

  # Copy the upper triangle to the lower triangle to make it symmetric
  symmetric_matrix <- symmetric_matrix + t(symmetric_matrix)

  return(symmetric_matrix)

}


#' Generate Random Samples from a Mixture of Gaussian Distributions
#'
#' This function generates random samples from a specified mixture of Gaussian distributions.
#' It allows for the specification of the means, standard deviations, and mixture probabilities
#' for each component of the mixture.
#'
#' @param n Integer, the number of samples to generate.
#' @param means Numeric vector, the means of each Gaussian component in the mixture.
#' @param sds Numeric vector, the standard deviations of each Gaussian component in the mixture.
#' @param probabilities Numeric vector, the mixing probabilities for each component.
#'        The probabilities must sum to 1.
#'
#' @return A numeric vector of random samples drawn from the specified mixture of Gaussian distributions.
#' @export
#'
#' @examples
#' set.seed(123)
#' n = 100
#' means = c(-2, 2)
#' sds = c(1, 1.5)
#' probabilities = c(0.3, 0.7)
#' samples = rnorm_mixture(n, means, sds, probabilities)
#' hist(samples, breaks = 30, col = "skyblue", main = "Samples from Mixture of Gaussians")
rnorm_mixture <- function(n, means, sds, probabilities) {
  if(length(means) != length(sds) || length(means) != length(probabilities)) {
    stop("Lengths of means, sds, and probabilities must be equal")
  }

  if(abs(sum(probabilities) - 1) > .Machine$double.eps^0.5) {
    stop("Probabilities must sum to 1")
  }

  # Generate component assignments based on mixture probabilities
  components <- sample(1:length(probabilities), size = n, replace = TRUE, prob = probabilities)

  # Generate samples
  samples <- mapply(function(component, n) {
    rnorm(n, mean = means[component], sd = sds[component])
  }, component = components, n = rep(1, n))

  return(samples)
}


#' Generate Samples from a Normal Mixture Model with a Base Measure
#'
#' This function generates samples from a Gaussian mixture model with a base measure,
#' accounting for centering and lambda parameters.
#'
#' @param n Integer, number of samples to generate.
#' @param probabilities Numeric vector, probabilities for each component of the mixture.
#' @param means Numeric vector, means of each component of the mixture.
#' @param sds Numeric vector, standard deviations of each component of the mixture.
#' @param lambda Numeric, parameter for the base measure.
#' @param centering_param Numeric, centering parameter for the base measure.
#'
#' @return A numeric vector of generated samples.
#' @export
#'
#' @examples
#' samples <- rnorm_mixture_wrt_base_measure(1000, c(0.3, 0.7), c(-2, 2), c(1, 1.5), 0.5, 1)
rnorm_mixture_wrt_base_measure <- function(n, probabilities, means, sds, lambda, centering_param) {

  require(distr)

  # Function to compute the PDF for the first integral
  pdf_function_for_integral_1 <- function(x, probabilities, means, sds, lambda, centering_param) {
    base_measure_gaussian_mixture <- 0
    for (i in seq_along(probabilities)) {
      base_measure_gaussian_mixture <- base_measure_gaussian_mixture +
        probabilities[i] / (sqrt(2 * pi) * sds[i]) * exp(- (x - means[i])^2 / (2 * sds[i]^2)) *
        exp(-lambda / 4 * (x^2 / centering_param + centering_param / 3))
    }
    return(base_measure_gaussian_mixture)
  }

  # Function to compute the PDF for the second integral
  pdf_function_for_integral_2 <- function(x, probabilities, means, sds, lambda, centering_param) {
    base_measure_gaussian_mixture <- 0
    for (i in seq_along(probabilities)) {
      base_measure_gaussian_mixture <- base_measure_gaussian_mixture +
        probabilities[i] / (sqrt(2 * pi) * sds[i]) * exp(- (x - means[i])^2 / (2 * sds[i]^2)) *
        exp(-lambda / 4 * (2 * abs(x) - 2 * centering_param / 3))
    }
    return(base_measure_gaussian_mixture)
  }

  # Integrate the PDF functions to normalize
  Integral_of_pdf_1 <- integrate(function(x) pdf_function_for_integral_1(x, probabilities,
                                                                         means, sds, lambda, centering_param),
                                 subdivisions = 10000, rel.tol = 1e-10,
                                 abs.tol = 1e-10, lower = -centering_param, upper = centering_param)$value

  Integral_of_pdf_2 <- integrate(function(x) pdf_function_for_integral_2(x, probabilities,
                                                                         means, sds, lambda, centering_param),
                                 subdivisions = 10000, rel.tol = 1e-10,
                                 abs.tol = 1e-10, lower = -Inf, upper = -centering_param)$value

  Integral_of_pdf_3 <- integrate(function(x) pdf_function_for_integral_2(x, probabilities,
                                                                         means, sds, lambda, centering_param),
                                 subdivisions = 10000, rel.tol = 1e-10,
                                 abs.tol = 1e-10, lower = centering_param, upper = Inf)$value

  Integral_of_pdf <- Integral_of_pdf_1 + Integral_of_pdf_2 + Integral_of_pdf_3

  # Define the final PDF
  pdf <- function(x) {
    v_x <- exp(-0.5 * kernel(x, x, type = "uniform_centered_brownian", length_scale = lambda, centering_param = centering_param))
    p_x <- 0
    for (i in seq_along(probabilities)) {
      p_x <- p_x +
        probabilities[i] * 1 / (sqrt(2 * pi) * sds[i]) * exp(-(x - means[i])^2 / (2 * sds[i]^2))
    }
    return(v_x * p_x / Integral_of_pdf)
  }

  # Vectorize the PDF function for efficiency
  pdf_vectorized <- Vectorize(pdf)

  # Create the distribution object
  dist <- AbscontDistribution(d = pdf_vectorized)

  # Generate samples from the distribution
  samples <- r(dist)(n)

  return(samples)
}



# Function to initialize the Excel file
#' Initialize Excel File
#'
#' This function creates an Excel file with specified columns to store
#' simulation parameters and error rates.
#'
#' @param file_name A string specifying the name of the Excel file to be created.
#'
#' @return None
#' @export
#'
#' @examples
#' initialize_excel_file("simulation_results.xlsx")
initialize_excel_file <- function(file_name) {
  # Create a data frame with the required columns
  df <- data.frame(n_fixed_u = numeric(),
                   n_fixed_x = numeric(),
                   sigma_w = numeric(),
                   lambda =  numeric(),
                   sigma_w_hat_n_power = numeric(),
                   lambda_hat = numeric(),
                   centering_dist = character(),
                   centering_param = numeric(),
                   natural_param_error = numeric())

  # Write the data frame to an Excel file
  write.xlsx(df, file = file_name, sheetName = "Results", rowNames = FALSE)
}

# Function to append data to the Excel file
#' Append Data to Excel File
#'
#' This function appends new simulation data to an existing Excel file.
#'
#' @param file_name A string specifying the name of the Excel file to append data to.
#' @param data A data frame containing the new simulation data to be appended.
#'
#' @return None
#' @export
#'
#' @examples
#' new_data <- data.frame(n = 10, lambda = 0.5, sigma = 0.1, error_rate = 0.02)
#' append_to_excel("simulation_results.xlsx", new_data)
append_to_excel <- function(file_name, data) {
  # Load the existing workbook
  wb <- loadWorkbook(file_name)

  # Read the existing data
  existing_data <- read.xlsx(file_name, sheet = "Results")

  # Combine the existing data with the new data
  combined_data <- rbind(existing_data, data)

  # Write the combined data back to the Excel file
  writeData(wb, sheet = "Results", x = combined_data, startCol = 1, startRow = 1, rowNames = FALSE)

  # Save the workbook
  saveWorkbook(wb, file = file_name, overwrite = TRUE)
}



#' Check if a Matrix is Singular
#'
#' This function checks if a given matrix is singular by calculating its determinant.
#' If the determinant is less than a specified tolerance, the matrix is considered singular.
#'
#' @param mat A numeric matrix to be checked for singularity.
#' @param tol A numeric value specifying the tolerance level for the determinant.
#'            Default is `.Machine$double.eps`, the smallest positive floating-point number.
#'
#' @return A logical value: `TRUE` if the matrix is singular, `FALSE` otherwise.
#' @export
#'
#' @examples
#' mat <- matrix(c(1, 2, 2, 4), nrow = 2)
#' is.singular(mat)  # Should return TRUE since the matrix is singular
#'
is.singular <- function(mat, tol = .Machine$double.eps) {
  if (!is.matrix(mat)) {
    stop("Input must be a matrix.")
  }
  det(mat) < tol
}
