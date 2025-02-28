% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/conv_gram_matrix.R
\name{conv_gram_matrix}
\alias{conv_gram_matrix}
\title{Compute Convolutional Gram matrix for a list of vectors using a specified kernel function.}
\usage{
conv_gram_matrix(
  vec_list,
  kernel_type = "rbf",
  kernel_params = list(length_scale = 1, degree = 2, free_add = 0, free_mult = 1,
    nu_matern = 1),
  dom_measure_type = "gaussian",
  dom_measure_parmas = list(eta = 1)
)
}
\arguments{
\item{vec_list}{A list of vectors for which the Gram matrix needs to be computed.}

\item{kernel_type}{The type of kernel function to be used (default: "rbf").}

\item{kernel_params}{A list of parameters for the kernel function (default values provided).}

\item{dom_measure_type}{Character string indicating the dominated finite measure for convolution (default: "gaussian").}

\item{dom_measure_parmas}{List of parameters for the dominated finite measure, with 'eta' as an example.}
}
\value{
A Convolutional Gram matrix computed using the specified kernel function.
}
\description{
Compute Convolutional Gram matrix for a list of vectors using a specified kernel function.
}
\examples{
vectors_list <- list(c(1, 2, 3), c(4, 5, 6), c(7, 8, 9))
result_conv_gram_matrix <- conv_gram_matrix(vectors_list, kernel_type = "linear")
print(result_conv_gram_matrix)
}
