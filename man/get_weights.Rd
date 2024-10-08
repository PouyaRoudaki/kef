% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_weights.R
\name{get_weights}
\alias{get_weights}
\title{Estimate Weights Using the Newton-Raphson Method}
\usage{
get_weights(
  lambda_hat,
  tau_hat,
  centered_kernel_mat_at_sampled,
  centered_kernel_mat_at_grid,
  centered_kernel_self_grid,
  sampled_x,
  x_grid,
  type_of_p_is_prob = TRUE,
  type_of_q_is_prob = TRUE,
  method_of_p_calculation = "ordinary"
)
}
\arguments{
\item{lambda_hat}{A scalar representing the estimated lambda parameter.}

\item{tau_hat}{A scalar representing the estimated tau parameter.}

\item{centered_kernel_mat_at_sampled}{A square matrix (n x n) representing the centered kernel matrix
evaluated at the sampled points.}

\item{centered_kernel_mat_at_grid}{A matrix (n x m) representing the centered kernel matrix evaluated
at the grid points, where n is the number of sampled points and m is the number of grid points.}

\item{centered_kernel_self_grid}{A vector of length m representing the diagonal of the centered kernel
matrix evaluated at the grid points.}

\item{x_grid}{A vector representing the grid points.}
}
\value{
A numeric vector of length n representing the estimated weight vector.
}
\description{
This function estimates the weight vector using an iterative Newton-Raphson method.
The weights are updated based on the provided kernel matrices, regularization parameters,
and grid points.
}
\examples{
# Example usage:
lambda_hat <- 1
tau_hat <- 0.5
n <- 10
m <- 20
centered_kernel_mat_at_sampled <- matrix(runif(n * n), n, n)
centered_kernel_mat_at_grid <- matrix(runif(n * m), n, m)
centered_kernel_self_grid <- runif(m)
x_grid <- seq(0, 1, length.out = m)
weights <- get_weights(lambda_hat, tau_hat, centered_kernel_mat_at_sampled,
                       centered_kernel_mat_at_grid, centered_kernel_self_grid, x_grid)
}
