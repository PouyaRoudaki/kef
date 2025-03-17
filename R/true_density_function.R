#' Compute True Density Based on Specified Distribution
#'
#' @param grid A numeric vector representing the grid points where the density should be computed.
#' @param density_characterization A list containing the density type and its parameters.
#'
#' @return A numeric vector containing the density values at the grid points.
#'
#' @examples
#' grid <- seq(-5, 5, length.out = 100)
#'
#' # Mixture Normal Example
#' params_mixture <- list(means = c(-2, 2), sds = c(1, 1.5), weights = c(0.4, 0.6))
#' density_characterization_mixture <- list(type = "mixture_normal", parameters = params_mixture)
#' true_density_mixture <- true_density_function(grid, density_characterization_mixture)
#' plot(grid, true_density_mixture, type = "l", main = "Mixture Normal Density", xlab = "x", ylab = "Density")
#'
#' # Uniform Example
#' params_uniform <- list(min = -3, max = 3)
#' density_characterization_uniform <- list(type = "uniform", parameters = params_uniform)
#' true_density_uniform <- true_density_function(grid, density_characterization_uniform)
#' plot(grid, true_density_uniform, type = "l", main = "Uniform Density", xlab = "x", ylab = "Density")
#'
#' # Poisson Example
#' grid_pos <- 0:10  # Poisson works with integer values
#' params_poisson <- list(lambda = 4)
#' density_characterization_poisson <- list(type = "poisson", parameters = params_poisson)
#' true_density_poisson <- true_density_function(grid_pos, density_characterization_poisson)
#' plot(grid_pos, true_density_poisson, type = "h", main = "Poisson Density", xlab = "x", ylab = "Density")
#'
#' # Beta Example
#' grid_beta <- seq(0, 1, length.out = 100)  # Beta is defined on [0,1]
#' params_beta <- list(shape1 = 2, shape2 = 5)
#' density_characterization_beta <- list(type = "beta", parameters = params_beta)
#' true_density_beta <- true_density_function(grid_beta, density_characterization_beta)
#' plot(grid_beta, true_density_beta, type = "l", main = "Beta Density", xlab = "x", ylab = "Density")
#'
#' @export
true_density_function <- function(grid, density_characterization) {

  x_grid <- grid  # Ensure the function uses the provided grid
  params <- density_characterization$parameters  # Extract parameters

  if (density_characterization$type == "mixture_normal") {
    means <- params$means
    sds <- params$sds
    mixture_weights <- params$weights

    if (length(means) != length(sds) || length(means) != length(mixture_weights)) {
      stop("Mismatch in mixture_normal parameters: means, sds, and weights must have the same length")
    }

    # Define a matrix of normal densities for each mean and standard deviation
    norm_density_matrix <- sapply(seq_along(means), function(i) {
      dnorm(x_grid, mean = means[i], sd = sds[i])
    })

    # Calculate the true density by taking the weighted sum of the columns
    true_density <- rowSums(norm_density_matrix %*% mixture_weights)

  } else if (density_characterization$type == "uniform") {
    unif_min <- params$min
    unif_max <- params$max
    true_density <- dunif(x_grid, min = unif_min, max = unif_max)

  } else if (density_characterization$type == "poisson") {
    lambda_poisson <- params$lambda
    true_density <- dpois(x_grid, lambda = lambda_poisson)

  } else if (density_characterization$type == "beta") {
    shape1_beta <- params$shape1
    shape2_beta <- params$shape2
    true_density <- dbeta(x_grid, shape1 = shape1_beta, shape2 = shape2_beta)
  } else {
    stop("Unsupported density type")
  }

  return(true_density)
}
