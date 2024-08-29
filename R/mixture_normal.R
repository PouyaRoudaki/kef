#' Generate Samples from a Normal Mixture Model
#'
#' This function generates random samples from a mixture of normal distributions.
#' The user specifies the number of samples, the means, standard deviations,
#' and mixture weights for each component of the mixture.
#'
#' @param n An integer specifying the number of samples to generate.
#' @param means A numeric vector of length 4 specifying the means of the normal components.
#' @param sds A numeric vector of length 4 specifying the standard deviations of the normal components.
#' @param mixture_weights A numeric vector of length 4 specifying the mixture weights for the components.
#'        The weights should sum to 1.
#'
#' @return A numeric vector of length `n` containing the generated samples from the normal mixture model.
#' @export
#'
#' @examples
#' # Example usage:
#' n <- 1000
#' means <- c(0, 5, 10, 15)
#' sds <- c(1, 1, 2, 2)
#' mixture_weights <- c(0.25, 0.25, 0.25, 0.25)
#' samples <- normal_mixture(n, means, sds, mixture_weights)
#' hist(samples, breaks = 30, main = "Histogram of Normal Mixture Samples")
normal_mixture <- function(n, means, sds, mixture_weights) {
  # Ensure that the mixture weights sum to 1 (optional check)
  if (sum(mixture_weights) != 1) {
    stop("The mixture weights must sum to 1.")
  }

  # Sample from the categorical distribution to choose the component for each sample
  component <- sample(1:length(mixture_weights), size = n, replace = TRUE, prob = mixture_weights)

  # Initialize the vector to store the generated samples
  samples <- numeric(n)

  # Generate the samples from the corresponding normal distributions
  for (i in 1:length(mixture_weights)) {
    samples[component == i] <- rnorm(sum(component == i), mean = means[i], sd = sds[i])
  }

  return(samples)
}


# Example: Generate 1000 samples from the mixture distribution
#set.seed(123) # For reproducibility
#samples <- normal_mixture(1000, means, sds, mixture_weights)

# Plot the histogram of the samples
#hist(samples, breaks = 50, main = "Histogram of Mixture Distribution Samples", xlab = "Value", probability = TRUE)

# Add density lines for each of the component distributions
#curve(mixture_weights[1] * dnorm(x, mean = means[1], sd = sds[1]), add = TRUE, col = "red")
#curve(mixture_weights[2] * dnorm(x, mean = means[2], sd = sds[2]), add = TRUE, col = "blue")
#curve(mixture_weights[3] * dnorm(x, mean = means[3], sd = sds[3]), add = TRUE, col = "green")
#curve(mixture_weights[4] * dnorm(x, mean = means[4], sd = sds[4]), add = TRUE, col = "purple")

#install.packages("ks")
#library(ks)

# Perform KDE with adaptive bandwidth
#kde_result <- kde(samples)

# Plot the KDE result
#plot(kde_result)


# Plot the histogram of the samples
#hist(samples, breaks = 50, probability = TRUE, main = "Histogram and KDE of Mixture Distribution Samples", xlab = "Value")

# Add density lines for each of the component distributions
#curve((mixture_weights[1] * dnorm(x, mean = means[1], sd = sds[1]) +
#         mixture_weights[2] * dnorm(x, mean = means[2], sd = sds[2]) +
#         mixture_weights[3] * dnorm(x, mean = means[3], sd = sds[3]) +
#         mixture_weights[4] * dnorm(x, mean = means[4], sd = sds[4])) , add = TRUE, col = "red",lwd = 2)


# Overlay the KDE plot on the histogram
#lines(kde_result$eval.points, kde_result$estimate, col = "blue", lwd  =2)

#lines(x_grid, probs$grid_x, col = "orange", lwd  =2)
