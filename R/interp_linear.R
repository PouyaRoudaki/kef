#' Fast Linear Interpolation and Extrapolation with Automatic Method Selection
#'
#' This function performs linear interpolation and extrapolation using either the
#' base R `approx()` function or a fast C++ implementation via Rcpp, depending
#' on the size of the dataset.
#' It automatically selects the most efficient method based on a threshold.
#'
#' @param x Numeric vector of x-values (must be sorted in ascending order).
#' @param y Numeric vector of corresponding y-values.
#' @param xnew Numeric vector of new x-values where interpolation is required.
#' @param threshold Integer specifying the dataset size threshold to switch
#'        from R’s `approx()` to the Rcpp-based implementation. Default is 10,000.
#'
#' @return A numeric vector of interpolated y-values corresponding to `xnew`.
#' @export
#'
#' @examples
#' # Example with small dataset (uses approx())
#' x <- 1:10
#' y <- (1:10)^2
#' xnew <- c(2.5, 4.5, 6.5)
#' interp_linear(x, y, xnew)
#'
#' # Example with large dataset (uses Rcpp)
#' set.seed(42)
#' x <- sort(runif(1e5, 0, 100))
#' y <- sin(x)
#' xnew <- sort(runif(1e4, 0, 100))
#' interp_linear(x, y, xnew)
interp_linear <- function(x, y, xnew, threshold = 10000) {
  # Check for duplicate x values
  if (any(duplicated(x))) {
    #warning("Duplicate values in 'x' detected. Removing duplicates and keeping the first occurrence.")
    unique_indices <- !duplicated(x)  # Keep only the first occurrence
    x <- x[unique_indices]
    y <- y[unique_indices]
  }
  if (length(x) > threshold) {
    # Use Rcpp for large datasets
    return(as.vector(interp_linear_cpp(x, y, xnew))) # O(log(n))
  } else {
    # Use R's approx() for small datasets
    return(approx(x, y, xout = xnew, method = "linear", rule = 2)$y) # O(n) rule 2 considers extrapolation as well
  }
}
