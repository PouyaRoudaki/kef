#' Create a Matrix of Ones
#'
#' This function generates a matrix filled with ones, with dimensions specified by the user.
#'
#' @param n Integer. The number of rows in the matrix.
#' @param m Integer. The number of columns in the matrix.
#'
#' @return A matrix with `n` rows and `m` columns, where all elements are set to 1.
#' @export
#'
#' @examples
#' # Create a 3x4 matrix of ones
#' ones(3, 4)
ones <- function(n, m) {
  matrix(1, nrow = n, ncol = m)
}



#' Format Number with Three Significant Digits
#'
#' This function formats a given numeric value `x` to three significant digits.
#' It uses scientific notation if the value is less than 1e-3 or greater than 1e3.
#'
#' @param x A numeric value that you want to format.
#'
#' @return A character string representing the formatted number with three significant digits.
#' The output will be in scientific notation if the value of `x` is very small (less than 1e-3)
#' or very large (greater than 1e3).
#'
#' @export
#'
#' @examples
#' format_number(3.456)
#' # Output: "3.46"
#'
#' format_number(1.127e-12)
#' # Output: "1.13e-12"
#'
#' format_number(123456)
#' # Output: "1.23e+05"
#'
format_number <- function(x) {
  format(x, digits = 3, scientific = ifelse(x < 1e-3 | x > 1e3, TRUE, FALSE))
}
