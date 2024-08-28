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


