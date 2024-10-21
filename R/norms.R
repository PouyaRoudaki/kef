#' Vector Norm Calculation
#'
#' This function calculates the vector norm of a numeric or complex vector.
#'
#' @param x Numeric or complex vector.
#' @param p Numeric value specifying the norm order (default is 2).
#'
#' @return The calculated vector norm.
#'
#' @export
#'
#' @examples
#' norm(c(3, 4)) # Calculates the Euclidean norm of the vector c(3, 4)
#' norm(c(1, 2, 3), p = 3) # Calculates the L3 norm of the vector c(1, 2, 3)
#'
#' @seealso
#' \code{\link{abs}}, \code{\link{sum}}, \code{\link{max}}, \code{\link{min}}
#'
#' @references
#' - Book: "Numerical Analysis" by Richard L. Burden and J. Douglas Faires.
#'
#' @author Pouya Roudaki \email{roudaki.pouya@@gmail.com}.
#'
#' @export
#'
#' @keywords norm, vector

norm <- function(x, p=2) {
  tryCatch({
    if (!is.numeric(x) && !is.complex(x))
      stop("Input 'x' must be numeric or complex.")

    if (!is.numeric(p) || length(p) != 1)
      stop("Input 'p' must be a numeric value with length 1.")

    if (p > -Inf && p < Inf) sum(abs(x)^p)^(1/p)
    else if (p ==  Inf) max(abs(x))
    else if (p == -Inf) min(abs(x))
    else return(NULL)
  }, error = function(e) {
    cat("Error: ", conditionMessage(e), "\n")
  })
}
