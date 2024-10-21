#' Computes the d_{i, lambda} vector for the i-th observation
#'
#' @param gram Matrix representing the Gram matrix.
#' @param regularized_inv The regularized inverse of the Gram matrix.
#' @param i Index of the observation.
#'
#' @return A numeric value representing the d_{i, lambda}
#' @export
d_i_lambda <- function(gram, regularized_inv, i) {
  n <- nrow(gram)
  k_i <- gram[, i]
  e_i <- rep(0, n)
  e_i[i] <- 1
  result <- as.numeric(t(k_i) %*% regularized_inv %*% e_i)
  return(result)
}

#' Computes the p_{i, lambda} vector for the i-th observation
#'
#' @param gram Matrix representing the Gram matrix.
#' @param regularized_inv The regularized inverse of the Gram matrix.
#' @param i Index of the observation.
#'
#' @return A vector representing the p_{i, lambda}
#' @export
p_i_lambda <- function(gram, regularized_inv, i) {
  n <- nrow(gram)
  k_i <- gram[, i]
  e_i <- rep(0, n)
  e_i[i] <- 1
  d_i <- d_i_lambda(gram, regularized_inv, i)
  result <- e_i %*% t(k_i) %*% regularized_inv / (1 - d_i)
  return(result)
}

#' Computes the c_{i, lambda} vector for the i-th observation
#'
#' @param gram Matrix representing the Gram matrix.
#' @param regularized_inv The regularized inverse of the Gram matrix.
#' @param i Index of the observation.
#'
#' @return A vector representing the c_{i, lambda}
#' @export
c_i_lambda <- function(gram, regularized_inv, i) {
  n <- nrow(gram)
  k_i <- gram[, i]
  e_i <- rep(0, n)
  e_i[i] <- 1
  p_i <- p_i_lambda(gram, regularized_inv, i)
  one <- rep(1, n)
  result <- gram %*% one - gram[, i] - e_i %*% t(k_i) %*% one + e_i * gram[i, i] +
    p_i %*% gram %*% one - p_i %*% gram[, i] - p_i %*% e_i %*% t(k_i) %*% one +
    p_i %*% e_i * gram[i, i]

  return(result)
}

#' Calculates the A matrix
#'
#' @param gram Matrix representing the Gram matrix.
#' @param regularized_inv The regularized inverse of the Gram matrix.
#'
#' @return The A matrix used in shrinkage estimation.
#' @export
A_lambda <- function(gram, regularized_inv) {
  n <- nrow(gram)

  c_i_vectors <- sapply(1:n, function(j) {
    c_i_lambda(gram = gram, regularized_inv = regularized_inv, i = j)
  })

  result <- Reduce(`+`, lapply(1:n, function(i) c_i_vectors[, i] %*% t(c_i_vectors[, i]))) / ((n - 1)^2)
  return(result)
}

#' Calculates the B matrix
#'
#' @param gram Matrix representing the Gram matrix.
#' @param regularized_inv The regularized inverse of the Gram matrix.
#'
#' @return The B matrix used in shrinkage estimation.
#' @export
B_lambda <- function(gram, regularized_inv) {
  n <- nrow(gram)

  c_i_vectors <- sapply(1:n, function(j) {
    c_i_lambda(gram = gram, regularized_inv = regularized_inv, i = j)
  })

  result <- Reduce(`+`, lapply(1:n, function(i) c_i_vectors[, i] %*% t(gram[, i]))) / (n - 1)
  return(result)
}

#' Leave-One-Out Cross-Validation (LOOCV) error for shrinkage estimator
#'
#' @param gram Matrix representing the Gram matrix.
#' @param lambda The shrinkage parameter.
#'
#' @return The LOOCV error estimate.
#' @export
loocv_shr <- function(gram, lambda, precomp_reg_inv = NULL) {
  n <- nrow(gram)

  lambda_n <- lambda * (n - 1)

  if(is.null(precomp_reg_inv)){
    regularized_inverse <- solve(gram + lambda_n * diag(n))
  } else {
    regularized_inverse <- precomp_reg_inv
  }


  A <- A_lambda(gram, regularized_inverse)
  B <- B_lambda(gram, regularized_inverse)

  result <- as.numeric(
    1 / n * sum(diag(regularized_inverse %*% gram %*% regularized_inverse %*% A)) -
      2 / n * sum(diag(regularized_inverse %*% B)) +
      1 / n * sum(diag(gram))
  )
  print("-")
  return(result)
}

