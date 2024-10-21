#' Calculate Betas for Different Estimator Types
#'
#' This function calculates the coefficients (betas) for each fold in a cross-validation
#' setup, based on the specified estimator type. It supports standard, regularized,
#' and shrinkage estimators by applying operations over each fold using `lapply`.
#'
#' @param gram A square symmetric matrix representing the Gram matrix of the data.
#' @param folds A list of folds, each containing a train and test dataset with ids.
#' @param lambda A regularization parameter used for regularized and shrinkage estimators.
#' @param estimator_type A character string indicating the type of estimator: "standard",
#'        "regularized", or "shrinkage".
#'
#' @return A list of calculated betas for each fold.
#' @export
beta_calculator <- function(gram, folds, lambda, estimator_type) {

  # Determine the number of observations
  n <- nrow(gram)

  # Calculate betas based on the estimator type for each training fold
  betas <- lapply(folds, function(fold) {

    # Extract length of train data
    train_length <- nrow(fold$train)

    if(estimator_type == "standard") {

      # return beta for "i"th partition
      return(rep(1 / train_length, train_length))

    } else if(estimator_type == "regularized") {

      # return beta for "i"th partition
      return(rep(1 / train_length, train_length) / (1 + lambda))

    } else if(estimator_type == "shrinkage") {

      # Extract train IDs from the fold
      train_ids <- fold$train$id

      # Extract gram matrix for the train data
      train_gram <- gram[train_ids, train_ids]

      # Extract regularized inverse matrix for the train data
      inverse_regularizer <- solve(train_gram + train_length * lambda * diag(train_length))

      # 1_n vector
      one_n <- rep(1, train_length) / train_length

      # return beta for "i"th partition
      return(inverse_regularizer %*% train_gram %*% one_n)
    }
  })

  return(betas)
}

#' Calculate Cross Validation Error for Each Fold
#'
#' This function calculates the cross-validation error for a given fold of data
#' by computing the sum of three summands based on the provided 'betas' and a
#' precomputed Gram matrix 'gram' indexed by train and test IDs.
#'
#' @param fold A list containing training and testing subsets identified by 'id'.
#' @param betas A numeric vector of coefficients to be applied during calculation.
#'
#' @return A single numeric value representing the cross-validation error for the fold.
#' @export
#'
#' @examples
#' # Assuming 'gram' and 'betas' are predefined and a 'fold' structure is available
#' # process_fold(fold, beta)
process_fold <- function(gram,fold, beta) {


  # Extract train and test IDs from the fold
  train_ids <- fold$train$id  # train_ids
  test_ids <- fold$test$id    # test_ids

  # Create a vector of ones for the test IDs
  one_test <- rep(1, length(test_ids))

  # Calculate the first summand
  first_summand <- as.numeric(1 / (length(test_ids)^2) * t(one_test) %*% gram[test_ids, test_ids] %*% one_test)

  # Calculate the second summand
  second_summand <- as.numeric(-2 / (length(test_ids)) * t(beta) %*% gram[train_ids, test_ids] %*% one_test)

  # Calculate the third summand
  third_summand <- as.numeric( t(beta) %*% gram[train_ids, train_ids] %*% beta)

  # Sum of the three summands represents the cross-validation error
  summands_sum <- first_summand + second_summand + third_summand

  return(summands_sum)
}

#' Perform Cross Validation
#'
#' This function performs cross-validation over a specified number of folds.
#' It divides the data into train and test sets, applies the 'process_fold'
#' function to each, and calculates the average error across all folds. Leave-One-Out
#' Cross Validation when number of folds is equal to the dimension of gram matrix.
#'
#' @param gram A square symmetric matrix representing the Gram matrix of the data.
#' @param betas A numeric vector of coefficients to be applied during calculation.
#' @param folds_number The number of folds to divide the data into for cross-validation. Defaults to the number of rows in 'gram'.
#'
#' @return A single numeric value representing the average cross-validation error.
#' @export
#'
#' @examples
#' # Assuming 'gram', 'betas', and an appropriate 'folds_number' are defined
#' # cross_validation(gram, betas, folds_number)
cross_validation <- function(gram, lambda, folds_number = nrow(gram), estimator_type = "standard") {

  # Total number of observations
  n <- nrow(gram)

  # Generate random fold assignments for each row in 'gram'
  fold_indices <- sample(rep(1:folds_number, length.out = n))

  # Create a dataframe mapping each row to its fold
  data <- data.frame(id = 1:n, fold = fold_indices)

  # Split data into folds
  folds <- lapply(1:folds_number, function(x) {
    list(
      train = data[data$fold != x, ],
      test = data[data$fold == x, ]
    )
  })

  betas <- beta_calculator(gram, folds, lambda, estimator_type)

  # Calculate the cross-validation error for each fold
  errors <- sapply(1:folds_number, function(i){
    process_fold(gram,folds[[i]], betas[[i]])})

  # Calculate and return the average error across all folds
  final_result <- mean(errors)

  print("-")
  return(final_result)
}



