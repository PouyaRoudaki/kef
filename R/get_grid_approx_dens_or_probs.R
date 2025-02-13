#' Get Middle Points of Grid
#'
#' This function calculates the midpoints between adjacent elements of a sample
#' and includes the minimum and maximum grid boundaries.
#'
#' @param min The minimum value of the grid.
#' @param samples A numeric vector representing the sample points within the grid.
#' @param max The maximum value of the grid.
#'
#' @return A numeric vector containing the minimum value, the midpoints of the sample, and the maximum value.
#' @export
#'
#' @examples
#' get_middle_points_grid(0, c(2, 4, 6), 8)
get_middle_points_grid_R <- function(min, samples, max){
  n <- length(samples)
  return(c( min , (samples[-n] + samples[-1])/2, max) )
}

#' Get Start and End Midpoints of Interval
#'
#' This function finds the start and end midpoints of an interval within a grid,
#' based on the sample midpoints.
#'
#' @param sample_mid_points A numeric vector of midpoints between the sample points.
#' @param grid A numeric vector representing the grid.
#' @param s_grid The start index of the grid interval.
#' @param e_grid The end index of the grid interval.
#'
#' @return A numeric vector containing the start and end midpoints of the interval.
#' @export
#'
#' @examples
#' get_s_e_mid_points_of_interval(c(1, 3, 5), c(0, 2, 4, 6), 2, 3)
get_s_e_mid_points_of_interval <- function(sample_mid_points, grid, s_grid, e_grid){
  s_mid_point <- max(sample_mid_points[which(sample_mid_points <= grid[s_grid])])
  e_mid_point <- min(sample_mid_points[which(sample_mid_points >= grid[e_grid])])
  return( c(s_mid_point,e_mid_point))
}

#' Get Corresponding Grids and Their Weights
#'
#' This function determines the grid intervals corresponding to each sample point
#' and calculates the weights for each grid point.
#'
#' @param sample A numeric vector representing the sample points.
#' @param grid A numeric vector representing the grid.
#'
#' @return A list containing two elements:
#' \item{grid_weights}{A data frame with the frequency and weight for each grid point.}
#' \item{corresponding_grid}{A data frame with the details of corresponding grid intervals.}
#' @export
#'
#' @examples
#' get_corresponding_grids_and_their_weights(c(1, 3, 5), c(0, 2, 4, 6))
get_corresponding_grids_and_their_weights <- function(sample, grid){

  # Find the middle points of samples
  sample_mid_points <- get_middle_points_grid(grid[1], sample, grid[length(grid)])

  # Find the interval positions within the grid
  mid_sample_positions_in_grid <- findInterval(sample_mid_points, grid)

  grid_idx <- vector()
  close_samples_flag <- vector()
  s_mid_sample_interval <- vector()
  e_mid_sample_interval <- vector()
  interval_length <- vector()

  # Loop over each sample point to determine the grid intervals and other properties
  for (i in 1:length(sample)) {

    if (i == 1) {
      s = 0
    } else {
      s = 1
    }

    grid_idx[i] <- toString(c((mid_sample_positions_in_grid[i]+s) : mid_sample_positions_in_grid[i+1]))

    if ((mid_sample_positions_in_grid[i]+s) <= mid_sample_positions_in_grid[i+1]) {
      close_samples_flag[i] <- FALSE
      s_mid_sample_interval[i] <- get_s_e_mid_points_of_interval(sample_mid_points,
                                                                 grid,
                                                                 mid_sample_positions_in_grid[i]+s,
                                                                 mid_sample_positions_in_grid[i+1])[1]
      e_mid_sample_interval[i] <- get_s_e_mid_points_of_interval(sample_mid_points,
                                                                 grid,
                                                                 mid_sample_positions_in_grid[i]+s,
                                                                 mid_sample_positions_in_grid[i+1])[2]
    } else {
      s_mid_sample_interval[i] <- grid[mid_sample_positions_in_grid[i+1]]
      e_mid_sample_interval[i] <- grid[(mid_sample_positions_in_grid[i]+s)]
      close_samples_flag[i] <- TRUE
    }

    interval_length[i] <- e_mid_sample_interval[i] - s_mid_sample_interval[i]
  }

  corresponding_grid_df <- data.frame(matrix(nrow = length(sample), ncol = 6))

  colnames(corresponding_grid_df) <- c("sample_idx", "grid", "close_sample_flag",
                                       "s_mid_sample_interval", "e_mid_sample_interval",
                                       "interval_length")

  corresponding_grid_df$sample_idx <- c(1:length(sample))

  corresponding_grid_df$grid <- grid_idx

  corresponding_grid_df$close_sample_flag <- close_samples_flag

  corresponding_grid_df$s_mid_sample_interval <- s_mid_sample_interval

  corresponding_grid_df$e_mid_sample_interval <- e_mid_sample_interval

  corresponding_grid_df$interval_length <- interval_length


  # Split the string into individual elements
  grid_index_with_duplication <- strsplit(toString(grid_idx), ",")[[1]]

  # Convert the string vector to a numeric vector
  grid_index_freq_of_duplication <- as.numeric(grid_index_with_duplication)

  # Create a data frame with the frequency of each grid index
  grid_index_freq_of_dup_df <- as.data.frame(table(grid_index_freq_of_duplication))

  # Calculate the weights for each grid point
  grid_index_freq_of_dup_df$weight <- (grid_index_freq_of_dup_df$Freq)^(-1)

  result <- list()

  result$grid_weights <- grid_index_freq_of_dup_df

  result$corresponding_grid <- corresponding_grid_df

  return(result)
}


#' Get Grid Approximate Densities or Probabilities
#'
#' This function computes the approximate densities or probabilities for each
#' sample point based on a grid and a list of densities.
#'
#' @param sample A numeric vector representing the sample points.
#' @param grid A numeric vector representing the grid.
#' @param dens_list A list containing the densities or probabilities at each grid point.
#'
#' @return A list containing two elements:
#' \item{prob}{A numeric vector of approximate probabilities for each sample point.}
#' \item{dens}{A numeric vector of approximate densities for each sample point.}
#' @export
#'
#' @examples
#' get_grid_approx_dens_or_probs(c(1, 3, 5), c(0, 2, 4, 6), list(grid_x = c(0.1, 0.4, 0.3, 0.2)))
get_grid_approx_dens_or_probs <- function(sample, grid, dens_list){

  list_of_corresponding_grid <- get_corresponding_grids_and_their_weights(sample, grid)

  q <- dens_list$grid_x

  # Normalize the densities or probabilities
  q_normalized <- q/sum(q)

  approx_prob <- numeric()

  # Loop over each sample to calculate the approximate probabilities
  for (i in 1:length(sample)) {
    approx_prob[i] <- 0
    corresponding_grid_index <- as.numeric(strsplit(list_of_corresponding_grid$corresponding_grid$grid[i], ",")[[1]])
    for (j in corresponding_grid_index) {
      approx_prob[i] <- approx_prob[i] + list_of_corresponding_grid$grid_weights$weight[j] * q_normalized[j]
    }
  }

  # Calculate the approximate densities based on interval lengths
  interval_lengths <- list_of_corresponding_grid$corresponding_grid$interval_length
  approx_dens <- approx_prob / interval_lengths

  approx_dens_or_prob <- list()

  approx_dens_or_prob$prob <- approx_prob / sum(approx_prob)
  approx_dens_or_prob$dens <- approx_dens

  return(approx_dens_or_prob)
}




#' Compute Approximate Probabilities and Densities on a Grid
#'
#' This function calculates approximate probabilities and densities for given sample points
#' based on a specified grid and a list of densities or probabilities on that grid.
#'
#' @param sample A numeric vector of sample points for which the probabilities or densities need to be approximated.
#' @param grid A numeric vector representing the grid points on which the densities or probabilities are defined.
#' @param dens_list A list containing the grid probabilities or densities.
#'        It must include an element `grid_x`, which is a numeric vector of probabilities or densities corresponding to the grid points.
#'
#' @return A list containing two elements:
#' \describe{
#'   \item{`prob`}{A numeric vector of approximate probabilities for the sample points, normalized to sum to 1.}
#'   \item{`dens`}{A numeric vector of approximate densities for the sample points, based on interval lengths.}
#' }
#'
#' @export
#'
#' @examples
#' # Example 1: Compute probabilities and densities
#' sample <- c(1, 3, 5)
#' grid <- c(0, 2, 4, 6)
#' dens_list <- list(grid_x = c(0.1, 0.4, 0.3, 0.2))
#'
#' result <- get_grid_approx_dens_or_probs(sample, grid, dens_list)
#'
#' # Access the approximate probabilities
#' print(result$prob)
#'
#' # Access the approximate densities
#' print(result$dens)
get_grid_approx_dens_or_probs_vectorized <- function(sample, grid, dens_list) {
  # Step 1: Get corresponding grid information
  list_of_corresponding_grid <- get_corresponding_grids_and_their_weights(sample, grid)

  # Step 2: Normalize the densities or probabilities
  q <- dens_list$grid_x
  q_normalized <- q / sum(q)

  # Step 3: Parse the grid indices for all samples
  grid_indices <- lapply(list_of_corresponding_grid$corresponding_grid$grid, function(x) as.numeric(unlist(strsplit(x, ","))))

  # Step 4: Compute weights for all grid indices
  weights <- list_of_corresponding_grid$grid_weights
  weight_map <- setNames(weights$weight, weights$grid_index_freq_of_duplication)

  # Vectorize probability calculation
  approx_prob <- sapply(grid_indices, function(indices) sum(weight_map[indices] * q_normalized[indices]))

  # Step 5: Calculate densities based on interval lengths
  interval_lengths <- list_of_corresponding_grid$corresponding_grid$interval_length
  approx_dens <- approx_prob / interval_lengths

  # Step 6: Normalize probabilities
  approx_prob_normalized <- approx_prob / sum(approx_prob)

  # Step 7: Construct the result list
  approx_dens_or_prob <- list(
    prob = approx_prob_normalized,
    dens = approx_dens
  )

  return(approx_dens_or_prob)
}

