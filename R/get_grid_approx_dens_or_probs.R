get_middle_points_grid <- function(min, samples, max){
  n <- length(samples)
  return(c( min , (samples[-n] + samples[-1])/2, max) )
}

get_s_e_mid_points_of_interval <- function(sample_mid_points, grid, s_grid, e_grid){
  s_mid_point <- max(sample_mid_points[which(sample_mid_points <= grid[s_grid])])
  e_mid_point <- min(sample_mid_points[which(sample_mid_points >= grid[e_grid])])
  return( c(s_mid_point,e_mid_point))
}

get_corresponding_grids_and_their_weights <- function(sample,grid){

  # Find the middle points of samples
  sample_mid_points <- get_middle_points_grid(grid[1],sample,grid[length(grid)])

  # Find the interval positions
  mid_sample_positions_in_grid <- findInterval(sample_mid_points, grid)

  grid_idx <- vector()
  close_samples_flag <- vector()
  s_mid_sample_interval <- vector()
  e_mid_sample_interval <- vector()
  interval_length <- vector()
  for (i in 1:length(sample)) {

    if( i == 1){
      s = 0
    }else{
      s = 1
    }

    grid_idx[i] <- toString(c((mid_sample_positions_in_grid[i]+s) : mid_sample_positions_in_grid[i+1]))

    if((mid_sample_positions_in_grid[i]+s) <= mid_sample_positions_in_grid[i+1]) {
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
      s_mid_sample_interval[i] <-  grid[mid_sample_positions_in_grid[i+1]]
      e_mid_sample_interval[i] <- grid[(mid_sample_positions_in_grid[i]+s)]
      close_samples_flag[i] <- TRUE
    }

    interval_length[i] <- e_mid_sample_interval[i] - s_mid_sample_interval[i]
  }

  corresponding_grid_df <- data.frame(matrix(nrow = length(sample),ncol = 6))

  colnames(corresponding_grid_df) <- c("sample_idx","grid","close_sample_flag",
                                       "s_mid_sample_interval","e_mid_sample_interval",
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

  # Print the numeric vector
  grid_index_freq_of_dup_df <- as.data.frame(table(grid_index_freq_of_duplication))

  grid_index_freq_of_dup_df$weight <- (grid_index_freq_of_dup_df$Freq)^(-1)

  result <- list()

  result$grid_weights <-  grid_index_freq_of_dup_df

  result$corresponding_grid <- corresponding_grid_df

  return(result)
}




get_grid_approx_dens_or_probs <- function(sample,grid, dens_list){

  list_of_corresponding_grid <- get_corresponding_grids_and_their_weights(sample,grid)

  q <- dens_list$grid_x

  q_normalized <- q/sum(q)

  approx_prob <- numeric()
  for (i in 1:length(sample)) {
    approx_prob[i] <- 0
    corresponding_grid_index <- as.numeric(strsplit(list_of_corresponding_grid$corresponding_grid$grid[i],",")[[1]])
    for (j in corresponding_grid_index) {
      approx_prob[i] <- approx_prob[i] + list_of_corresponding_grid$grid_weights$weight[j] * q_normalized[j]
    }
  }

  interval_lengths <- list_of_corresponding_grid$corresponding_grid$interval_length
  approx_dens <- approx_prob / interval_lengths

  approx_dens_or_prob <- list()

  approx_dens_or_prob$prob <- approx_prob
  approx_dens_or_prob$dens <- approx_dens

  return(approx_dens_or_prob)
}
