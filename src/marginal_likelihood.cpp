#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>  // For parallel computing
#endif
#include "density.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

// Function to compute the marginal log likelihood for a pair of hyperparameters lambda and tau
// [[Rcpp::export]]
double marginal_log_likelihood(
    const arma::mat& centered_kernel_mat_at_sampled,
    const arma::vec& sampled_x,
    double min_x, double max_x,
    arma::vec p_vec,  // Ensure p_vec is passed by value, not reference
    double lambda, double tau,
    const arma::mat& std_rnorm_matrix,
    int MC_iterations,
    bool parallel_computing = true) {

  int n = centered_kernel_mat_at_sampled.n_rows;

  //  Print dimensions of input matrices
  //Rcpp::Rcout << "Dim(centered_kernel_mat_at_sampled): "
  //            << centered_kernel_mat_at_sampled.n_rows << " x " << centered_kernel_mat_at_sampled.n_cols << std::endl;
  //Rcpp::Rcout << "Dim(sampled_x): " << sampled_x.n_elem << " (vector)" << std::endl;
  //Rcpp::Rcout << "Dim(std_rnorm_matrix): "
  //            << std_rnorm_matrix.n_rows << " x " << std_rnorm_matrix.n_cols << std::endl;
  //Rcpp::Rcout << "MC_iterations: " << MC_iterations << ", n: " << n << std::endl;

  arma::mat w_sampled(MC_iterations, n);
  p_vec = arma::vectorise(p_vec);  // Ensure column vector

  //  Print dimensions of transformed p_vec
  //Rcpp::Rcout << "Dim(p_vec): " << p_vec.n_elem << " (vector)" << std::endl;

  // Parallel Monte Carlo weight sampling
#pragma omp parallel for if(parallel_computing) num_threads(omp_get_max_threads())
  for (int i = 0; i < MC_iterations; i++) {
    for (int j = 0; j < n; j++) {
      w_sampled(i, j) = std_rnorm_matrix(i, j) * std::sqrt(p_vec(j) / tau); //std::sqrt(p_vec(j) / tau)
    }
  }

  //  Print dimensions of generated w_sampled
  //Rcpp::Rcout << "Dim(w_sampled): " << w_sampled.n_rows << " x " << w_sampled.n_cols << std::endl;

  // Compute densities for sampled weights
  arma::mat densities_for_given_weights(MC_iterations, n);
  bool nan_detected = false;

#pragma omp parallel for if(parallel_computing) num_threads(omp_get_max_threads())
  for (int i = 0; i < MC_iterations; i++) {
    arma::vec dens_row = get_dens_wo_grid(
      centered_kernel_mat_at_sampled, min_x, max_x, sampled_x, lambda, w_sampled.row(i).t()
    );

    if (!dens_row.is_finite()) {
#pragma omp critical
{
  Rcpp::Rcout << "Error: NaN values in get_dens_wo_grid at iteration " << i << std::endl;
}
      nan_detected = true;
    }

    densities_for_given_weights.row(i) = dens_row.t();
  }

  if (nan_detected) {
    return -arma::datum::inf;
  }

  //  Print dimensions of densities_for_given_weights
  //Rcpp::Rcout << "Dim(densities_for_given_weights): "
  //           << densities_for_given_weights.n_rows << " x " << densities_for_given_weights.n_cols << std::endl;

  // Compute the likelihood vector by taking the mean across Monte Carlo iterations
  arma::vec likelihood_vector = arma::mean(densities_for_given_weights, 0).t();  // Transpose to column vector

  // Print dimensions of likelihood_vector
  //Rcpp::Rcout << "Dim(likelihood_vector): " << likelihood_vector.n_elem << " (vector)" << std::endl;

  // Check for NaN in likelihood_vector
  if (!likelihood_vector.is_finite()) {
    Rcpp::Rcout << "Error: NaN detected in likelihood_vector" << std::endl;
    return -arma::datum::inf;
  }

  // Compute the normalizing constant using trapezoidal integration
  double normalizing_cte = arma::as_scalar(arma::trapz(sampled_x, likelihood_vector));

  if (normalizing_cte == 0 || std::isnan(normalizing_cte) || std::isinf(normalizing_cte)) {
    Rcpp::Rcout << "Error: Normalizing constant is zero, NaN, or Inf" << std::endl;
    return -arma::datum::inf;
  }

  // Normalize the likelihood vector
  arma::vec norm_likelihood_vector = likelihood_vector / normalizing_cte;

  //  Print dimensions of normalized likelihood vector
  //Rcpp::Rcout << "Dim(norm_likelihood_vector): " << norm_likelihood_vector.n_elem << " (vector)" << std::endl;

  // Check for log(0) or negative values
  if (!norm_likelihood_vector.is_finite() || arma::any(norm_likelihood_vector <= 0)) {
    Rcpp::Rcout << "Error: norm_likelihood_vector contains zeros or negatives." << std::endl;
    return -arma::datum::inf;
  }

  // Compute the log of the marginal likelihood
  double marginal_log_likelihood_result = arma::sum(arma::log(norm_likelihood_vector));

  return marginal_log_likelihood_result;

}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

// Function to call R's `get_weights_wo_grid_BBsolve` from C++ using `.Call()`
// [[Rcpp::export]]
arma::vec call_get_weights_wo_grid_BBsolve(double lambda, double tau,
                                           const arma::mat& centered_kernel_mat_at_sampled,
                                           const arma::vec& sampled_x,
                                           double min_x, double max_x,
                                           const arma::vec& p_vec,
                                           bool print_trace = false) {
  // Calling the R function
  SEXP res = Rcpp::Function("get_weights_wo_grid_BBsolve")(lambda, tau,
                            centered_kernel_mat_at_sampled,
                            sampled_x, min_x, max_x,
                            p_vec, print_trace);
  return Rcpp::as<arma::vec>(res);
}

// Function to compute the marginal likelihood over a grid of hyperparameters
// [[Rcpp::export]]
arma::mat compute_marginal_likelihood_grid_parallel(
    const arma::mat& centered_kernel_mat_at_sampled,
    double min_x, double max_x,
    const arma::vec& sampled_x,
    const arma::mat& hyperparam_grid,
    double initial_lambda,
    const arma::vec& initial_w,
    int MC_iterations,
    int max_iterations,
    bool parallel_computing = true) {

  int n = sampled_x.n_elem;
  arma::mat results(hyperparam_grid.n_rows, 3);  // Store lambda, tau, and MLL
  arma::mat std_rnorm_matrix = arma::randn(MC_iterations, n);  // Monte Carlo samples

  // Initialize lambda, w_vec, and p_vec
  double lambda = initial_lambda;
  arma::vec w_vec = initial_w;

  arma::vec dens_vec = get_dens_wo_grid(centered_kernel_mat_at_sampled, min_x, max_x,
                                        sampled_x, lambda, w_vec);
  arma::vec p_vec = dens_vec / arma::sum(dens_vec);

  double max_mll = -arma::datum::inf;

  // Iterative optimization loop
  for (int t = 1; t <= max_iterations; t++) {
    arma::mat local_results(hyperparam_grid.n_rows, 3);

    // Parallel grid search over hyperparameters
#pragma omp parallel for if(parallel_computing) num_threads(omp_get_max_threads())
    for (size_t i = 0; i < hyperparam_grid.n_rows; i++) {
      double candidate_lambda = hyperparam_grid(i, 0);
      double candidate_tau = hyperparam_grid(i, 1);

      double mll = marginal_log_likelihood(
        centered_kernel_mat_at_sampled, sampled_x, min_x, max_x,
        p_vec, candidate_lambda, candidate_tau, std_rnorm_matrix, MC_iterations, parallel_computing
      );

      local_results(i, 0) = candidate_lambda;
      local_results(i, 1) = candidate_tau;
      local_results(i, 2) = mll;
    }
    results = local_results;
    // Find the best hyperparameters
    arma::uword best_idx;
    max_mll = local_results.col(2).max(best_idx);
    lambda = local_results(best_idx, 0);
    double tau = local_results(best_idx, 1);

    // Call the R function to update weights using BBsolve
    w_vec = call_get_weights_wo_grid_BBsolve(lambda, tau, centered_kernel_mat_at_sampled,
                                             sampled_x, min_x, max_x, p_vec, false);
    dens_vec = get_dens_wo_grid(centered_kernel_mat_at_sampled, min_x, max_x,
                                sampled_x, lambda, w_vec);
    p_vec = dens_vec / arma::sum(dens_vec);

    // Print iteration results
    Rcpp::Rcout << "\nIteration: " << t << ", lambda_hat: " << lambda
                << ", tau_hat: " << tau << ", max_mll: " << max_mll << std::endl;


  }


  return results;
}
