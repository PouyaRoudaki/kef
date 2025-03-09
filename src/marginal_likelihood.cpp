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



  arma::mat w_sampled(MC_iterations, n); // Make the w_sampled matrix
  p_vec = arma::vectorise(p_vec);  // Ensure column vector

  // Parallel Monte Carlo weight sampling
#pragma omp parallel for if(parallel_computing) num_threads(omp_get_max_threads())
  for (int i = 0; i < MC_iterations; i++) {
    for (int j = 0; j < n; j++) {
      w_sampled(i, j) = std_rnorm_matrix(i, j) * std::sqrt(p_vec(j) / tau); //std::sqrt(p_vec(j) / tau)
    }
  }

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

  // Compute the product along each column (equivalent to apply(..., 2, prod) in R)
  arma::vec monte_carlo_likelihood = arma::prod(densities_for_given_weights, 1);


  //Rcpp::Rcout << "p_vec size: " << monte_carlo_likelihood.n_rows << " x " << monte_carlo_likelihood.n_cols << std::endl;

  // Compute the marginal likelihood (mean of all joint conditional densities)
  double marginal_likelihood = arma::mean(monte_carlo_likelihood);

  double log_marginal_likelihood = std::log(marginal_likelihood);

  // Handle -inf case
  if (!std::isfinite(log_marginal_likelihood)) {
    log_marginal_likelihood = -10000;
  }

  //Rcpp::Rcout << "-----------------------------------------------------" << std::endl;
  //Rcpp::Rcout << "lambda: " << lambda << ", tau: " << tau << ", mll: " << log_marginal_likelihood <<std::endl;
  //Rcpp::Rcout << "-----------------------------------------------------" << std::endl;

  // Return the marginal_log_likelihood
  return log_marginal_likelihood;
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
