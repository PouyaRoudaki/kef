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
    arma::vec p_vec,
    double lambda, double tau,
    const arma::mat& std_rnorm_matrix,
    int MC_iterations,
    bool parallel_computing = true) {

  int n = centered_kernel_mat_at_sampled.n_rows;
  arma::mat w_sampled(MC_iterations, n);

  // Parallel Monte Carlo weight sampling
  #pragma omp parallel for if(parallel_computing) num_threads(omp_get_max_threads())
  for (int i = 0; i < MC_iterations; i++) {
    for (int j = 0; j < n; j++) {
      w_sampled(i, j) = std_rnorm_matrix(i, j) * (p_vec(j) / tau);
    }
  }

  // Compute densities for sampled weights
  arma::vec likelihood_vector(MC_iterations);

  #pragma omp parallel for if(parallel_computing) num_threads(omp_get_max_threads())
  for (int i = 0; i < MC_iterations; i++) {
    likelihood_vector(i) = arma::sum(get_dens_wo_grid(centered_kernel_mat_at_sampled,
                                     min_x, max_x, sampled_x,
                                     lambda, w_sampled.row(i).t()));
  }

  // Normalize likelihood using Monte Carlo integration
  double normalizing_cte = arma::sum(likelihood_vector) / MC_iterations;
  if (normalizing_cte == 0) return -arma::datum::inf;  // Prevent log(0)

  arma::vec norm_likelihood_vector = likelihood_vector / normalizing_cte;

  // Compute and return log marginal likelihood
  return arma::sum(arma::log(norm_likelihood_vector));
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

    results.row(t - 1) = local_results.row(best_idx);
  }

  return results;
}
