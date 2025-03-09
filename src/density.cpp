#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// Density at sampled points
// [[Rcpp::export]]
arma::vec density_at_sampled_x(const arma::mat& centered_kernel_mat_at_sampled,
                                   double lambda_hat,
                                   const arma::vec& weight_hat_vec) {
  arma::vec diag_vals = centered_kernel_mat_at_sampled.diag();  // Extract diagonal elements
  arma::vec p_vec = exp(lambda_hat * (centered_kernel_mat_at_sampled.t() * weight_hat_vec - 0.5 * diag_vals));
  return p_vec;
}

// Density at grid points
// [[Rcpp::export]]
arma::vec density_at_grid(const arma::mat& centered_kernel_mat_at_grid,
                              const arma::vec& centered_kernel_self_grid,
                              double lambda_hat,
                              const arma::vec& weight_hat_vec) {
  arma::vec p_vec = exp(lambda_hat * (centered_kernel_mat_at_grid.t() * weight_hat_vec - 0.5 * centered_kernel_self_grid));
  return p_vec;
}

// Compute densities for sampled points without a grid
// [[Rcpp::export]]
arma::vec get_dens_wo_grid(const arma::mat& centered_kernel_mat_at_sampled,
                           double min_x,
                           double max_x,
                           const arma::vec& sampled_x,
                           double lambda_hat,
                           const arma::vec& weight_hat_vec) {

  arma::vec f_x = lambda_hat * (centered_kernel_mat_at_sampled.t() * weight_hat_vec - 0.5 * centered_kernel_mat_at_sampled.diag());
  double max_f_x = arma::max(f_x);  // Get max value

  // Set a safety threshold to detect overflow risk (700 is a reasonable limit)
  const double EXP_THRESHOLD = 700;

  arma::vec p_vec;

  if (max_f_x > EXP_THRESHOLD) {
    //Rcpp::Rcout << "Warning: Exponential function at risk of overflow. Applying log-scaling.\n";
    p_vec = arma::exp(f_x - max_f_x);  // Apply log-stabilization
  } else {
    p_vec = arma::exp(f_x);  // Compute normally if safe
  }

  // Compute normalizing constant safely using trapezoidal integration
  double normalizing_cte = arma::as_scalar(arma::trapz(sampled_x, p_vec));

  // Prevent division by zero
  if (normalizing_cte == 0 || std::isnan(normalizing_cte) || std::isinf(normalizing_cte)) {
    Rcpp::Rcout << "Error: Normalizing constant is zero, NaN, or Inf\n";
    return arma::vec(sampled_x.n_elem, arma::fill::zeros);  // Return a zero vector to avoid NaN propagation
  }

  return p_vec / normalizing_cte;  // Return normalized density
}


// Compute densities for sampled and grid points
// [[Rcpp::export]]
List get_dens_or_prob(const arma::mat& centered_kernel_mat_at_sampled,
                          const arma::mat& centered_kernel_mat_at_grid,
                          const arma::vec& centered_kernel_self_grid,
                          const arma::vec& sampled_x,
                          const arma::vec& x_grid,
                          double lambda_hat,
                          const arma::vec& weight_hat_vec,
                          bool type_of_p_is_prob = true,
                          bool type_of_q_is_prob = true) {

  arma::vec dens_sampled_x = density_at_sampled_x(centered_kernel_mat_at_sampled, lambda_hat, weight_hat_vec);
  arma::vec dens_grid = density_at_grid(centered_kernel_mat_at_grid, centered_kernel_self_grid, lambda_hat, weight_hat_vec);

  double normalizing_cte = arma::as_scalar(trapz(x_grid, dens_grid));
  arma::vec dens_sampled_x_norm = dens_sampled_x / normalizing_cte;
  arma::vec dens_grid_norm = dens_grid / normalizing_cte;

  arma::vec prob_sampled_x = dens_sampled_x_norm / sum(dens_sampled_x_norm);
  arma::vec prob_grid_x = dens_grid_norm / sum(dens_grid_norm);

  return List::create(
    Named("sampled_x") = (type_of_p_is_prob ? prob_sampled_x : dens_sampled_x_norm),
    Named("grid_x") = (type_of_q_is_prob ? prob_grid_x : dens_grid_norm)
  );
}
