#include <RcppArmadillo.h>
#include "density.h"
#include "get_middle_points_grid.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// Compute the function s(weight_hat_vec)
// [[Rcpp::export]]
arma::vec get_s_function(const arma::vec& weight_hat_vec,
                         double lambda_hat,
                         double tau_hat,
                         const arma::mat& centered_kernel_mat_at_sampled,
                         const arma::vec& sampled_x,
                         double min_x,
                         double max_x,
                         Rcpp::Nullable<arma::vec> prior_variance_p_vector = R_NilValue) {

  // Compute densities
  arma::vec dens = get_dens_wo_grid(centered_kernel_mat_at_sampled, min_x, max_x, sampled_x, lambda_hat, weight_hat_vec);

  // Compute base measure
  arma::vec sample_mid_points = get_middle_points_grid(min_x, sampled_x, max_x);
  arma::vec base_measure_weights = sample_mid_points.subvec(1, sample_mid_points.n_elem - 1) -
    sample_mid_points.subvec(0, sample_mid_points.n_elem - 2);

  // Compute density-based probabilities
  arma::vec dens_sampled_base = dens % base_measure_weights;
  arma::vec prob_sampled_base = dens_sampled_base / sum(dens_sampled_base);

  // Check if prior_variance_p_vector is provided
  arma::vec prob_sampled;
  if (prior_variance_p_vector.isNotNull()) {
    prob_sampled = Rcpp::as<arma::vec>(prior_variance_p_vector);
  } else {
    prob_sampled = dens / sum(dens);
  }

  // Compute the function s(weight_hat_vec)
  arma::vec s = lambda_hat * (sum(centered_kernel_mat_at_sampled, 1) -
    centered_kernel_mat_at_sampled * prob_sampled_base * sampled_x.n_elem) -
    tau_hat * (weight_hat_vec / prob_sampled);

  return s;
}

