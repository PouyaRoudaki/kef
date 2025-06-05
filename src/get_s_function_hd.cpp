#include <RcppArmadillo.h>
#include "density.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// Compute the function s(weight_hat_vec)
// [[Rcpp::export]]
arma::vec get_s_function_hd(const arma::vec& weight_hat_vec,
                            double lambda_hat,
                            double tau_hat,
                            const arma::mat& centered_kernel_mat_at_sampled,
                            const arma::mat& sampled_x,  // Now using matrix for general d-dim
                            Rcpp::Nullable<arma::vec> prior_variance_p_vector = R_NilValue,
                            bool with_prob_in_var = true,
                            bool normalised_q = true,
                            bool normalised_p = true,
                            bool p_with_base = false) {

  int n = sampled_x.n_rows;

  // Call the generate_voronoi() function defined in Rcpp
  Rcpp::Function generate_voronoi("generate_voronoi");
  Rcpp::NumericMatrix input_points = Rcpp::wrap(sampled_x); // arma::mat to R matrix

  Rcpp::List vor = generate_voronoi(input_points, -10, 10, -10, 10);
  arma::vec voronoi_weights = Rcpp::as<arma::vec>(vor["polygon_areas"]);

  // Compute densities
  arma::vec dens = get_dens_wo_grid(centered_kernel_mat_at_sampled, 0, 1, arma::vec(n, fill::zeros), lambda_hat, weight_hat_vec);

  // Base-weighted density q(x) * dμ(x)
  arma::vec dens_sampled_base = dens % voronoi_weights;

  // Normalize q
  arma::vec prob_sampled_base;
  if (normalised_q) {
    prob_sampled_base = dens_sampled_base / sum(dens_sampled_base);
  } else {
    prob_sampled_base = dens_sampled_base;
  }

  // Define p(x)
  arma::vec prob_sampled;
  if (prior_variance_p_vector.isNotNull()) {
    prob_sampled = Rcpp::as<arma::vec>(prior_variance_p_vector);
  } else if (normalised_p) {
    prob_sampled = dens / sum(dens);
  } else {
    prob_sampled = dens;
  }

  if (!with_prob_in_var) {
    prob_sampled = arma::ones<arma::vec>(n);
  }

  if (p_with_base) {
    prob_sampled = prob_sampled_base;
  }

  // Final score function
  arma::vec s = lambda_hat * (sum(centered_kernel_mat_at_sampled, 1) -
                                centered_kernel_mat_at_sampled * prob_sampled_base * n) -
    tau_hat * (weight_hat_vec / prob_sampled);

  return s;
}
