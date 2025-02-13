#ifndef DENSITY_H  // Prevents multiple inclusions
#define DENSITY_H

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// Function declaration (prototype)
arma::vec get_dens_wo_grid(const arma::mat& centered_kernel_mat_at_sampled,
                           double min_x,
                           double max_x,
                           const arma::vec& sampled_x,
                           double lambda_hat,
                           const arma::vec& weight_hat_vec);

#endif  // DENSITY_H
