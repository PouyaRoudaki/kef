// src/centered_kernel_matrix_hd.cpp
#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::export]]
arma::mat centered_kernel_matrix_hd(const arma::mat& first_mat_kernel,
                                    const arma::mat& second_mat_kernel,
                                    const arma::mat& centering_mat_grid,
                                    double hurst_coef) {
  int n0 = first_mat_kernel.n_rows;
  int n1 = second_mat_kernel.n_rows;
  int n2 = centering_mat_grid.n_rows;

  arma::mat term1_matrix(n0, n1);
  arma::mat term2_matrix(n0, n2);
  arma::mat term3_matrix(n1, n2);
  arma::mat term4_matrix(n2, n2);

  // Term1: pairwise between first_mat_kernel and second_mat_kernel
#pragma omp parallel for collapse(2)
  for (int i = 0; i < n0; i++) {
    for (int j = 0; j < n1; j++) {
      term1_matrix(i, j) = std::pow(arma::norm(first_mat_kernel.row(i) - second_mat_kernel.row(j), 2), 2 * hurst_coef);
    }
  }

  // Term2: first_mat_kernel vs centering_mat_grid
#pragma omp parallel for collapse(2)
  for (int i = 0; i < n0; i++) {
    for (int j = 0; j < n2; j++) {
      term2_matrix(i, j) = std::pow(arma::norm(first_mat_kernel.row(i) - centering_mat_grid.row(j), 2), 2 * hurst_coef);
    }
  }

  // Term3: second_mat_kernel vs centering_mat_grid
#pragma omp parallel for collapse(2)
  for (int i = 0; i < n1; i++) {
    for (int j = 0; j < n2; j++) {
      term3_matrix(i, j) = std::pow(arma::norm(second_mat_kernel.row(i) - centering_mat_grid.row(j), 2), 2 * hurst_coef);
    }
  }

  // Term4: centering_mat_grid vs centering_mat_grid
#pragma omp parallel for collapse(2)
  for (int i = 0; i < n2; i++) {
    for (int j = 0; j < n2; j++) {
      term4_matrix(i, j) = std::pow(arma::norm(centering_mat_grid.row(i) - centering_mat_grid.row(j), 2), 2 * hurst_coef);
    }
  }

  arma::vec term2_means = arma::mean(term2_matrix, 1);  // n0 × 1
  arma::vec term3_means = arma::mean(term3_matrix, 1);  // n1 × 1
  double term4_mean = arma::mean(arma::vectorise(term4_matrix));

  arma::mat result_matrix = -0.5 * (term1_matrix -
    arma::repmat(term2_means, 1, n1) -
    arma::repmat(term3_means.t(), n0, 1) +
    term4_mean);

  return result_matrix;
}
