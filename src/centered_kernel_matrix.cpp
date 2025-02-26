// src/centered_kernel_matrix.cpp
#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::export]]
arma::mat centered_kernel_matrix(const arma::vec& first_vec_kernel,
                                 const arma::vec& second_vec_kernel,
                                 const arma::vec& centering_grid,
                                 double hurst_coef) {
  int n0 = first_vec_kernel.n_elem;
  int n1 = second_vec_kernel.n_elem;
  int n2 = centering_grid.n_elem;

  arma::mat term1_matrix(n0, n1);
  arma::mat term2_matrix(n0, n2);
  arma::mat term3_matrix(n1, n2);
  arma::mat term4_matrix(n2, n2);

  // Parallel computation of term1_matrix
#pragma omp parallel for collapse(2)
  for (int i = 0; i < n0; i++) {
    for (int j = 0; j < n1; j++) {
      term1_matrix(i, j) = std::pow(std::abs(first_vec_kernel(i) - second_vec_kernel(j)), 2 * hurst_coef);
    }
  }

  // Parallel computation of term2_matrix
#pragma omp parallel for collapse(2)
  for (int i = 0; i < n0; i++) {
    for (int j = 0; j < n2; j++) {
      term2_matrix(i, j) = std::pow(std::abs(first_vec_kernel(i) - centering_grid(j)), 2 * hurst_coef);
    }
  }

  // Parallel computation of term3_matrix
#pragma omp parallel for collapse(2)
  for (int i = 0; i < n1; i++) {
    for (int j = 0; j < n2; j++) {
      term3_matrix(i, j) = std::pow(std::abs(second_vec_kernel(i) - centering_grid(j)), 2 * hurst_coef);
    }
  }

  // Parallel computation of term4_matrix
#pragma omp parallel for collapse(2)
  for (int i = 0; i < n2; i++) {
    for (int j = 0; j < n2; j++) {
      term4_matrix(i, j) = std::pow(std::abs(centering_grid(i) - centering_grid(j)), 2 * hurst_coef);
    }
  }

  // Compute mean values (these parts don't need parallelization)
  arma::vec term2_means = arma::mean(term2_matrix, 1);
  arma::vec term3_means = arma::mean(term3_matrix, 1);
  double term4_mean = arma::mean(arma::vectorise(term4_matrix));

  // Compute final matrix using broadcasting
  arma::mat result_matrix = -0.5 * (term1_matrix - arma::repmat(term2_means, 1, n1) -
    arma::repmat(term3_means.t(), n0, 1) + term4_mean);

  return result_matrix;
}

