#include "get_middle_points_grid.h"

// [[Rcpp::export]]
arma::vec get_middle_points_grid(double min, const arma::vec& samples, double max) {
  int n = samples.n_elem;
  arma::vec midpoints(n + 1);

  midpoints[0] = min;
  for (int i = 1; i < n; i++) {
    midpoints[i] = (samples[i - 1] + samples[i]) / 2.0;
  }
  midpoints[n] = max;

  return midpoints;
}
