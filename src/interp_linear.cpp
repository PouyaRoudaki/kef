#include <RcppArmadillo.h>
#include <algorithm> // for std::lower_bound
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::vec interp_linear_cpp(const arma::vec& x, const arma::vec& y, const arma::vec& xnew) {
  size_t n = x.n_elem;
  size_t m = xnew.n_elem;

  // Ensure valid input
  if (n == 0 || m == 0) {
    Rcpp::stop("Error: Input vectors must not be empty.");
  }
  if (x.n_elem != y.n_elem) {
    Rcpp::stop("Error: x and y must have the same length.");
  }

  arma::vec ynew(m);  // Output vector

  for (size_t i = 0; i < m; ++i) {
    auto it = std::lower_bound(x.begin(), x.end(), xnew[i]);
    size_t idx = std::distance(x.begin(), it);

    if (idx == 0) {
      // Left extrapolation (linear)
      double slope = (y[1] - y[0]) / (x[1] - x[0]);
      ynew[i] = y[0] + slope * (xnew[i] - x[0]);
    } else if (idx >= n) {
      // Right extrapolation (linear)
      double slope = (y[n-1] - y[n-2]) / (x[n-1] - x[n-2]);
      ynew[i] = y[n-1] + slope * (xnew[i] - x[n-1]);
    } else {
      // Standard linear interpolation
      double x1 = x[idx - 1], x2 = x[idx];
      double y1 = y[idx - 1], y2 = y[idx];
      ynew[i] = y1 + (y2 - y1) * (xnew[i] - x1) / (x2 - x1);
    }
  }

  return ynew;
}
