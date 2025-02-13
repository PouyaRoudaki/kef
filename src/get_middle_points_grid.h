#ifndef GET_MIDDLE_POINTS_GRID_H
#define GET_MIDDLE_POINTS_GRID_H

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// Function declaration
arma::vec get_middle_points_grid(double min, const arma::vec& samples, double max);

#endif // GET_MIDDLE_POINTS_GRID_H
