#ifndef MATRIXNORM_H
#define MATRIXNORM_H
#include "shared.h"
#include <RcppArmadillo.h>
#define DEBUG_MATRIXNORM false

inline arma::mat rmatrixnorm(const arma::mat &M, const arma::mat &U,
                             const arma::mat &V) {
  int n = M.n_rows;
  int p = M.n_cols;

  arma::mat X = Rcpp::rnorm(n * p);
  X.set_size(n, p);

  mat A = arma::chol(U, "lower");
  mat B = arma::chol(V, "upper");

  if (DEBUG_MATRIXNORM) {
    Rcpp::Rcout << "A = " << std::endl;
    A.print(Rcpp::Rcout);
    Rcpp::Rcout << "B = " << std::endl;
    B.print(Rcpp::Rcout);
  }

  return M + A * X * B;
}

inline double dmatrixnorm(const arma::mat &x, const arma::mat &M,
                          const arma::mat &U, const arma::mat &V,
                          const bool log_p = false) {
  arma::uword n = x.n_rows;
  arma::uword p = x.n_cols;
  double det_U = arma::det(U);
  double det_V = arma::det(V);
  arma::mat U_inv = arma::inv_sympd(U);
  arma::mat V_inv = arma::inv_sympd(V);

  arma::mat dev = x - M; // Taille (n x p)

  double exponent = -0.5 * arma::accu((U_inv * dev) % (dev * V_inv));

  double out;

  double log_norm_const = -0.5 * n * p * log(2 * M_PI) - 0.5 * n * log(det_V) -
                          0.5 * p * log(det_U);

  if (log_p) {
    out = log_norm_const + exponent;
  } else {
    out = std::exp(log_norm_const + exponent);
  }
  return out;
}

#endif
