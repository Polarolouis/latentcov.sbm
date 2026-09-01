#include "utils.h"

//' Sample a Matrix Normal with parameters M, U and V
//'
//' @param M a matrix of size n,p the mean of the matrix normal
//' @param U a variance-covariance matrix for the rows of Y
//' @param V a variance-covariance matrix for the columns of Y
//'
//' @returns Y a matrix of size n,p with the above mean and covariances
// [[Rcpp::export]]
arma::mat rmatrixnormal_cpp(const arma::mat &M, const arma::mat &U,
                            const arma::mat &V) {
  return rmatrixnorm(M, U, V);
}

//' Compute density for X of a Matrix Normal with parameters M, U and V
//'
//' @param X a matrix of size n,p for which to compute the density
//' @param M a matrix of size n,p the mean of the matrix normal
//' @param U a variance-covariance matrix for the rows of Y
//' @param V a variance-covariance matrix for the columns of Y
//'
//' @returns the density a scalar
// [[Rcpp::export]]
double dmatrixnormal_cpp(const arma::mat &X, const arma::mat &M,
                         const arma::mat &U, const arma::mat &V,
                         bool log_p = false) {
  return dmatrixnorm(X, M, U, V, log_p);
}

// [[Rcpp::export]]
arma::mat compute_Uinv(const arma::mat &U) { return inv_sympd(U); }

//' Inverse pivot coordinate transformation (C++)
//'
//' Performs the inverse pivot coordinate transformation for the latent
//' covariance stochastic block model, mapping the \eqn{K-1} latent
//' coordinates of a row back to the \eqn{K}-simplex of membership
//' probabilities.
//' @param x Input matrix (\eqn{n \times K-1}) of latent coordinates
//' @param norm Normalization type ("orthogonal" or "orthonormal")
//' @param log if \code{TRUE}, the log-probabilities are returned
//' @return A matrix (\eqn{n \times K}) of (log) simplex probabilities
//'
//' @export
// [[Rcpp::export]]
arma::mat pivot_coord_inv(arma::mat &x, std::string norm = "orthonormal",
                          bool log = false) {
  // Mirror precisely the R implementation: x <- -x and then operate on that
  arma::mat xneg = -x;
  arma::mat xback;
  arma::mat y(x.n_rows, x.n_cols + 1, arma::fill::zeros);
  int D = x.n_cols + 1;
  double first_fill = 1.0;

  if (norm != "orthogonal" && norm != "orthonormal") {
    Rcpp::stop("Norm %s not implemented !", norm);
  }

  if (norm == "orthonormal") {
    first_fill = -std::sqrt((double)(D - 1) / (double)D);
  } else if (norm == "orthogonal") {
    first_fill = 1.0;
  }

  y.col(0) = first_fill * xneg.col(0);

  for (int i = 1; i < (int)y.n_cols; ++i) {
    for (int j = 0; j < i; ++j) {
      unsigned int ull_i = static_cast<unsigned int>(i);
      unsigned int ull_j = static_cast<unsigned int>(j);
      double denom = 1.0;
      if (norm == "orthonormal") {
        denom = std::sqrt((double)(D - ull_j) * (double)(D - ull_j - 1.0));
      }
      y.col(ull_i) += xneg.col(ull_j) / denom;
    }
  }

  for (int i = 1; i < (int)y.n_cols - 1; ++i) {
    unsigned int ull_i = static_cast<unsigned int>(i);
    double multip = 1.0;
    if (norm == "orthonormal") {
      multip = std::sqrt((double)(D - i - 1) / (double)(D - i));
    }
    y.col(ull_i) -= xneg.col(ull_i) * multip;
  }

  arma::vec max_rows = arma::max(y, 1);
  arma::mat yexp = arma::exp(y - arma::repmat(max_rows, 1, y.n_cols));
  xback = yexp / arma::repmat(arma::sum(yexp, 1), 1, y.n_cols);

  if (log) {
    xback = arma::log(xback);
  }

  return xback;
}