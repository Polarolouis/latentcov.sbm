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

// [[Rcpp::export]]
arma::mat default_Psi_function_cpp(int K) {
  arma::mat Psi(K - 1, K, arma::fill::zeros);
  for (int j = 0; j < K - 1; ++j) {
    double coef =
        std::sqrt(static_cast<double>(K - j - 1) / static_cast<double>(K - j));
    for (int c = j + 1; c < K; ++c) {
      Psi(j, c) = coef / static_cast<double>(K - j - 1);
    }
    Psi(j, j) = -coef;
  }
  return Psi;
}

//' Inverse pivot coordinate transformation (C++)
//'
//' Performs the inverse pivot coordinate transformation for the latent
//' covariance stochastic block model, mapping the \eqn{K-1} latent
//' coordinates of a row back to the \eqn{K}-simplex of membership
//' probabilities.
//' @param x Input matrix (\eqn{n \times K-1}) of latent coordinates
//' @param basis The basis to use
//' @param log if \code{TRUE}, the log-probabilities are returned
//' @return A matrix (\eqn{n \times K}) of (log) simplex probabilities
//'
//' @export
// [[Rcpp::export]]
arma::mat ilrInv_cpp(const arma::mat &z, const arma::mat &basis,
                     bool log = false) {
  arma::mat clr = z * basis;

  arma::vec max_clr = arma::max(clr, 1);
  arma::mat xexp = arma::exp(clr - arma::repmat(max_clr, 1, clr.n_cols));

  arma::mat x = xexp / arma::repmat(arma::sum(xexp, 1), 1, xexp.n_cols);

  if (log) {
    x = arma::log(x);
  }

  return x;
}