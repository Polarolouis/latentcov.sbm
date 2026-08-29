#include "matrixnorm.h"
#include "shared.h"
#include <RcppArmadillo.h>

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