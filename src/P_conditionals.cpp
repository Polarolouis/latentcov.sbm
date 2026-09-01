#include "P_conditionals.h"

//' Mean of Pi given P minus i and Sigma (C++)
//'
//' Computes the mean of the conditional distribution of a single row of P
//' given the rest of P, Sigma, and sigma2. This is the C++ version of
//' \code{cond_Pi_given_P_min_i_sigma}.
//' @param P a matrix of size \eqn{n_1 \times K-1} specifying latent
//'   positions
//' @param Sigma a covariance matrix between the row nodes
//' @param sigma2 variance parameter indicating the variance between the
//'   K-1 columns of P
//' @param i 0-based index of the row to compute the conditional for
//' @return the conditional mean, a row vector of length \eqn{K-1}
// [[Rcpp::export]]
arma::rowvec mean_of_Pi_given_P_min_i_sigma(arma::mat &P, arma::mat &Sigma,
                                            double sigma2, int i) {
  int n = P.n_rows;

  // Indices minus_i = setdiff(1:n, i)  (0-based: tous sauf i)
  arma::uvec minus_i = arma::regspace<arma::uvec>(0, n - 1);
  minus_i.shed_row(i);

  // Si = Sigma[i, minus_i] %*% solve(Sigma[minus_i, minus_i])
  arma::rowvec Sigma_i_minus = Sigma.row(i);
  Sigma_i_minus.shed_col(i); // Sigma[i, minus_i]
  arma::mat Sigma_minus_minus = Sigma.submat(minus_i, minus_i);
  // Resolve the linear system directly instead of computing the explicit
  // inverse: more stable numerically and closer to R's solve().
  arma::vec Si_t = arma::solve(Sigma_minus_minus, Sigma_i_minus.t());
  arma::rowvec Si = Si_t.t();

  // mean_i = Si %*% P[minus_i, ]
  arma::rowvec mean_i = Si * P.rows(minus_i);
  return mean_i;
}

//' Covariance of Pi given P minus i and Sigma (C++)
//'
//' Computes the covariance of the conditional distribution of a single row
//' of P given the rest of P, Sigma, and sigma2. This is the C++ version of
//' \code{cond_Pi_given_P_min_i_sigma}.
//' @param P a matrix of size \eqn{n_1 \times K-1} specifying latent
//'   positions
//' @param Sigma a covariance matrix between the row nodes
//' @param sigma2 variance parameter indicating the variance between the
//'   K-1 columns of P
//' @param i 0-based index of the row to compute the conditional for
//' @return the conditional covariance, a matrix of size
//'   \eqn{K-1 \times K-1}
// [[Rcpp::export]]
arma::mat cov_of_Pi_given_P_min_i_sigma(arma::mat &P, arma::mat &Sigma,
                                        double sigma2, int i) {
  int n = P.n_rows;
  int K_minus_1 = P.n_cols;

  // Indices minus_i = setdiff(1:n, i)  (0-based: tous sauf i)
  arma::uvec minus_i = arma::regspace<arma::uvec>(0, n - 1);
  minus_i.shed_row(i);

  // Si = Sigma[i, minus_i] %*% solve(Sigma[minus_i, minus_i])
  arma::rowvec Sigma_i_minus = Sigma.row(i);
  Sigma_i_minus.shed_col(i); // Sigma[i, minus_i]
  arma::mat Sigma_minus_minus = Sigma.submat(minus_i, minus_i);
  // Resolve the linear system directly instead of computing the explicit
  // inverse: more stable numerically and closer to R's solve().
  arma::vec Si_t = arma::solve(Sigma_minus_minus, Sigma_i_minus.t());
  arma::rowvec Si = Si_t.t();

  // cov_i = sigma2 * (Sigma[i,i] - Si %*% Sigma[-i, i]) * I
  arma::vec Sigma_minus_i = Sigma.col(i);
  Sigma_minus_i.shed_row(i); // Sigma[-i, i]
  double scalar = Sigma(i, i) - arma::as_scalar(Si * Sigma_minus_i);
  arma::mat cov_i =
      (sigma2 * scalar) * arma::eye<arma::mat>(K_minus_1, K_minus_1);

  return cov_i;
}

//' Mean of P_b given P minus b and Sigma (C++)
//'
//' Computes the mean of the conditional distribution of a block of P
//' given the rest of P, Sigma, and sigma2. This is the C++ version of
//' \code{cond_Pi_given_P_min_i_sigma}.
//' @param P a matrix of size \eqn{n_1 \times K-1} specifying latent
//'   positions
//' @param Sigma a covariance matrix between the row nodes
//' @param sigma2 variance parameter indicating the variance between the
//'   K-1 columns of P
//' @param indiv_indicees 0-based indices of the rows to compute the conditional
// for ' @return the conditional mean, a row vector of length \eqn{K-1}
// [[Rcpp::export]]
arma::mat block_mean_of_Pi_given_P_min_i_sigma(arma::mat &P, arma::mat &Sigma,
                                               double sigma2,
                                               arma::uvec indiv_indices) {
  int n = P.n_rows;

  // Indices minus_i = setdiff(1:n, i)  (0-based: tous sauf i)
  arma::uvec minus_i = arma::regspace<arma::uvec>(0, n - 1);
  minus_i.shed_rows(indiv_indices);

  // Si = Sigma[i, minus_i] %*% solve(Sigma[minus_i, minus_i])
  arma::mat Sigma_i_minus = Sigma.rows(indiv_indices);
  Sigma_i_minus.shed_cols(indiv_indices); // Sigma[i, minus_i]
  arma::mat Sigma_minus_minus = Sigma.submat(minus_i, minus_i);
  // Resolve the linear system directly instead of computing the explicit
  // inverse: more stable numerically and closer to R's solve().
  arma::mat Si_t = arma::solve(Sigma_minus_minus, Sigma_i_minus.t());
  arma::mat Si = Si_t.t();

  // mean_i = Si %*% P[minus_i, ]
  arma::mat mean_i = Si * P.rows(minus_i);
  return mean_i;
}
