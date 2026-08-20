#include "P_samplers.h"

// [[Rcpp::depends(RcppArmadillo)]]

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

using namespace arma;
// [[Rcpp::export]]
arma::mat sample_P_classical(arma::mat &P, arma::mat &Z, arma::mat &Sigma,
                             double sigma2, bool minibatch = true,
                             int niter_metropolis = 50, double rho = 1.0) {
  mat running_P = P;
  int n = running_P.n_rows;

  vec rho_values = {1.0, 0.1, 10};

  uvec row_order;

  if (minibatch) {
    row_order = arma::randperm(n);
  } else {
    row_order = arma::regspace<uvec>(0, n - 1);
  }

  // Loop over the individuals
  uvec::iterator ind_end = row_order.end();
  for (uvec::iterator ind_id = row_order.begin(); ind_id != ind_end; ++ind_id) {
    // Sample the variance of the normal
    int indiv_idx = *ind_id;
    double rho_iter =
        rho * Rcpp::RcppArmadillo::sample(rho_values, 1, false)(0);
    ;

    vec Pi_candidate = running_P.row(indiv_idx);
    arma::rowvec noise =
        arma::randn<arma::rowvec>(P.n_cols) * std::sqrt(rho_iter);

    Pi_candidate += noise;

    double log_u = log(Rcpp::runif(1)(0));

    double log_accept = 0;
  }

  return running_P;
};

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