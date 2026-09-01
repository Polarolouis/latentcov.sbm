#ifndef P_CONDITIONALS_H
#define P_CONDITIONALS_H

#include "mvnorm.h"
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

// Mean of Pi given P_min_i and Sigma
arma::rowvec mean_of_Pi_given_P_min_i_sigma(arma::mat &P, arma::mat &Sigma,
                                            double sigma2, int i);

// Covariance of Pi given P_min_i and Sigma
arma::mat cov_of_Pi_given_P_min_i_sigma(arma::mat &P, arma::mat &Sigma,
                                        double sigma2, int i);
// Mean of Pi given P_min_i and Sigma
arma::mat block_mean_of_Pi_given_P_min_i_sigma(arma::mat &P, arma::mat &Sigma,
                                               double sigma2,
                                               arma::uvec indiv_indices);

#endif // P_CONDITIONALS_H