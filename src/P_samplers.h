#ifndef P_SAMPLERS_H
#define P_SAMPLERS_H

#include "mvnorm.h"
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>

// [[Rcpp::depends(RcppArmadillo)]]

// Pivot coordinate inversion
arma::mat pivot_coord_inv(arma::mat &x, std::string norm, bool log);

// Classical sampling of P
arma::mat sample_P_classical(arma::mat &P, arma::mat &Z, arma::mat &Sigma,
                             double sigma2, bool minibatch,
                             int niter_metropolis, double rho);

// Mean of Pi given P_min_i and Sigma
arma::rowvec mean_of_Pi_given_P_min_i_sigma(arma::mat &P, arma::mat &Sigma,
                                            double sigma2, int i);

// Covariance of Pi given P_min_i and Sigma
arma::mat cov_of_Pi_given_P_min_i_sigma(arma::mat &P, arma::mat &Sigma,
                                        double sigma2, int i);

#endif // P_SAMPLERS_H