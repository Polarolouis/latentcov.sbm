#ifndef P_SAMPLERS_H
#define P_SAMPLERS_H
#define DEBUG_SAMPLE false

#include "P_conditionals.h"
#include "mvnorm.h"
#include "utils.h"
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>

// [[Rcpp::depends(RcppArmadillo)]]

// Classical sampling of P
arma::mat sample_P_metropolis_classical_cpp(arma::mat &P, arma::mat &Z,
                                            arma::mat &Sigma, double sigma2,
                                            bool minibatch,
                                            int niter_metropolis, double rho);

// Sampling of P using smart kernel
arma::mat sample_P_metropolis_trick_cpp(arma::mat &P, arma::mat &Z,
                                        arma::mat &Sigma, double sigma2,
                                        bool minibatch, int niter_metropolis,
                                        double rho);

// Sampling of P using smart kernel with block updates
arma::mat sample_P_metropolis_trick_cpp_block(arma::mat &P, arma::mat &Z,
                                              arma::mat &Sigma, double sigma2,
                                              bool minibatch,
                                              int niter_metropolis,
                                              double rho, int block_size);

// Sampling of P using smart kernel and the precision matrix Theta (fast)
arma::mat sample_P_metropolis_trick_cpp_fast(arma::mat &P, arma::mat &Z,
                                             arma::mat &Theta, double sigma2,
                                             bool minibatch,
                                             int niter_metropolis);

// Sampling of P under soft memberships using the precision matrix Theta (fast)
arma::mat sample_P_metropolis_trick_soft_cpp_fast(arma::mat &P, arma::mat &Z,
                                                  arma::mat &Theta,
                                                  double sigma2, bool minibatch,
                                                  int niter_metropolis);

#endif // P_SAMPLERS_H