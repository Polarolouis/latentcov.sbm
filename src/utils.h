#ifndef UTILS_H
#define UTILS_H

#include "matrixnorm.h"
#include <RcppArmadillo.h>

arma::mat rmatrixnormal_cpp(const arma::mat &M, const arma::mat &U,
                            const arma::mat &V);

double dmatrixnormal_cpp(const arma::mat &X, const arma::mat &M,
                         const arma::mat &U, const arma::mat &V, bool log_p);

arma::mat pivot_coord_inv(arma::mat &x, std::string norm, bool log);

#endif // UTILS_H