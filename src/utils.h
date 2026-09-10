#ifndef UTILS_H
#define UTILS_H

#include "matrixnorm.h"
#include <RcppArmadillo.h>

arma::mat rmatrixnormal_cpp(const arma::mat &M, const arma::mat &U,
                            const arma::mat &V);

double dmatrixnormal_cpp(const arma::mat &X, const arma::mat &M,
                         const arma::mat &U, const arma::mat &V, bool log_p);

arma::mat default_Psi_function_cpp(int K);

arma::mat ilrInv_cpp(const arma::mat &z, const arma::mat &basis, bool log);

#endif // UTILS_H