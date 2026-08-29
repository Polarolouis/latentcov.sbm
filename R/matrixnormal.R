#' Compute density for X of a Matrix Normal with parameters M, U and V
#'
#' @param X a matrix of size n,p for which to compute the density
#' @param M a matrix of size n,p the mean of the matrix normal
#' @param U a variance-covariance matrix for the rows of Y
#' @param V a variance-covariance matrix for the columns of Y
#' @param log_p a boolean indicating if the log density should be returned
#'
#' @returns the density a scalar
dmatrixnormal_plain <- function(X, M, U, V, log_p = FALSE) {
  n <- nrow(X)
  p <- ncol(X)

  dev <- X - M
  U_inv <- solve(U)
  V_inv <- solve(V)
  det_U <- det(U)
  det_V <- det(V)

  exponent <- -0.5 * sum(diag(V_inv %*% t(dev) %*% U_inv %*% dev))

  log_norm_const <- -0.5 * (n * p * log(2 * pi) + n * log(det_V) + p * log(det_U))

  if (log_p) {
    out <- log_norm_const + exponent
  } else {
    out <- exp(log_norm_const + exponent)
  }
  return(out)
}

#' Compute density for X of a Matrix Normal with parameters M, U and V
#'
#' @param X a matrix of size n,p for which to compute the density
#' @param M a matrix of size n,p the mean of the matrix normal
#' @param U a variance-covariance matrix for the rows of Y
#' @param V a variance-covariance matrix for the columns of Y
#' @param log_p a boolean indicating if the log density should be returned
#'
#' @returns the density a scalar
dmatrixnormal_schur <- function(X, M, U, V, log_p = FALSE) {
  n <- nrow(X)
  p <- ncol(X)

  dev <- X - M
  U_inv <- solve(U)
  V_inv <- solve(V)
  det_U <- det(U)
  det_V <- det(V)

  exponent <- -0.5 * sum((U_inv %*% dev) * (dev %*% V_inv))

  log_norm_const <- -0.5 * (n * p * log(2 * pi) + n * log(det_V) + p * log(det_U))

  if (log_p) {
    out <- log_norm_const + exponent
  } else {
    out <- exp(log_norm_const + exponent)
  }
  return(out)
}


rmatrixnormal <- function(M, U, V) {
  n <- nrow(M)
  p <- ncol(M)

  X <- matrix(rnorm(n = n * p), n, p)
  A <- t(chol(U))
  B <- chol(V)

  Y <- M + A %*% X %*% B
  return(Y)
}
