#### Stick-breaking

stick_breaking <- function(x) {
  sapply(seq_len(length(x) - 1), function(idx) {
    if (idx == 1) {
      return(exp(x[idx]) / sum(exp(x)))
    }
    return(.logit(x[idx]) * sum(.logit(x[seq(idx - 1)])))
  })
}


#### Change-of-basis

#' Default Psi contrast matrix
#'
#' Builds the default \eqn{(K-1) \times K} contrast (basis) matrix used by the
#' isometric log-ratio transformation of the block
#' proportions.
#'
#' @param K Number of groups (blocks).
#'
#' @return A matrix of size \eqn{(K-1)\times K}.
default_Psi_function <- function(K) {
  simple_Psi_function(K)
}

#' Build the simplest change-of-basis matrix of size \eqn{(K-1)\times K}
#' by outputting a diagonal with one and a last column of zeros
#'
#' @inheritParams default_Psi_function
simple_Psi_function <- function(K) {
  diag(K)[-K, ]
}

#' Builds the change-of-basis from Egozcue et al. 2003,2005 or Hron et al. 2010
#' @inheritParams default_Psi_function
egozcue_Psi_function <- function(K) {
  t(sapply(seq(K - 1), function(j) {
    sqrt((K - j) / (K - j + 1)) * c(rep(0, j - 1), -1, rep(1 / (K - j), K - j))
  }))
}

##### Log ratio transforms
clr <- function(x) {
  lx <- log(x)
  lx_bar <- apply(lx, 1, mean)
  return(lx - lx_bar)
}

ilr <- function(x, basis = default_Psi_function(ncol(x))) {
  if (is.vector(x)) {
    cli::cli_warn("{.arg x} is a vector, casting it to a one row matrix.")
    x <- matrix(x, nrow = 1)
  }
  clr(x) %*% t(basis)
}

#' Inverse isometric log-ratio transformation
#'
#' Performs the inverse isometric log-ratio transformation for the latent covariance stochastic block model
#' @param z Input matrix
#' @param basis The basis to use
#' @return Transformed matrix
ilrInv <- function(z, basis = default_Psi_function(ncol(z) + 1)) {
  clr <- z %*% basis
  e <- exp(clr - apply(clr, 1, max))
  e / rowSums(e)
}

#### log ilrInv gradient, hessian

#' Log of the density of \eqn{p(Z\mid P)}
#'
#' @param P the latent position for \eqn{n_1} row nodes
#' @param Z the hard or soft memberships of the row nodes
#' @param basis the \eqn{(K-1)\times K} change-of-basis matrix used to link \eqn{\Delta_K} to \eqn{\mathbf{R}^{K-1}}
#' @param transformation the transformation to apply to the P to
#' obtain the probabilities
#'
#'
log_Z_given_P <- function(
  P,
  Z,
  basis = default_Psi_function(ncol(P) + 1),
  transformation = ilrInv
) {
  log_probs <- log(transformation(P, basis = basis))
  return(sum(Z * log_probs))
}

#' Compute \eqn{\nabla_P \operatorname{ilr}^{-1}(P)}
#'
#' @inheritParams log_Z_given_P
#' @return an \eqn{n_1 \times (K-1)} gradient for the \eqn{\operatorname{ilr}^{-1}(P)}
#' @keywords internal
log_Z_given_P_grad <- function(
  P,
  Z,
  basis = default_Psi_function(ncol(P) + 1),
  transformation = ilrInv
) {
  (Z - transformation(P, basis = basis)) %*% t(basis)
}

#' Compute \eqn{\nabla_P^2 \operatorname{ilr}^{-1}(P)}
#'
#' @inheritParams log_Z_given_P_grad
#'
#' @return an \eqn{n_1(K-1)\times n_1(K-1)} Hessian for the \eqn{\operatorname{ilr}^{-1}(P)}
#' @keywords internal
log_Z_given_P_hessian_row <- function(
  P,
  basis = default_Psi_function(ncol(P) + 1),
  transformation = ilrInv
) {
  probs <- transformation(P)
  n1 <- nrow(P)
  K <- ncol(P) + 1

  block_i <- function(prob_i) {
    basis %*% (diag(prob_i) - prob_i %*% t(prob_i)) %*% t(basis)
  }
  hessian <- matrix(0, nrow = n1 * (K - 1), ncol = n1 * (K - 1))
  for (i in seq(n1)) {
    b_coord <- seq((i - 1) * (K - 1) + 1, i * (K - 1))
    hessian[b_coord, b_coord] <- -block_i(probs[i, ])
  }
  return(hessian)
}


#### matrix normal prior gradient, hessian
#' Computes the log prior density for the matrix normal distribution
#' with mean 0, row covariance \eqn{\Sigma} and column covariance
#' \eqn{\sigma^2 I_{K-1}}
#'
#' @param Sigma the variance-covariance matrix between the rows of P
#' @param sigma2 the variance between columns
#' @inheritParams log_Z_given_P_grad
log_dmatrixnormal <- function(P, Sigma, sigma2) {
  dmatrixnormal_cpp(
    X = P,
    M = matrix(0, nrow(P), ncol(P)),
    U = Sigma,
    V = sigma2 * diag(ncol(P)),
    log_p = TRUE
  )
}

#' Compute \eqn{\nabla_P \mathcal{N}(0, \Sigma,\sigma^2 I_{K-1}(P))}
#'
#' @param Theta the precision matrix \eqn{\Sigma^{-1}}
#' @param sigma2 the variance between columns
#' @inheritParams log_dmatrixnormal
#'
#' @return an \eqn{n_1 \times (K-1)} gradient for the \eqn{\mathcal{N}(0, \Sigma,\sigma^2 I_{K-1}(P)}
#' @keywords internal
log_dmatrixnormal_grad <- function(
  P,
  Theta,
  sigma2
) {
  return(-1 / sigma2 * Theta %*% P)
}

#' Compute \eqn{\nabla_{\operatorname{vec}_r(P)}^2 \operatorname{ilr}^{-1}(P)}
#'
#' @inheritParams log_dmatrixnormal_grad
#'
#' @return an \eqn{n_1(K-1)\times n_1(K-1)} Hessian for the \eqn{\mathcal{N}(0, \Sigma,\sigma^2 I_{K-1}(P)}
#' @keywords internal
log_dmatrixnormal_hessian_vec_row <- function(
  Theta,
  sigma2,
  K
) {
  hessian <- -(1 / sigma2) * (kronecker(Theta, diag(K - 1)))

  return(hessian)
}

#' Compute \eqn{\nabla_{\operatorname{vec}(P)}^2 \operatorname{ilr}^{-1}(P)}
#'
#' @inheritParams log_dmatrixnormal_grad
#'
#' @return an \eqn{n_1(K-1)\times n_1(K-1)} Hessian for the \eqn{\mathcal{N}(0, \Sigma,\sigma^2 I_{K-1}(P)}
#' @keywords internal
log_dmatrixnormal_hessian_vec_column <- function(
  Theta,
  sigma2,
  K
) {
  hessian <- -(1 / sigma2) * (kronecker(diag(K - 1), Theta))

  return(hessian)
}

#### log prior + log ilrInv

#' Compute the log prior + log ilrInv
#'
#' In details \eqn{\log [p(P\mid \Sigma, \sigma^2)p(Z\mid P)
#' \propto  p(P\mid Z)]}, called \eqn{h(P,Z)} in the paper
#'
#' @inheritParams log_Z_given_P
#' @inheritParams log_dmatrixnormal
#'
#' @return the combined log probabilities
log_prior_plus_log_Z_given_P <- function(
  P,
  Z,
  Sigma,
  sigma2,
  basis = default_Psi_function(ncol(P) + 1),
  transformation = ilrInv
) {
  log_dmatrixnormal(P, Sigma, sigma2) +
    log_Z_given_P(P, Z, basis, transformation)
}

#' Compute the gradient of log prior + log p(Z|P)
#'
#' @inheritParams log_prior_plus_log_Z_given_P
#'
#' @return the combined log probabilities
log_prior_plus_log_Z_given_P_grad <- function(
  P,
  Z,
  Sigma,
  sigma2,
  basis = default_Psi_function(ncol(P) + 1),
  transformation = ilrInv
) {
  log_dmatrixnormal_grad(P, Sigma, sigma2) +
    log_Z_given_P_grad(P, Z, basis, transformation)
}

#' Compute the Hessian of log prior + log p(Z|P)
#'
#' @inheritParams log_prior_plus_log_Z_given_P
#' @inheritParams log_dmatrixnormal_hessian
#'
#' @return the combined log probabilities
log_prior_plus_log_Z_given_P_hessian_vec_row <- function(
  P,
  Sigma,
  sigma2,
  basis = default_Psi_function(ncol(P) + 1),
  transformation = ilrInv
) {
  log_dmatrixnormal_hessian_vec_row(solve(Sigma), sigma2, K = ncol(P) + 1) +
    log_Z_given_P_hessian_row(P, basis, transformation)
}
