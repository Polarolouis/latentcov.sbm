# Sufficient statistics of a Poisson LBM iteration
#
# The per-iteration aggregates consumed by the `alpha`, `W` and `Z` full
# conditionals. Computed once per iteration and reused across conditionals
# instead of being recomputed inside each one. A single source of truth for
# both hard (one-hot) and soft (probability) memberships; always mask-aware:
# missing cells are never counted as observed.

#' Row sufficient statistics of a Poisson LBM iteration
#'
#' Computes the row-side aggregates \eqn{R_W = Y_{obs} W} and
#' \eqn{N_W = \text{mask} W} used by the row-membership (`Z`) and connectivity
#' (`alpha`) full conditionals.
#'
#' @param Y the observed count matrix (\eqn{n_1 \times n_2}), with or without NAs
#' @param mask a matrix of size \eqn{n_1 \times n_2} with 1 for observed values
#'   and 0 for NAs
#' @param W the column memberships as an (\eqn{n_2 \times R}) one-hot matrix
#'   (hard) or a probability matrix (soft)
#'
#' @return a list with `R_W` (\eqn{n_1 \times R}) and `N_W` (\eqn{n_1 \times R})
#' @keywords internal
sufficient_stats_rows <- function(Y, mask, W) {
  Y_obs <- mask_observed(Y, mask)
  list(
    R_W = Y_obs %*% W,
    N_W = mask %*% W
  )
}

#' Column sufficient statistics of a Poisson LBM iteration
#'
#' Computes the column-side aggregates \eqn{R_Z = Y_{obs}^\top Z} and
#' \eqn{N_Z = \text{mask}^\top Z} used by the column-membership (`W`) and
#' connectivity (`alpha`) full conditionals.
#'
#' @param Y the observed count matrix (\eqn{n_1 \times n_2}), with or without NAs
#' @param mask a matrix of size \eqn{n_1 \times n_2} with 1 for observed values
#'   and 0 for NAs
#' @param Z the row memberships as an (\eqn{n_1 \times K}) one-hot matrix
#'   (hard) or a probability matrix (soft)
#'
#' @return a list with `R_Z` (\eqn{n_2 \times K}) and `N_Z` (\eqn{n_2 \times K})
#' @keywords internal
sufficient_stats_cols <- function(Y, mask, Z) {
  Y_obs <- mask_observed(Y, mask)
  list(
    R_Z = t(Y_obs) %*% Z,
    N_Z = t(mask) %*% Z
  )
}

#' Combined sufficient statistics of a Poisson LBM iteration
#'
#' Convenience wrapper computing both the row- and column-side sufficient
#' statistics in one call.
#'
#' @inheritParams sufficient_stats_rows
#' @inheritParams sufficient_stats_cols
#'
#' @return a list with `R_W`, `N_W`, `R_Z`, `N_Z`
#' @keywords internal
sufficient_stats <- function(Y, mask, Z, W) {
  cols <- sufficient_stats_cols(Y, mask, Z)
  rows <- sufficient_stats_rows(Y, mask, W)
  c(cols, rows)
}

#' Posterior Gamma parameters of the connectivity matrix `alpha`
#'
#' Derives the shape and rate of the full conditional
#' \eqn{\alpha \mid Y, Z, W} from the column sufficient statistics, avoiding a
#' third independent aggregation pass over `Y`.
#'
#' @param a0 prior shape parameter of the Gamma
#' @param b0 prior rate parameter of the Gamma
#' @param W the column memberships as an (\eqn{n_2 \times R}) one-hot matrix
#' @param stats_cols a list returned by [sufficient_stats_cols()]
#'
#' @return a list with `shape` and `rate` matrices (\eqn{K \times R})
#' @keywords internal
alpha_params_from_stats_cols <- function(a0, b0, W, stats_cols) {
  list(
    shape = a0 + t(stats_cols$R_Z) %*% W,
    rate = b0 + t(stats_cols$N_Z) %*% W
  )
}

#' Multinomial probabilities of the row memberships
#'
#' Computes the normalized conditional probabilities of each row-node block
#' membership from the row sufficient statistics and a per-row prior over the
#' blocks.
#'
#' @param stats_rows a list returned by [sufficient_stats_rows()]
#' @param alpha a matrix (\eqn{K \times R}) of connectivity coefficients
#' @param row_prior the row block proportions: either a vector of size K
#'   (`pi`, broadcast over rows) or an (\eqn{n_1 \times K}) matrix of per-row
#'   probabilities (`ilrInv(P)`)
#' @param tol Numeric; a small number used to clamp probabilities to
#'   `[tol, 1-tol]`. Defaults to `TOL`.
#'
#' @return a matrix (\eqn{n_1 \times K}) of normalized membership probabilities
#' @keywords internal
multinom_probs_Z <- function(stats_rows, alpha, row_prior, tol = TOL) {
  if (is.matrix(row_prior)) {
    prior_term <- log(row_prior)
  } else {
    prior_term <- matrix(1, nrow = nrow(stats_rows$R_W)) %*% log(row_prior)
  }
  unormalized_log_probs <- prior_term +
    stats_rows$R_W %*% log(t(alpha)) -
    stats_rows$N_W %*% t(alpha)

  row_normalize_matrix(unormalized_log_probs, tol = tol)
}

#' Multinomial probabilities of the column memberships
#'
#' Computes the normalized conditional probabilities of each column-node block
#' membership from the column sufficient statistics and the column block
#' proportions.
#'
#' @param stats_cols a list returned by [sufficient_stats_cols()]
#' @param alpha a matrix (\eqn{K \times R}) of connectivity coefficients
#' @param rho a vector of size R containing the column block proportions
#' @param tol Numeric; a small number used to clamp probabilities to
#'   `[tol, 1-tol]`. Defaults to `TOL`.
#'
#' @return a matrix (\eqn{n_2 \times R}) of normalized membership probabilities
#' @keywords internal
multinom_probs_W <- function(stats_cols, alpha, rho, tol = TOL) {
  unormalized_log_probs <- matrix(1, nrow = nrow(stats_cols$R_Z)) %*%
    log(rho) +
    stats_cols$R_Z %*% log(alpha) -
    stats_cols$N_Z %*% alpha

  row_normalize_matrix(unormalized_log_probs, tol = tol)
}
