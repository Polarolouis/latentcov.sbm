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
#' @param emission a string indicating the emission distribution for
#'  which to compute sufficient statistics. Defaults to "poisson"
#'
#' @return a list with `R_W` (\eqn{n_1 \times R}) and `N_W` (\eqn{n_1 \times R})
#' @keywords internal
sufficient_stats_rows <- function(Y, mask, W, emission = "poisson") {
  Y_obs <- mask_observed(Y, mask)
  check_emission_distribution(emission)
  switch(
    emission,
    "poisson" = {
      out <- list(
        R_W = Y_obs %*% W,
        N_W = mask %*% W
      )
    },
    "bernoulli" = {
      out <- list(R_W = Y_obs %*% W, F_W = (1 - Y_obs) %*% W)
    }
  )
  return(out)
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
#' @param emission a string indicating the emission distribution for
#'  which to compute sufficient statistics. Defaults to "poisson"
#'
#'
#' @return a list with `R_Z` (\eqn{n_2 \times K}) and `N_Z` (\eqn{n_2 \times K})
#' @keywords internal
sufficient_stats_cols <- function(Y, mask, Z, emission = "poisson") {
  Y_obs <- mask_observed(Y, mask)
  check_emission_distribution(emission)
  switch(
    emission,
    "poisson" = {
      out <- list(
        R_Z = t(Y_obs) %*% Z,
        N_Z = t(mask) %*% Z
      )
    },
    "bernoulli" = {
      out <- list(
        R_Z = t(Y_obs) %*% Z,
        F_Z = (1 - t(Y_obs)) %*% Z
      )
    }
  )
  return(out)
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
sufficient_stats <- function(Y, mask, Z, W, emission = "poisson") {
  cols <- sufficient_stats_cols(Y, mask, Z, emission = emission)
  rows <- sufficient_stats_rows(Y, mask, W, emission = emission)
  c(cols, rows)
}

#' Posterior parameters of the connectivity matrix `alpha`
#'
#' Derives the shape and rate of the full conditional
#' \eqn{\alpha \mid Y, Z, W} from the column sufficient statistics
#'
#' @param ... the prior parameters for the alpha distribution
#' @param W the column memberships as an (\eqn{n_2 \times R}) one-hot matrix
#' @param stats_cols a list returned by [sufficient_stats_cols()]
#'
#' @details On the different emission distributions:
#' \itemize{
#'   \item{`bernoulli:`}{takes `a0` and `b0` as priors for the Beta distribution
#' and outputs `shape1` and `shape2` the parameters for the
#' a posteriori Beta distribution.}
#'   \item{`poisson:`}{takes `a0` and `b0` as priors for the Gamma distribution
#' and outputs `shape` and `rate` the parameters for the
#' a posteriori Gamma distribution.}
#' }
#'
#' @return a list containing the updated parameters for the `emission`
#' distribution provided. See details for the exact names
#' @keywords internal
alpha_params_from_stats_cols <- function(
  ...,
  W,
  stats_cols,
  emission = "poisson"
) {
  rlang::check_required(W)
  rlang::check_required(stats_cols)
  check_emission_distribution(emission)
  args <- rlang::list2(...)
  # TODO Choose a priori parameters names that are more coherent with emission distribution
  switch(
    emission,
    "poisson" = {
      a0 <- args[["a0"]]
      b0 <- args[["b0"]]
      out <- list(
        shape = a0 + t(stats_cols$R_Z) %*% W,
        rate = b0 + t(stats_cols$N_Z) %*% W
      )
    },
    "bernoulli" = {
      a0 <- args[["a0"]]
      b0 <- args[["b0"]]
      out <- list(
        shape1 = a0 + t(stats_cols$R_Z) %*% W,
        shape2 = b0 + t(stats_cols$F_Z) %*% W
      )
    }
  )
  return(out)
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
alpha_params_from_stats_cols_poisson <- function(a0, b0, W, stats_cols) {
  alpha_params_from_stats_cols(
    a0 = a0,
    b0 = b0,
    W = W,
    stats_cols = stats_cols,
    emission = "poisson"
  )
}

#' Posterior Beta parameters of the connectivity matrix `alpha`
#'
#' Derives the shape and rate of the full conditional
#' \eqn{\alpha \mid Y, Z, W} from the column sufficient statistics.
#'
#' @param a0 prior shape1 parameter of the Beta
#' @param b0 prior shape2 parameter of the Beta
#' @param W the column memberships as an (\eqn{n_2 \times R}) one-hot matrix
#' @param stats_cols a list returned by [sufficient_stats_cols()]
#'
#' @return a list with `shape1` and `shape2` matrices (\eqn{K \times R})
#' @keywords internal
alpha_params_from_stats_cols_bernoulli <- function(a0, b0, W, stats_cols) {
  alpha_params_from_stats_cols(
    a0 = a0,
    b0 = b0,
    W = W,
    stats_cols = stats_cols,
    emission = "bernoulli"
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
#' @param emission a string indicating the emission distribution for which to
#'   compute the multinomial probabilities. Defaults to "poisson"
#' @param tol Numeric; a small number used to clamp probabilities to
#'   `[tol, 1-tol]`. Defaults to `TOL`.
#'
#' @return a matrix (\eqn{n_1 \times K}) of normalized membership probabilities
#' @keywords internal
multinom_probs_Z <- function(
  stats_rows,
  alpha,
  row_prior,
  emission = "poisson",
  tol = TOL
) {
  check_emission_distribution(emission)

  if (is.matrix(row_prior)) {
    prior_term <- .log(row_prior)
  } else {
    prior_term <- matrix(1, nrow = nrow(stats_rows$R_W)) %*% .log(row_prior)
  }

  loglik_term <- switch(
    emission,
    "poisson" = stats_rows$R_W %*% .log(t(alpha)) - stats_rows$N_W %*% t(alpha),
    "bernoulli" = stats_rows$R_W %*%
      .log(t(alpha)) +
      stats_rows$F_W %*% .log(1 - t(alpha))
  )

  row_normalize_matrix(prior_term + loglik_term, tol = tol)
}

#' @rdname multinom_probs_Z
#' @keywords internal
multinom_probs_Z_poisson <- function(stats_rows, alpha, row_prior, tol = TOL) {
  multinom_probs_Z(
    stats_rows = stats_rows,
    alpha = alpha,
    row_prior = row_prior,
    emission = "poisson",
    tol = tol
  )
}

#' @rdname multinom_probs_Z
#' @keywords internal
multinom_probs_Z_bernoulli <- function(
  stats_rows,
  alpha,
  row_prior,
  tol = TOL
) {
  multinom_probs_Z(
    stats_rows = stats_rows,
    alpha = alpha,
    row_prior = row_prior,
    emission = "bernoulli",
    tol = tol
  )
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
#' @param emission a string indicating the emission distribution for which to
#'   compute the multinomial probabilities. Defaults to "poisson"
#' @param tol Numeric; a small number used to clamp probabilities to
#'   `[tol, 1-tol]`. Defaults to `TOL`.
#'
#' @return a matrix (\eqn{n_2 \times R}) of normalized membership probabilities
#' @keywords internal
multinom_probs_W <- function(
  stats_cols,
  alpha,
  rho,
  emission = "poisson",
  tol = TOL
) {
  check_emission_distribution(emission)

  prior_term <- matrix(1, nrow = nrow(stats_cols$R_Z)) %*% .log(rho)

  loglik_term <- switch(
    emission,
    "poisson" = stats_cols$R_Z %*% .log(alpha) - stats_cols$N_Z %*% alpha,
    "bernoulli" = stats_cols$R_Z %*%
      .log(alpha) +
      stats_cols$F_Z %*% .log(1 - alpha)
  )

  row_normalize_matrix(prior_term + loglik_term, tol = tol)
}

#' @rdname multinom_probs_W
#' @keywords internal
multinom_probs_W_poisson <- function(stats_cols, alpha, rho, tol = TOL) {
  multinom_probs_W(
    stats_cols = stats_cols,
    alpha = alpha,
    rho = rho,
    emission = "poisson",
    tol = tol
  )
}

#' @rdname multinom_probs_W
#' @keywords internal
multinom_probs_W_bernoulli <- function(stats_cols, alpha, rho, tol = TOL) {
  multinom_probs_W(
    stats_cols = stats_cols,
    alpha = alpha,
    rho = rho,
    emission = "bernoulli",
    tol = tol
  )
}
