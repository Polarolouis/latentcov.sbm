#' A Metropolis-Hasting sampler for with standard random walk
#'
#' @importFrom stats rnorm runif
#'
#' @param P a matrix of size \eqn{n_1 \times K-1} specifying latent
#' position of the nodes probabilities of membership in the latent
#' space
#' @param Z a matrix of size \eqn{n_1\times K} with a single 1 per line
#' indicating the membership of row node \eqn{i}, \eqn{Z_{i,k} = 1} if \eqn{i} is in group \eqn{k} 0 else
#' @param Sigma a covariance matrix describing the covariance between row nodes \eqn{(n_1\times n_1)}
#' @param sigma2 a variance parameter indicating the variance between the K-1 columns of P
#' @param minibatch a boolean indicating wether to update the rows in the lexical order or to sample at each iteration. Default to TRUE
#' @param niter_metropolis an integer specifying the number of metropolis iterations to perform. Defaults to 50.
#' @param rho a double indicating the variance of the gaussian proposal distribution.
#' @param verbose a boolean indicating if verbose message should be outputted.
#'
#' @export
sample_P_metropolis_classical <- function(P, Z, Sigma, sigma2, minibatch = TRUE, niter_metropolis = 50L, rho = 1, verbose = getOption("latentcov.sbm.verbose", default = FALSE)) {
  n <- nrow(P)
  out_P <- array(P, dim = dim(P), dimnames = list("Individual" = paste0("P", seq_len(nrow(P))), "Coordinates" = seq_len(ncol(P))))
  # Updating P
  ## Sampling an order of update for the Pis
  accepted_count <- 0
  row_order <- if (minibatch) sample(x = n, size = n) else seq(1, n)
  for (ind_iter in seq(n)) {
    # Metropolis
    i <- row_order[ind_iter]
    for (iter_metro in seq(niter_metropolis)) {
      rho_iter <- rho * sample(c(1, 1 / 10, 10), size = 1)
      noise <- rnorm(n = ncol(P), sd = sqrt(rho_iter))
      Pi_candidate <- P[i, ] + noise
      log_u <- log(runif(n = 1))

      # Conditional parameters
      P_candidate <- P
      P_candidate[i, ] <- Pi_candidate

      params <- cond_Pi_given_P_min_i_sigma(P_candidate, Sigma, sigma2, i)

      log_accept <- log(cat_dist_ilr_given_Pi(Zi = which.max(Z[i, ]), Pi = Pi_candidate)) - log(cat_dist_ilr_given_Pi(Zi = which.max(Z[i, ]), Pi = P[i, ])) + mvtnorm::dmvnorm(x = Pi_candidate, mean = params$mean, sigma = params$cov, log = TRUE) - mvtnorm::dmvnorm(x = P[i, ], mean = params$mean, sigma = params$cov, log = TRUE)

      if (log_u < log_accept) {
        accepted_count <- accepted_count + 1
        out_P[i, ] <- Pi_candidate
        P[i, ] <- Pi_candidate
      }
    }
  }
  if (verbose) {
    message("Mean accepted rate ", accepted_count / (n * niter_metropolis))
  }
  return(out_P)
}

#' A Metropolis-Hasting sampler for with a clever proposition
#'
#' @inheritParams sample_P_metropolis_classical
#'
#' @param ... Used to pass various args
#' @export
sample_P_metropolis_trick <- function(P, Z, Sigma, sigma2, minibatch = TRUE, niter_metropolis = 50L, ...) {
  n <- nrow(P)
  out_P <- array(P, dim = dim(P), dimnames = list("Individual" = paste0("P", seq_len(nrow(P))), "Coordinates" = seq_len(ncol(P))))
  # Updating P
  ## Sampling an order of update for the Pis
  row_order <- if (minibatch) sample(x = n, size = n) else seq(1, n)
  for (ind_iter in seq(n)) {
    # Metropolis
    i <- row_order[ind_iter]
    for (iter_metro in seq(niter_metropolis)) {
      Pi_candidate <- sample_Pi_given(P = P, Sigma = Sigma, sigma2 = sigma2, i = i)

      log_u <- log(runif(n = 1))

      log_accept <- log(cat_dist_ilr_given_Pi(Zi = which.max(Z[i, ]), Pi = Pi_candidate)) - log(cat_dist_ilr_given_Pi(Zi = which.max(Z[i, ]), Pi = P[i, ]))

      # log_accept <- log(cat_dist_ilr_given_Pi(Zi = Z[i], Pi = Pi_candidate)) - log(cat_dist_ilr_given_Pi(Zi = Z[i], Pi = current_P[i, ]))

      if (log_u < log_accept) {
        out_P[i, ] <- Pi_candidate
        P[i, ] <- Pi_candidate
      }
    }
  }
  return(out_P)
}

#' Soft Metropolis sampler for P under mixed membership
#'
#' Metropolis-Hastings sampler for the latent positions `P` with the clever
#' proposal (as in [sample_P_metropolis_trick()]) but under a
#' *mixed-membership* likelihood: instead of conditioning on the single
#' hard assignment `which.max(Z[i, ])`, the acceptance ratio weights the
#' categorical log-likelihood over all blocks by the soft memberships of row
#' `i`.
#'
#' @inheritParams sample_P_metropolis_trick
#' @param Z a matrix (\eqn{n_1 \times K}) of **soft** membership
#'   probabilities (rows sum to 1) in lieu of a hard indicator matrix
#'
#' @return an updated matrix of latent positions, of the same size as `P`
#' @export
sample_P_metropolis_trick_soft <- function(P, Z, Sigma, sigma2, minibatch = TRUE, niter_metropolis = 50L, ...) {
  n <- nrow(P)
  out_P <- array(P, dim = dim(P), dimnames = list("Individual" = paste0("P", seq_len(nrow(P))), "Coordinates" = seq_len(ncol(P))))
  row_order <- if (minibatch) sample(x = n, size = n) else seq(1, n)
  for (ind_iter in seq(n)) {
    i <- row_order[ind_iter]
    for (iter_metro in seq(niter_metropolis)) {
      Pi_candidate <- sample_Pi_given(P = P, Sigma = Sigma, sigma2 = sigma2, i = i)

      # Mixed-membership log-likelihood ratio:
      #   sum_k Z_soft[i,k] * (log ilrInv(P_i^cand)[k] - log ilrInv(P_i)[k])
      log_accept <- sum(Z[i, ] * (log(ilrInv(matrix(Pi_candidate, nrow = 1))) - log(ilrInv(matrix(P[i, ], nrow = 1)))))
      log_u <- log(runif(n = 1))

      if (log_u < log_accept) {
        out_P[i, ] <- Pi_candidate
        P[i, ] <- Pi_candidate
      }
    }
  }
  return(out_P)
}
