# Tool constants
TOL <- 1e-6


pivotCoordInv <- function(x, norm = "orthonormal") {
  if (!(norm %in% c("orthogonal", "orthonormal"))) stop("only orthogonal and orthonormal is allowd for norm")
  x <- -x
  y <- matrix(0, nrow = nrow(x), ncol = ncol(x) + 1)
  D <- ncol(x) + 1
  if (norm == "orthonormal") y[, 1] <- -sqrt((D - 1) / D) * x[, 1] else y[, 1] <- x[, 1]
  for (i in 2:ncol(y)) {
    for (j in 1:(i - 1)) {
      y[, i] <- y[, i] + x[, j] / if (norm == "orthonormal") sqrt((D - j + 1) * (D - j)) else 1
    }
  }
  for (i in 2:(ncol(y) - 1)) {
    y[, i] <- y[, i] - x[, i] * if (norm == "orthonormal") sqrt((D - i) / (D - i + 1)) else 1
  }
  yexp <- exp(y - apply(y, 1, max))
  x.back <- yexp / apply(yexp, 1, sum) # * rowSums(derOriginaldaten)
  if (is.data.frame(x)) x.back <- data.frame(x.back)
  return(x.back)
  # return(yexp)
}

simulate_P_and_Z <- function(K, Sigma, sigma2, rep_Z) {
  P <- t(mvtnorm::rmvnorm(n = K - 1, mean = rep(0, nrow(Sigma)), sigma = sigma2 * Sigma))

  probs <- pivotCoordInv(P)

  Z <- sapply(seq(rep_Z), function(Zidx) {
    sapply(seq_len(nrow(probs)), function(i) seq_len(K)[(rmultinom(n = 1, size = 1, prob = probs[i, ]) == 1)])
  }) |> t()
  return(list(Z = Z, P = P, probs = probs))
}

cond_Pi_fast <- function(P, Theta, sigma2, i) {
  n <- nrow(P)
  K1 <- ncol(P)

  idx <- setdiff(1:n, i)

  Theta_ii <- Theta[i, i]
  Theta_i_rest <- Theta[i, idx]
  P_rest <- P[idx, , drop = FALSE]

  mean_i <- -(1 / Theta_ii) * (Theta_i_rest %*% P_rest)
  cov_i <- (sigma2 / Theta_ii) * diag(K1)

  list(mean = as.vector(mean_i), cov = cov_i)
}

cond_Pi_given_P_min_i_sigma <- function(P, Sigma, sigma2, i) {
  n <- nrow(P)
  K_minus_1 <- ncol(P)
  minus_i <- setdiff(1:n, i)

  Si <- Sigma[i, minus_i, drop = FALSE] %*% solve(Sigma[minus_i, minus_i, drop = FALSE])

  mean_i <- Si %*% P[minus_i, , drop = FALSE]

  cov_i <- sigma2 * as.numeric(Sigma[i, i] - Si %*% Sigma[-i, i]) * diag(1, nrow = K_minus_1)

  return(list(mean = mean_i, cov = cov_i))
}

sample_Pi_given <- function(P, Sigma, sigma2, i) {
  res <- cond_Pi_given_P_min_i_sigma(P, Sigma, sigma2, i)

  # Draw from multivariate normal
  Pi_sample <- MASS::mvrnorm(
    n = 1,
    mu = res$mean,
    Sigma = res$cov
  )

  return(Pi_sample)
}

cat_dist_ilr_given_Pi <- function(Zi, Pi) {
  probs <- pivotCoordInv(matrix(Pi, nrow = 1))
  return(probs[Zi])
}

posterior_param_inv_gamma <- param_sigma2_given_P <- function(alpha_0, beta_0, P, Theta) {
  return(list(alpha = alpha_0 + (nrow(P) / 2), beta = beta_0 + 0.5 * sum(diag(t(P) %*% Theta %*% P))))
}

# Explicit inverse-gamma sampler using the identity:
# if X ~ IG(shape = a, rate = b) then 1 / X ~ Gamma(shape = a, rate = b)
sample_inv_gamma_rate <- sample_sigma2_given_P <- function(shape, rate) {
  return(1 / rgamma(n = 1, shape = shape, rate = rate))
}

## CLASSICAL LBM Z and pi


# pi | Z

param_pi_given_Z <- function(etas, Z) {
  etas + colSums(Z)
}

sample_pi_given_Z <- function(etas_post) {
  as.vector(MCMCpack::rdirichlet(n = 1, alpha = etas_post))
}

# Z | alpha, pi, Y, W
# Multinomial prob to normalize

param_multinom_probs_Z_poisson <- function(Y, alpha, W, pi, tol = TOL) {
  R_W <- Y %*% W
  N_W <- diag(colSums(W))
  unormalized_log_probs <- matrix(1, nrow = nrow(Y)) %*% log(pi) + R_W %*% log(t(alpha)) - matrix(1, nrow = nrow(Y), ncol = ncol(W)) %*% N_W %*% t(alpha)

  return(row_normalize_matrix(unormalized_log_probs, tol = tol))
}

sample_Z_given_alpha_pi_Y_W <- function(probs) {
  sapply(seq_len(nrow(probs)), function(i) {
    sample.int(n = ncol(probs), size = 1, replace = TRUE, prob = probs[i, ])
  })
}

## END OF CLASSICAL LBM

# rho | W

#' Compute the posterior parameters for the Dirichlet of row groups
#'
#' @param gammas a vector of size R (the number of column groups),
#' the prior of the Dirichlet
#' @param W a matrix of size \eqn{n_2 \times R} with a single 1 per line
#' indicating the membership of column node \eqn{j}, \eqn{W_{j,r} = 1} if \eqn{j} is in group \eqn{r} 0 else
#'
#' @return a vector of size R with the updated Dirichlet parameters
param_rho_given_W <- function(gammas, W) {
  gammas + colSums(W)
}

#' Sample from the posterior \eqn{\rho\mid W}
#'
#' A function to sample from the posterior distribution \eqn{\rho\mid W} which is a Dirichlet of posterior parameters
#'
#' @seealso [param_rho_given_W()] for the computations of the posterior parameters
sample_rho_given_W <- function(gammas_post) {
  as.vector(MCMCpack::rdirichlet(n = 1, alpha = gammas_post))
}

# Z | alpha, P, Y, W
# Multinomial prob to normalize

param_multinom_probs_Z_cov_poisson <- function(Y, alpha, W, P, tol = TOL) {
  R_W <- Y %*% W
  N_W <- diag(colSums(W))
  unormalized_log_probs <- log(pivotCoordInv(P)) + R_W %*% log(t(alpha)) - matrix(1, nrow = nrow(Y), ncol = ncol(W)) %*% N_W %*% t(alpha)

  return(row_normalize_matrix(unormalized_log_probs, tol = tol))
}

sample_Z_given_alpha_P_Y_W <- function(probs) {
  sapply(seq_len(nrow(probs)), function(i) {
    sample.int(n = ncol(probs), size = 1, replace = TRUE, prob = probs[i, ])
  })
}

# W | alpha, rho, Y, Z

param_multinom_probs_W_poisson <- function(Y, alpha, Z, rho, tol = TOL) {
  R_Z <- t(Y) %*% Z
  N_Z <- diag(colSums(Z))
  unormalized_log_probs <- matrix(1, nrow = ncol(Y)) %*% log(rho) + R_Z %*% log(alpha) - matrix(1, nrow = ncol(Y), ncol = ncol(Z)) %*% N_Z %*% alpha

  return(row_normalize_matrix(unormalized_log_probs, tol = tol))
}

sample_W_given_alpha_rho_Y_Z <- function(probs) {
  sapply(seq_len(nrow(probs)), function(i) {
    sample.int(n = ncol(probs), size = 1, replace = TRUE, prob = probs[i, ])
  })
}


# alpha | Y,Z,W

param_alpha_given_Y_Z_W_poisson <- function(a0, b0, Y, Z, W) {
  return(list(shape = a0 + t(Z) %*% Y %*% W, rate = b0 + t(Z) %*% matrix(1, nrow = nrow(Z), ncol = nrow(W)) %*% W))
}

sample_alpha_given_Y_Z_W_poisson <- function(shape, rate) {
  matrix(rgamma(length(shape), shape = shape, rate = rate), nrow = nrow(shape))
}

# Gibbs samplers

## Knowing Z labels
gibbs_sampling <- function(init_sigma2 = NULL, init_P = NULL, Sigma = NULL, Z = NULL, priors_hyper_params = list(alpha_0 = 3, beta_0 = 0.1), niter = 50L, niter_metropolis = 10L, minibatch = TRUE, sigma_fixed = NULL, K = NULL) {
  alpha_0 <- priors_hyper_params[["alpha_0"]]
  beta_0 <- priors_hyper_params[["beta_0"]]
  if (is.null(Sigma) || is.null(Z)) {
    stop("Sigma and Z must be provided")
  }
  if (is.null(K)) {
    K <- if (!is.null(init_P)) ncol(init_P) + 1L else max(Z)
  }
  if (is.null(init_sigma2)) {
    init_sigma2 <- sample_inv_gamma_rate(shape = alpha_0, rate = beta_0)
  }

  if (is.null(init_P)) {
    init_P <- t(mvtnorm::rmvnorm(n = K - 1, mean = rep(0, nrow(Sigma)), sigma = init_sigma2 * Sigma))
  }

  n <- nrow(Sigma)
  Theta <- solve(Sigma)

  post_sigma2 <- rep(NA, niter)
  post_P <- array(NA, dim = c(dim(init_P), niter))

  current_sigma2 <- ifelse(is.null(sigma_fixed), init_sigma2, sigma_fixed)
  current_P <- init_P

  post_sigma2[1] <- current_sigma2
  post_P[, , 1] <- current_P
  p <- progressr::progressor(steps = 2 * niter)
  metro_iters <- 0L
  accepted_metro <- 0L
  for (iter in seq(2, niter)) {
    p(sprintf("Updating sigma2"), amount = 0)
    # Updating sigma2
    if (is.null(sigma_fixed)) {
      current_sigma_params <- posterior_param_inv_gamma(alpha_0 = alpha_0, beta_0, current_P, Theta)
      current_sigma2 <- sample_inv_gamma_rate(shape = current_sigma_params$alpha, rate = current_sigma_params$beta)
    } else {
      current_sigma2 <- sigma_fixed
    }

    post_sigma2[iter] <- current_sigma2
    p(sprintf("Iter %d: current sigma2: %f", iter, current_sigma2), class = if (iter %% 10 == 0) "sticky")

    # Updating P
    ## Sampling an order of update for the Pis
    row_order <- if (minibatch) sample(x = n, size = n) else seq(1, n)
    for (ind_iter in seq(n)) {
      p(sprintf("Updating individual %i on %i", ind_iter, n), amount = 0)
      # Metropolis
      i <- row_order[ind_iter]
      for (iter_metro in seq(0, niter_metropolis - 1)) {
        Pi_candidate <- sample_Pi_given(P = current_P, Sigma = Sigma, sigma2 = current_sigma2, i = i)

        log_u <- log(runif(n = 1))

        log_accept <- sum(sapply(seq_len(nrow(Z)), function(Zidx) {
          log(cat_dist_ilr_given_Pi(Zi = Z[Zidx, i], Pi = Pi_candidate)) - log(cat_dist_ilr_given_Pi(Zi = Z[Zidx, i], Pi = current_P[i, ]))
        }))

        # log_accept <- log(cat_dist_ilr_given_Pi(Zi = Z[i], Pi = Pi_candidate)) - log(cat_dist_ilr_given_Pi(Zi = Z[i], Pi = current_P[i, ]))

        post_P[, , iter] <- current_P
        metro_iters <- metro_iters + 1L
        if (log_u < log_accept) {
          accepted_metro <- accepted_metro + 1L
          post_P[i, , iter] <- Pi_candidate
          current_P[i, ] <- Pi_candidate
        }
      }
    }
    p(sprintf("Iter %d: updated P", iter), class = if (iter %% 10 == 0) "sticky")
  }
  return(list(sigma2 = post_sigma2, P = post_P, P_accept = accepted_metro / metro_iters))
}

sample_P_metropolis_classical <- function(P, Z, Sigma, sigma2, minibatch = TRUE, niter_metropolis = 50L, rho = 1) {
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

      Pi_candidate <- P[i, ] + rnorm(n = ncol(P), sd = rho_iter)

      log_u <- log(runif(n = 1))

      # Conditional parameters
      P_candidate <- P
      P_candidate[i, ] <- Pi_candidate

      params_candidate <- cond_Pi_given_P_min_i_sigma(P_candidate, Sigma, sigma2, i)

      params_old <- cond_Pi_given_P_min_i_sigma(P, Sigma, sigma2, i)


      log_accept <- log(cat_dist_ilr_given_Pi(Zi = which.max(Z[i, ]), Pi = Pi_candidate)) - log(cat_dist_ilr_given_Pi(Zi = which.max(Z[i, ]), Pi = P[i, ])) + mvtnorm::dmvnorm(x = Pi_candidate, mean = params_candidate$mean, sigma = params_candidate$cov, log = TRUE) - mvtnorm::dmvnorm(x = P[i, ], mean = params_old$mean, sigma = params_old$cov, log = TRUE)


      # log_accept <- log(cat_dist_ilr_given_Pi(Zi = Z[i], Pi = Pi_candidate)) - log(cat_dist_ilr_given_Pi(Zi = Z[i], Pi = current_P[i, ]))

      if (log_u < log_accept) {
        accepted_count <- accepted_count + 1
        out_P[i, ] <- Pi_candidate
        P[i, ] <- Pi_candidate
      }
    }
  }
  message("Mean accepted rate ", accepted_count / (n * niter_metropolis))
  return(out_P)
}


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
## Classical LBM with Poisson

#' Gibbs sampler for Poisson LBM
#'
#' Runs a Gibbs sampler for a Poisson latent block model.
#'
#' @param Y Non-negative integer matrix of observed
#' counts.
#' @param init_Z Initial row-membership indicator
#' matrix (`nrow(Y)` x `K`).
#' @param init_W Initial column-membership indicator
#' matrix (`ncol(Y)` x `R`).
#' @param K Number of row groups.
#' @param R Number of column groups.
#' @param niter Number of Gibbs iterations.
#' @param priors_hyper_params List of hyperparameters:
#' `alpha_0`, `beta_0`,
#'   `gammas_0`, `a0`, `b0`.
#' @param tol Numeric; a small number indicating the
#' tolerance at which one wants to clamp the values,
#' i.e., if `value<tol` it gets set to tol and if
#' `value>1-tol` it gets set to `1-tol`. Only used for
#' probabilities.
#' @param verbose Logical; if `TRUE` the iterations
#' will print a message each 10 iterations
#' @param prefix Optional character string; a prefix to
#' append before logs. Will not print anything if
#' `verbose=FALSE`. Defaults to '' (empty string)
#'
#' @return A list with sampled trajectories:
#'   `sigma2_array`, `P_array`, `W_array`, `Z_array`, `rho_array`, `alpha_array`.
gibbs_sampling_lbm_poisson <- function(
  Y, init_Z, init_W, K, R,
  niter = 50L,
  priors_hyper_params = list(
    etas_0 = rep(2, K),
    gammas_0 = rep(2, R),
    a0 = 1, b0 = 1
  ),
  rho = 1,
  known_alpha = NULL,
  known_W = NULL,
  known_Z = NULL,
  known_pi = NULL,
  tol = TOL,
  verbose = FALSE,
  prefix = "",
  force_order = FALSE
) {
  # Forcing future exports
  invisible(c(TOL))
  # Initialize the whole arrays of variables
  rho_array <- array(NA, dim = c(niter, R), dimnames = list("Iteration" = seq(niter), "Parameter" = paste0("rho.", seq(1, R))))

  pi_array <- array(NA, dim = c(niter, K), dimnames = list("Iteration" = seq(niter), "Parameter" = paste0("pi.", seq(1, K))))

  alpha_array <- array(NA, dim = c(niter, K, R), dimnames = list("Iteration" = seq(niter), "RowGroup" = paste0("RowGroup", seq(1, K)), "ColGroup" = paste0("ColGroup", seq(1, R))))

  Z_array <- array(NA, dim = c(niter, nrow(Y)), dimnames = list("Iteration" = seq(niter), "Parameter" = paste0("Z.", seq_len(nrow(Y)))))

  W_array <- array(NA, dim = c(niter, ncol(Y)), dimnames = list("Iteration" = seq(niter), "Parameter" = paste0("W.", seq_len(ncol(Y)))))

  # Initialization
  ## Hyperparameters
  ### pi
  etas_0 <- priors_hyper_params[["etas_0"]]

  ### rho
  gammas_0 <- priors_hyper_params[["gammas_0"]]

  ### alpha
  a0 <- priors_hyper_params[["a0"]]
  b0 <- priors_hyper_params[["b0"]]

  ### Passing Z and W
  if (is.null(init_W)) {
    current_rho <- as.vector(MCMCpack::rdirichlet(n = 1, alpha = gammas_0))
    if (force_order) {
      col_order <- order(current_rho, decreasing = TRUE)
      current_rho <- current_rho[col_order]
    }
    W <- t(sapply(sample.int(n = R, size = ncol(Y), prob = current_rho, replace = TRUE), function(W_label) as.integer(seq(R) == W_label)))
  } else {
    W <- init_W
  }

  ### Z
  if (is.null(init_Z)) {
    current_pi <- as.vector(MCMCpack::rdirichlet(n = 1, alpha = etas_0))
    if (force_order) {
      row_order <- order(current_pi, decreasing = TRUE)
      current_pi <- current_pi[row_order]
    }
    Z <- t(sapply(sample.int(n = K, size = nrow(Y), prob = current_pi, replace = TRUE), function(Z_label) as.integer(seq(K) == Z_label)))
  } else {
    Z <- init_Z
  }

  for (iter in seq(niter)) {
    if (iter %% 10 == 0) {
      message(prefix, "Iter : ", iter, " on ", niter)
    }
    ### rho | W
    current_gammas <- param_rho_given_W(gammas = gammas_0, W)
    current_rho <- sample_rho_given_W(gammas_post = current_gammas)
    if (force_order) {
      col_order <- order(current_rho, decreasing = TRUE)
      current_rho <- current_rho[col_order]
    }
    rho_array[iter, ] <- current_rho

    ### pi | Z
    if (is.null(known_pi)) {
      current_etas <- param_pi_given_Z(etas = etas_0, Z)
      current_pi <- sample_pi_given_Z(etas_post = current_etas)
    } else {
      current_pi <- known_pi
    }

    if (force_order) {
      row_order <- order(current_pi, decreasing = TRUE)
      current_pi <- current_pi[row_order]
    }

    pi_array[iter, ] <- current_pi

    ### alpha | Y, Z, W
    if (is.null(known_alpha)) {
      alpha_params <- param_alpha_given_Y_Z_W_poisson(a0 = a0, b0 = b0, Y = Y, Z = Z, W = W)
      current_alpha <- sample_alpha_given_Y_Z_W_poisson(shape = alpha_params[["shape"]], rate = alpha_params[["rate"]])
    } else {
      current_alpha <- known_alpha
    }

    alpha_array[iter, , ] <- current_alpha


    ### W | Z,Y,rho,alpha
    if (is.null(known_W)) {
      W_post_probs <- param_multinom_probs_W_poisson(Y, current_alpha, Z, rho = current_rho, tol = tol)
      current_W_memb <- sample_W_given_alpha_rho_Y_Z(probs = W_post_probs)
    } else {
      current_W_memb <- known_W
    }
    W <- t(sapply(current_W_memb, function(W_label) {
      (seq(R) == W_label) * 1
    }))
    W_array[iter, ] <- current_W_memb

    ### Z | pi,W,Y,alpha
    if (is.null(known_Z)) {
      Z_post_probs <- param_multinom_probs_Z_poisson(Y = Y, alpha = current_alpha, W = W, pi = current_pi, tol = tol)
      current_Z_memb <- sample_Z_given_alpha_P_Y_W(probs = Z_post_probs)
    } else {
      current_Z_memb <- known_Z
    }
    Z <- t(sapply(current_Z_memb, function(Z_label) {
      (seq(K) == Z_label) * 1
    }))
    Z_array[iter, ] <- current_Z_memb
  }
  return(list(W_array = W_array, Z_array = Z_array, rho_array = rho_array, pi_array = pi_array, alpha_array = alpha_array))
}

#' @inheritParams gibbs_sampling_lbm_poisson
chains_gibbs_sampling_lbm_poisson <- function(nchains, ...) {
  lapply(seq(nchains), function(i) {
    gibbs_sampling_lbm_poisson(..., prefix = paste0("Chain ", i, " - "))
  }) |> futurize::futurize(seed = TRUE)
}

## Full LBM with latent phylo Poisson

#' Gibbs sampler for latent phylogenetic Poisson LBM
#'
#' Runs a Gibbs sampler for a Poisson latent block model with latent row effects
#' `P` structured by the covariance matrix `Sigma`.
#'
#' @param Sigma Numeric covariance matrix for row
#' latent effects.
#' @param Y Non-negative integer matrix of observed
#' counts.
#' @param init_Z Initial row-membership indicator
#' matrix (`nrow(Y)` x `K`).
#' @param init_W Initial column-membership indicator
#' matrix (`ncol(Y)` x `R`).
#' @param K Number of row groups.
#' @param R Number of column groups.
#' @param niter Number of Gibbs iterations.
#' @param niter_metropolis Number of Metropolis updates
#' per row for `P`.
#' @param priors_hyper_params List of hyperparameters:
#' `alpha_0`, `beta_0`,
#'   `gammas_0`, `a0`, `b0`.
#' @param sigma2_fixed Logical; if `TRUE`, keeps
#' `sigma2 = 1`, otherwise samples
#'   `sigma2` from its inverse-gamma full conditional.
#' @param minibatch Logical; if `TRUE` the update order
#' for P will change on each Gibbs iteration
#' @param tol Numeric; a small number indicating the
#' tolerance at which one wants to clamp the values,
#' i.e., if `value<tol` it gets set to tol and if
#' `value>1-tol` it gets set to `1-tol`. Only used for
#' probabilities.
#' @param verbose Logical; if `TRUE` the iterations
#' will print a message each 10 iterations
#' @param prefix Optional character string; a prefix to
#' append before logs. Will not print anything if
#' `verbose=FALSE`. Defaults to '' (empty string)
#'
#' @return A list with sampled trajectories:
#'   `sigma2_array`, `P_array`, `W_array`, `Z_array`, `rho_array`, `alpha_array`.
gibbs_sampling_lbm_cov_poisson <- function(
  Sigma, Y, init_Z, init_W, K, R,
  niter = 50L, niter_metropolis = 10L,
  priors_hyper_params = list(alpha_0 = 1, beta_0 = 1, gammas_0 = rep(2, R), a0 = 1, b0 = 1),
  rho = 1,
  sigma2_fixed = TRUE,
  known_alpha = NULL,
  known_P = NULL,
  known_W = NULL,
  known_Z = NULL,
  P_sampler = sample_P_metropolis_trick,
  minibatch = TRUE,
  tol = TOL,
  verbose = FALSE,
  prefix = ""
) {
  # Forcing future exports
  invisible(c(pivotCoordInv, cat_dist_ilr_given_Pi, sample_Pi_given, TOL))
  # Initialize the whole arrays of variables
  sigma2_array <- array(NA, dim = c(niter, 1), dimnames = list("Iteration" = seq(niter), "Parameter" = "sigma2"))

  P_array <- array(NA, dim = c(niter, nrow(Y), K - 1), dimnames = list("Iteration" = seq(niter), "Individual" = paste0("P", seq_len(nrow(Y))), "Coordinates" = seq(1, K - 1)))

  rho_array <- array(NA, dim = c(niter, R), dimnames = list("Iteration" = seq(niter), "Parameter" = paste0("rho.", seq(1, R))))

  alpha_array <- array(NA, dim = c(niter, K, R), dimnames = list("Iteration" = seq(niter), "RowGroup" = paste0("RowGroup", seq(1, K)), "ColGroup" = paste0("ColGroup", seq(1, R))))

  Z_array <- array(NA, dim = c(niter, nrow(Y)), dimnames = list("Iteration" = seq(niter), "Parameter" = paste0("Z.", seq_len(nrow(Y)))))

  W_array <- array(NA, dim = c(niter, ncol(Y)), dimnames = list("Iteration" = seq(niter), "Parameter" = paste0("W.", seq_len(ncol(Y)))))

  # Initialization

  Theta <- solve(Sigma)

  ## Hyperparameters
  ### sigma2
  alpha_0 <- priors_hyper_params[["alpha_0"]]
  beta_0 <- priors_hyper_params[["beta_0"]]

  ### rho
  gammas_0 <- priors_hyper_params[["gammas_0"]]

  ### alpha
  a0 <- priors_hyper_params[["a0"]]
  b0 <- priors_hyper_params[["b0"]]

  ### Passing Z and W
  if (is.null(init_W)) {
    current_rho <- as.vector(MCMCpack::rdirichlet(n = 1, alpha = gammas_0))
    W <- sapply(seq_len(ncol(Y)), function(j) {
      (seq(R) == sample.int(n = R, size = 1, replace = TRUE, prob = current_rho)) * 1
    }) |> t()
  } else {
    W <- init_W
  }


  ### sigma2
  if (!sigma2_fixed) {
    current_sigma2 <- sample_inv_gamma_rate(shape = alpha_0, rate = beta_0)
  } else {
    current_sigma2 <- 1L
  }

  ### P
  current_P <- t(mvtnorm::rmvnorm(n = K - 1, mean = rep(0, nrow(Sigma)), sigma = current_sigma2 * Sigma))
  dimnames(current_P) <- list("Individual" = paste0("P", seq_len(nrow(Y))), "Coordinates" = seq(1, K - 1))

  ### Z
  if (is.null(init_Z)) {
    Z <- sapply(seq_len(nrow(Y)), function(j) {
      (seq(K) == sample.int(n = K, size = 1, replace = TRUE, prob = pivotCoordInv(current_P)[j, ])) * 1
    }) |> t()
  } else {
    Z <- init_Z
  }

  for (iter in seq(niter)) {
    if (iter %% 10 == 0) {
      message(prefix, "Iter : ", iter, " on ", niter)
    }
    ### rho | W
    current_gammas <- param_rho_given_W(gammas = gammas_0, W)
    current_rho <- sample_rho_given_W(gammas_post = current_gammas)
    rho_array[iter, ] <- current_rho

    ### alpha | Y, Z, W
    if (is.null(known_alpha)) {
      alpha_params <- param_alpha_given_Y_Z_W_poisson(a0 = a0, b0 = b0, Y = Y, Z = Z, W = W)
      current_alpha <- sample_alpha_given_Y_Z_W_poisson(shape = alpha_params[["shape"]], rate = alpha_params[["rate"]])
    } else {
      current_alpha <- known_alpha
    }

    alpha_array[iter, , ] <- current_alpha

    ### P | sigma2, Z
    if (is.null(known_P)) {
      current_P <- P_sampler(P = current_P, Z = Z, Sigma = Sigma, sigma2 = current_sigma2, minibatch = minibatch, niter_metropolis = niter_metropolis, rho = rho)
    } else {
      current_P <- known_P
    }
    P_array[iter, , ] <- current_P

    ### sigma2 | P (not running currently)

    if (!sigma2_fixed) {
      sigma2_post_params <- param_sigma2_given_P(alpha_0 = alpha_0, beta_0, P = current_P, Theta = Theta)
      current_sigma2 <- sample_sigma2_given_P(shape = sigma2_post_params[["alpha"]], rate = sigma2_post_params[["beta"]])
    } else {
      current_sigma2 <- 1
    }
    sigma2_array[iter, ] <- current_sigma2

    ### W | Z,Y,rho,alpha
    if (is.null(known_W)) {
      W_post_probs <- param_multinom_probs_W_poisson(Y, current_alpha, Z, rho = current_rho, tol = tol)
      current_W_memb <- sample_W_given_alpha_rho_Y_Z(probs = W_post_probs)
    } else {
      current_W_memb <- known_W
    }
    W <- t(sapply(current_W_memb, function(W_label) {
      (seq(R) == W_label) * 1
    }))
    W_array[iter, ] <- current_W_memb

    ### Z | P,W,Y,alpha
    if (is.null(known_Z)) {
      Z_post_probs <- param_multinom_probs_Z_cov_poisson(Y = Y, alpha = current_alpha, W = W, P = current_P, tol = tol)
      current_Z_memb <- sample_Z_given_alpha_P_Y_W(probs = Z_post_probs)
    } else {
      current_Z_memb <- known_Z
    }
    Z <- t(sapply(current_Z_memb, function(Z_label) {
      (seq(K) == Z_label) * 1
    }))
    Z_array[iter, ] <- current_Z_memb
  }
  return(list(sigma2_array = sigma2_array, P_array = P_array, W_array = W_array, Z_array = Z_array, rho_array = rho_array, alpha_array = alpha_array))
}

#' @inheritParams gibbs_sampling_lbm_cov_poisson
chains_gibbs_sampling_lbm_cov_poisson <- function(nchains, ...) {
  lapply(seq(nchains), function(i) {
    gibbs_sampling_lbm_cov_poisson(..., prefix = paste0("Chain ", i, " - "))
  }) |> futurize::futurize(seed = TRUE)
}

mse <- function(x, y) sum((x - y)^2)
