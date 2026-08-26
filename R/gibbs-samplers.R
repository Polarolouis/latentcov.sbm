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
#' @param known_alpha NULL (default) or a matrix (\eqn{K\times R}) indicating the known \eqn{\alpha} connectivity matrix
#' @param known_pi NULL (default) or a vector (\eqn{K})
#' indicating the known row block proportions
#' indicating the known columns block memberships
#' @param known_Z NULL (default) or a matrix (\eqn{n_1\times K})
#' indicating the known columns block memberships
#' @param known_W NULL (default) or a matrix (\eqn{n_2\times R})
#' indicating the known columns block memberships
#' @param force_order a boolean indicating if the order the groups
#' \eqn{\{1,\dots,K\}} and \eqn{\{1,\dots,R\}} should be the
#' decreasing marginals order(pi%*%alpha,decreasing=TRUE) and
#' order(alpha %*% rho,decreasing=TRUE) to
#' alleviate label-switching
#'
#' @inheritParams sample_P_metropolis_classical rho
#'
#' @return A list with sampled trajectories:
#'   `sigma2_array`, `P_array`, `W_array`, `Z_array`, `rho_array`, `alpha_array`.
#' @export
gibbs_sampling_lbm_poisson <- function(
  Y, K, R,
  niter = 50L,
  priors_hyper_params = list(
    etas_0 = rep(2, K),
    gammas_0 = rep(2, R),
    a0 = 1, b0 = 1
  ),
  rho = 1,
  init_Z = NULL,
  init_W = NULL,
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

  pb <- progressr::progressor(niter)

  for (iter in seq(niter)) {
    if (iter %% 10 == 0) {
      message(prefix, "Iter : ", iter, " on ", niter)
    }
    pb(sprintf("%sIter : %d on %d", prefix, iter, niter), class = if (iter %% 10 == 0) "sticky", amount = 0)
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
    pb()
  }
  out_list <- list(W_array = W_array, Z_array = Z_array, rho_array = rho_array, pi_array = pi_array, alpha_array = alpha_array)

  return(posterior::as_draws_array(list_arrays_to_stan(array_list = out_list)))
}

#' Run nchains of LBM Poisson Gibbs Sampler
#'
#' Runs several independent chains of [gibbs_sampling_lbm_poisson()]
#' concurrently, prefixing each chain logs with its index.
#'
#' @param nchains the number of chains to run concurrently
#' @inheritDotParams gibbs_sampling_lbm_poisson
#'
#' @return A list of length `nchains` where each element is the output
#'   of a single call to [gibbs_sampling_lbm_poisson()].
#' @export
chains_gibbs_sampling_lbm_poisson <- function(nchains, ...) {
  out_list <- lapply(seq(nchains), function(i) {
    gibbs_sampling_lbm_poisson(..., prefix = paste0("Chain ", i, " - "))
  }) |> futurize::futurize(seed = TRUE)
  out_draws <- out_list[[1]]
  for (idx in seq_along(out_list[-1])) {
    out_draws <- posterior::bind_draws(out_draws, out_list[[idx + 1]], along = "chain")
  }
  return(out_draws)
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
#' @param known_P NULL (default) or a matrix (\eqn{n_1\times K-1}) indicating the known latent correlated position used to compute probabilities of blocks
#' @param P_sampler a function indicating the sampler to use for the P
#' variable. Default to sample_P_metropolis_trick
#'
#' @inheritParams sample_P_metropolis_classical
#' @inheritParams sample_P_metropolis_trick
#' @inheritParams gibbs_sampling_lbm_poisson
#'
#'
#' @return A list with sampled trajectories:
#'   `sigma2_array`, `P_array`, `W_array`, `Z_array`, `rho_array`, `alpha_array`.
#' @export
gibbs_sampling_lbm_cov_poisson <- function(
  Sigma, Y, init_Z, init_W, K, R,
  niter = 50L, niter_metropolis = 1L,
  priors_hyper_params = list(alpha_0 = 1, beta_0 = 1, gammas_0 = rep(2, R), a0 = 1, b0 = 1),
  rho = 1,
  sigma2_fixed = TRUE,
  known_alpha = NULL,
  known_P = NULL,
  known_W = NULL,
  known_Z = NULL,
  P_sampler = sample_P_metropolis_trick_cpp,
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
  } else if (is.numeric(sigma2_fixed)) {
    message("Using sigma2=", sigma2_fixed)
    current_sigma2 <- sigma2_fixed
  } else {
    current_sigma2 <- 1.0
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
  pb <- progressr::progressor(niter)

  for (iter in seq(niter)) {
    if (iter %% 10 == 0) {
      message(prefix, "Iter : ", iter, " on ", niter)
    }
    pb(sprintf("%sIter : %d on %d", prefix, iter, niter), class = if (iter %% 10 == 0) "sticky", amount = 0)
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
    } else if (is.numeric(sigma2_fixed)) {
      current_sigma2 <- sigma2_fixed
    } else {
      current_sigma2 <- 1.0
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
    pb()
  }
  out_list <- list(sigma2_array = sigma2_array, P_array = P_array, W_array = W_array, Z_array = Z_array, rho_array = rho_array, alpha_array = alpha_array)
  return(posterior::as_draws_array(list_arrays_to_stan(array_list = out_list)))
}

#' Run nchains of the latent phylogenetic Poisson LBM Gibbs sampler
#'
#' Runs several independent chains of [gibbs_sampling_lbm_cov_poisson()]
#' concurrently, prefixing each chain logs with its index.
#'
#' @param nchains the number of chains to run concurrently
#' @inheritDotParams gibbs_sampling_lbm_cov_poisson
#'
#' @return A list of length `nchains` where each element is the output
#'   of a single call to [gibbs_sampling_lbm_cov_poisson()].
#' @export
chains_gibbs_sampling_lbm_cov_poisson <- function(nchains, ...) {
  out_list <- lapply(seq(nchains), function(i) {
    gibbs_sampling_lbm_cov_poisson(..., prefix = paste0("Chain ", i, " - "))
  }) |> futurize::futurize(seed = TRUE)
  out_draws <- out_list[[1]]
  for (idx in seq_along(out_list[-1])) {
    out_draws <- posterior::bind_draws(out_draws, out_list[[idx + 1]], along = "chain")
  }
  return(out_draws)
}
