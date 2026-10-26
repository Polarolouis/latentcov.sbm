# Shared orchestration helpers of the Gibbs samplers
#
# Pure helpers extracted from the duplicated setup/loop/tail blocks of the
# three samplers in gibbs-samplers.R. Each sampler keeps its own explicit loop
# and only delegates the repetitive parts to these functions.

#' Mask-aware data preparation for a Gibbs run
#'
#' Builds the observation mask and replaces the NAs of `Y` with 0 once, so
#' that downstream aggregations (`sufficient_stats_*`) never multiply `NA`
#' values.
#'
#' @param Y the observed count matrix (\eqn{n_1 \times n_2}), with or without NAs
#'
#' @return a list with `Y` (NA-free) and `mask` (\eqn{n_1 \times n_2} 0/1 matrix)
#' @keywords internal
setup_gibbs_data <- function(Y) {
  list(
    Y = replace(Y, is.na(Y), 0),
    mask = (!is.na(Y)) * 1L
  )
}

#' Allocate a trajectory array for a scalar/membership parameter
#'
#' @param niter number of Gibbs iterations
#' @param n number of nodes (rows or columns) tracked
#' @param what a character string naming the tracked parameter (`"Z"`, `"W"`,
#'   `"rho"`, `"pi"`); used to build the column names
#'
#' @return an array of dimension \eqn{niter \times n} filled with `NA`
#' @keywords internal
new_membership_array <- function(niter, n, what = c("Z", "W", "rho", "pi")) {
  what <- match.arg(what)
  label <- switch(what, Z = "Z.", W = "W.", rho = "rho.", pi = "pi.")
  array(
    NA,
    dim = c(niter, n),
    dimnames = list(
      "Iteration" = seq(niter),
      "Parameter" = paste0(label, seq_len(n))
    )
  )
}

#' Allocate a trajectory array for the connectivity matrix `alpha`
#'
#' @param niter number of Gibbs iterations
#' @param K number of row groups
#' @param R number of column groups
#'
#' @return an array of dimension \eqn{niter \times K \times R} filled with `NA`
#' @keywords internal
new_alpha_array <- function(niter, K, R) {
  array(
    NA,
    dim = c(niter, K, R),
    dimnames = list(
      "Iteration" = seq(niter),
      "RowGroup" = paste0("RowGroup", seq_len(K)),
      "ColGroup" = paste0("ColGroup", seq_len(R))
    )
  )
}

#' Allocate a trajectory array for `sigma2`
#'
#' @param niter number of Gibbs iterations
#'
#' @return an array of dimension \eqn{niter \times 1} filled with `NA`
#' @keywords internal
new_sigma2_array <- function(niter) {
  array(
    NA,
    dim = c(niter, 1),
    dimnames = list("Iteration" = seq(niter), "Parameter" = "sigma2")
  )
}

#' Allocate a trajectory array for the latent positions `P`
#'
#' @param niter number of Gibbs iterations
#' @param n1 number of row nodes
#' @param K number of row groups
#'
#' @return an array of dimension \eqn{niter \times n_1 \times (K-1)} filled with `NA`
#' @keywords internal
new_P_array <- function(niter, n1, K) {
  array(
    NA,
    dim = c(niter, n1, K - 1),
    dimnames = list(
      "Iteration" = seq(niter),
      "Individual" = paste0("P", seq_len(n1)),
      "Coordinates" = seq_len(K - 1)
    )
  )
}

#' Allocate a trajectory array for the soft memberships `Z_post_probs`
#'
#' @param niter number of Gibbs iterations
#' @param n1 number of row nodes
#' @param K number of row groups
#'
#' @return an array of dimension \eqn{niter \times n_1 \times K} filled with `NA`
#' @keywords internal
new_Z_post_probs_array <- function(niter, n1, K) {
  array(
    NA,
    dim = c(niter, n1, K),
    dimnames = list(
      "Iteration" = seq(niter),
      "Individual" = paste0("Zpp.", seq_len(n1)),
      "RowGroup" = paste0("RowGroup", seq_len(K))
    )
  )
}

#' Extract the hyperparameters from the prior list
#'
#' @param priors_hyper_params a list of hyperparameters
#' @param latent logical; whether the latent structure priors (`alpha_0`,
#'   `beta_0`, no `etas_0`) are expected
#'
#' @return a list with the hyperparameters relevant to the model flavour
#' @keywords internal
extract_priors <- function(priors_hyper_params, latent = FALSE) {
  out <- list(
    gammas_0 = priors_hyper_params[["gammas_0"]],
    a0 = priors_hyper_params[["a0"]],
    b0 = priors_hyper_params[["b0"]]
  )
  if (latent) {
    out$alpha_0 <- priors_hyper_params[["alpha_0"]]
    out$beta_0 <- priors_hyper_params[["beta_0"]]
  } else {
    out$etas_0 <- priors_hyper_params[["etas_0"]]
  }
  out
}

#' Report the current iteration to the console and progressor
#'
#' @param iter the current iteration
#' @param niter the total number of iterations
#' @param prefix a prefix to prepend to log messages
#' @param pb a progressr progressor created with [progressr::progressor()]
#' @keywords internal
report_iteration <- function(iter, niter, prefix, pb) {
  if (iter %% 10 == 0) {
    message(prefix, "Iter : ", iter, " on ", niter)
  }
  pb(
    sprintf("%sIter : %d on %d", prefix, iter, niter),
    class = if (iter %% 10 == 0) "sticky",
    amount = 0
  )
}


#' Run several independent Gibbs chains concurrently
#'
#' Runs `nchains` independent calls of a single-chain Gibbs sampler
#' concurrently, prefixing each chain logs with its index and saving each
#' chain to its own `auto_save_path`, then binds the chains along the `chain`
#' dimension.
#'
#' @param fun a single-chain Gibbs sampler such as
#'   [gibbs_sampling_lbm_poisson()]
#' @param nchains the number of chains to run concurrently
#' @param auto_save_path the base name from which each chain will be named
#'   `auto_save_path` + `_i.Rds`
#' @param verbose a boolean indicating if the function should print messages.
#'   Defaults to the value of `option("latentcov.sbm.verbose")`.
#' @param ... further arguments passed to `fun`
#'
#' @return a `posterior` draws array with the concatenated chains results
#' @keywords internal
run_chains <- function(
  fun,
  nchains,
  auto_save_path = tempfile("lbm_chain_"),
  verbose = getOption("latentcov.sbm.verbose", default = FALSE),
  ...
) {
  if (verbose) {
    cli::cli_inform(
      "Chains will be saved as {.file {paste0(auto_save_path, '_*.Rds')}}"
    )
  }
  out_list <- lapply(seq(nchains), function(i) {
    fun(
      ...,
      prefix = paste0("Chain ", i, " - "),
      auto_save_path = paste0(auto_save_path, "_", i, ".Rds")
    )
  }) |>
    futurize::futurize(seed = TRUE)
  out_draws <- out_list[[1]]
  for (idx in seq_along(out_list[-1])) {
    out_draws <- posterior::bind_draws(
      out_draws,
      out_list[[idx + 1]],
      along = "chain"
    )
  }
  return(out_draws)
}
