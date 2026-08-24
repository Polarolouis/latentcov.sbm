#' Plot the recovery of the connectivity matrix alpha
#'
#' For each MCMC chain, plots histograms of the sampled `alpha` entries
#' against the true values (dashed vertical lines), after aligning the block
#' labels with [find_permutation_alphas()].
#'
#' @param draws A [posterior::draws_array] with dimensions
#'   iteration \eqn{\times} chain \eqn{\times} parameter, typically produced
#'   by [lbm_results_to_stan_draws()].
#' @param true_alpha The true connectivity matrix (\eqn{K \times R}).
#' @param find_permutation A function estimating row/column permutations
#'   between two alpha matrices. Defaults to [find_permutation_alphas()].
#'
#' @return A [patchwork] plot object.

alpha_recovery_plots <- function(draws, true_alpha, find_permutation = find_permutation_alphas) {
  K <- nrow(true_alpha)
  R <- ncol(true_alpha)
  alpha_var_idx <- grep("^alpha", dimnames(draws)[["variable"]])
  alpha_ref <- matrix(apply(draws[, 1, alpha_var_idx], 3, mean), K, R)
  perms_ref <- find_permutation(alpha_ref, true_alpha)

  true_alpha_df <- true_alpha[perms_ref[["row_perm"]], perms_ref[["col_perm"]]] |>
    reshape2::melt(value.name = "true") |>
    dplyr::mutate(Parameter = paste("alpha[", Var1, ",", Var2, "]", sep = ""))
  old_colors <- unname(unlist(bayesplot::color_scheme_get()))

  all_plots <- lapply(seq(dim(draws)[2]), function(chain_idx) {
    bayesplot::color_scheme_set(scheme = old_colors[c(chain_idx, seq(length(old_colors))[-chain_idx])])
    chain_draws <- draws[, chain_idx, ]
    p <- bayesplot::mcmc_hist(chain_draws, pars = true_alpha_df$Parameter, facet_args = list(nrow = 1)) & ggplot2::lims(x = c(min(true_alpha) - 5, max(true_alpha) + 5))

    p +
      ggplot2::geom_vline(
        data = true_alpha_df,
        ggplot2::aes(xintercept = true, colour = "true", ),
        linetype = "dashed"
      ) + ggplot2::labs(colour = "Parameter type") +
      ggplot2::ggtitle(paste0("Chain #", chain_idx))
  }) |> patchwork::wrap_plots(nrow = posterior::nchains(draws)) + patchwork::plot_annotation(subtitle = "alpha histograms")

  bayesplot::color_scheme_set(scheme = old_colors)
  return(all_plots)
}

#' Plot the recovery of the column block proportions rho
#'
#' For each MCMC chain, plots histograms of the sampled `rho` entries against
#' the true values and the simulated column-block frequencies (dashed vertical
#' lines), after aligning the block labels.
#'
#' @param draws A [posterior::draws_array] with dimensions
#'   iteration \eqn{\times} chain \eqn{\times} parameter.
#' @param true_alpha The true connectivity matrix (\eqn{K \times R}), used to
#'   align the column labels.
#' @param true_rho The true column block proportions (vector of size \eqn{R}).
#' @param true_W_ind A matrix (\eqn{n_2 \times R}) of true column memberships
#'   (one-hot encoding).
#' @param find_permutation A function estimating row/column permutations
#'   between two alpha matrices. Defaults to [find_permutation_alphas()].
#'
#' @return A [patchwork] plot object.

rho_recovery_plots <- function(draws, true_alpha, true_rho, true_W_ind, find_permutation = find_permutation_alphas) {
  R <- ncol(true_alpha)

  alpha_var_idx <- grep("^alpha", dimnames(draws)[["variable"]])
  alpha_ref <- matrix(apply(draws[, 1, alpha_var_idx], 3, mean), ncol = R)
  perms_ref <- find_permutation(alpha_ref, true_alpha)
  true_rho_df <- true_rho[perms_ref[["col_perm"]]] |>
    reshape2::melt(value.name = "true") |>
    dplyr::mutate(Parameter = paste("rho[", seq(R), "]", sep = ""), Type = "true")
  sim_rho_df <- colMeans(true_W_ind)[perms_ref[["col_perm"]]] |>
    reshape2::melt(value.name = "true") |>
    dplyr::mutate(Parameter = paste("rho[", seq(R), "]", sep = ""), Type = "simulated")
  true_rho_df <- rbind(true_rho_df, sim_rho_df)

  old_colors <- unname(unlist(bayesplot::color_scheme_get()))


  all_plots <- lapply(seq_len(dim(draws)[2]), function(chain_idx) {
    bayesplot::color_scheme_set(scheme = old_colors[c(chain_idx, seq(length(old_colors))[-chain_idx])])
    chain_draws <- draws[, chain_idx, ]
    p <- bayesplot::mcmc_hist(chain_draws, pars = true_rho_df$Parameter) & ggplot2::lims(x = c(0, 1))

    p +
      ggplot2::geom_vline(
        data = true_rho_df,
        ggplot2::aes(
          xintercept = true,
          colour = Type
        ),
        linetype = "dashed"
      ) + ggplot2::labs(colour = "Parameter type") +
      ggplot2::ggtitle(paste0("Chain #", chain_idx))
  }) |> patchwork::wrap_plots(nrow = posterior::nchains(draws)) + patchwork::plot_annotation(subtitle = "rho histograms")
  bayesplot::color_scheme_set(old_colors)

  return(all_plots)
}


#' Plot the recovery of the row block proportions pi
#'
#' For each MCMC chain, plots histograms of the sampled `pi` entries against
#' the simulated row-block frequencies (dashed vertical lines), after aligning
#' the block labels.
#'
#' @param draws A [posterior::draws_array] with dimensions
#'   iteration \eqn{\times} chain \eqn{\times} parameter.
#' @param true_alpha The true connectivity matrix (\eqn{K \times R}), used to
#'   align the row labels.
#' @param true_Z_ind A matrix (\eqn{n_1 \times K}) of true row memberships
#'   (one-hot encoding).
#' @param find_permutation A function estimating row/column permutations
#'   between two alpha matrices. Defaults to [find_permutation_alphas()].
#'
#' @return A [patchwork] plot object.
#' @seealso [pi_recovery_plots()], [mean_pi_recovery_plots()]

pi_recovery_plots <- function(draws, true_alpha, true_Z_ind, find_permutation = find_permutation_alphas) {
  K <- nrow(true_alpha)


  alpha_var_idx <- grep("^alpha", dimnames(draws)[["variable"]])
  alpha_ref <- matrix(apply(draws[, 1, alpha_var_idx], 3, mean), nrow = K)
  perms_ref <- find_permutation(alpha_ref, true_alpha)

  sim_pi_df <- colMeans(true_Z_ind)[perms_ref[["row_perm"]]] |>
    reshape2::melt(value.name = "true") |>
    dplyr::mutate(Parameter = paste("pi[", seq(K), "]", sep = ""), Type = "simulated")


  lapply(seq_len(dim(draws)[2]), function(chain_idx) {
    chain_draws <- draws[, chain_idx, ]
    p <- bayesplot::mcmc_hist(chain_draws, pars = sim_pi_df$Parameter) & ggplot2::lims(x = c(0, 1))

    p +
      ggplot2::geom_vline(
        data = sim_pi_df,
        ggplot2::aes(
          xintercept = true,
          colour = Type
        ),
        linetype = "dashed"
      ) + ggplot2::labs(colour = "Parameter type") +
      ggplot2::ggtitle(paste0("Chain #", chain_idx))
  }) |> patchwork::wrap_plots(nrow = posterior::nchains(draws)) + patchwork::plot_annotation(subtitle = "pi histograms")
}


#' Plot the recovery of the posterior mean of pi
#'
#' For each MCMC chain, plots histograms of the posterior mean row-block
#' proportions (computed from the sampled labels) against the simulated
#' row-block frequencies (dashed vertical lines), after aligning the block
#' labels.
#'
#' @param draws A [posterior::draws_array] with dimensions
#'   iteration \eqn{\times} chain \eqn{\times} parameter.
#' @param true_alpha The true connectivity matrix (\eqn{K \times R}), used to
#'   align the row labels.
#' @param true_Z_ind A matrix (\eqn{n_1 \times K}) of true row memberships
#'   (one-hot encoding).
#' @param find_permutation A function estimating row/column permutations
#'   between two alpha matrices. Defaults to [find_permutation_alphas()].
#'
#' @return A [patchwork] plot object.
#' @seealso [pi_recovery_plots()]

mean_pi_recovery_plots <- function(draws, true_alpha, true_Z_ind, find_permutation = find_permutation_alphas) {
  K <- nrow(true_alpha)

  alpha_var_idx <- grep("^alpha", dimnames(draws)[["variable"]])
  Z_var_idx <- grep("^Z", dimnames(draws)[["variable"]])
  alpha_ref <- matrix(apply(draws[, 1, alpha_var_idx], 3, mean), nrow = K)
  perms_ref <- find_permutation(alpha_ref, true_alpha)

  sim_pi_df <- colMeans(true_Z_ind)[perms_ref[["row_perm"]]] |>
    reshape2::melt(value.name = "true") |>
    dplyr::mutate(Parameter = paste("mean_pi", seq(K), sep = "."), Type = "simulated")


  lapply(seq_len(dim(draws)[2]), function(chain_idx) {
    chain_draws <- draws[, chain_idx, ]
    mean_pi_draws <- apply(chain_draws[, , Z_var_idx], 1:2, function(row) {
      Z_fact <- factor(row, levels = seq(K))
      Z_ind <- matrix(0, length(Z_fact), ncol = K)
      Z_ind[cbind(seq_along(Z_fact), as.integer(Z_fact))] <- 1
      colMeans(Z_ind)
    })
    mean_pi_draws <- aperm(mean_pi_draws, c(2, 3, 1))
    dimnames(mean_pi_draws)[[3]] <- paste0("mean_pi.", seq(dim(mean_pi_draws)[3]))
    names(dimnames(mean_pi_draws)) <- c("Iteration", "Chain", "Parameter")
    p <- bayesplot::mcmc_hist(mean_pi_draws, pars = sim_pi_df$Parameter) & ggplot2::lims(x = c(0, 1))

    p +
      ggplot2::geom_vline(
        data = sim_pi_df,
        ggplot2::aes(
          xintercept = true,
          colour = Type
        ),
        linetype = "dashed"
      ) + ggplot2::labs(colour = "Parameter type") +
      ggplot2::ggtitle(paste0("Chain #", chain_idx))
  }) |> patchwork::wrap_plots(nrow = posterior::nchains(draws)) + patchwork::plot_annotation(subtitle = "mean pi histograms")
}

#' Plot the recovery of the latent positions P
#'
#' For each MCMC chain, plots histograms of the sampled latent position
#' coordinates alongside the simulated column-block frequencies.
#'
#' @param draws A [posterior::draws_array] with dimensions
#'   iteration \eqn{\times} chain \eqn{\times} parameter.
#' @param true_alpha The true connectivity matrix (\eqn{K \times R}), used to
#'   align the block labels.
#' @param find_permutation A function estimating row/column permutations
#'   between two alpha matrices. Defaults to [find_permutation_alphas()].
#'
#' @return A [patchwork] plot object.

P_recovery_plots <- function(draws, true_alpha, find_permutation = find_permutation_alphas) {
  K <- nrow(true_alpha)
  alpha_ref <- matrix(apply(draws[, 1, alpha_var_idx], 3, mean), nrow = K)
  perms_ref <- find_permutation(alpha_ref, true_alpha)


  true_P_df <- true_rho[perms_ref[["col_perm"]]] |>
    reshape2::melt(value.name = "true") |>
    dplyr::mutate(Parameter = paste("rho", seq(R), sep = "."), Type = "true")
  sim_rho_df <- colMeans(true_W_ind)[perms_ref[["col_perm"]]] |>
    reshape2::melt(value.name = "true") |>
    dplyr::mutate(Parameter = paste("rho", seq(R), sep = "."), Type = "simulated")
  true_rho_df <- rbind(true_rho_df, sim_rho_df)

  lapply(seq_len(dim(draws)[2]), function(chain_idx) {
    chain_draws <- draws[, chain_idx, ]
    p <- bayesplot::mcmc_hist(chain_draws, pars = true_rho_df$Parameter)

    p +
      ggplot2::geom_vline(
        data = true_rho_df,
        ggplot2::aes(
          xintercept = true,
          colour = Type
        ),
        linetype = "dashed"
      ) + ggplot2::labs(colour = "Parameter type") +
      ggplot2::ggtitle(paste0("Chain #", chain_idx))
  }) |> patchwork::wrap_plots(nrow = nchains)
}
