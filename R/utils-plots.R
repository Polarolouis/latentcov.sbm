alpha_recovery_plots <- function(draws, true_alpha, find_permutation = find_permutation_alphas) {
  K <- nrow(true_alpha)
  R <- ncol(true_alpha)
  alpha_var_idx <- grep("^alpha", dimnames(draws)[["variable"]])
  alpha_ref <- matrix(apply(draws[, 1, alpha_var_idx], 3, mean), K, R)
  perms_ref <- find_permutation(alpha_ref, true_alpha)

  true_alpha_df <- true_alpha[perms_ref[["row_perm"]], perms_ref[["col_perm"]]] |>
    reshape2::melt(value.name = "true") |>
    dplyr::mutate(Parameter = paste("alpha[", Var1, ",", Var2, "]", sep = ""))

  lapply(seq(dim(draws)[2]), function(chain_idx) {
    chain_draws <- draws[, chain_idx, ]
    p <- bayesplot::mcmc_hist(chain_draws, pars = true_alpha_df$Parameter, facet_args = list(nrow = 1)) & ggplot2::lims(x = c(min(true_alpha) - 5, max(true_alpha) + 5))

    p +
      ggplot2::geom_vline(
        data = true_alpha_df,
        ggplot2::aes(xintercept = true, colour = "true", ),
        linetype = "dashed"
      ) + ggplot2::labs(colour = "Parameter type") +
      ggplot2::ggtitle(paste0("Chain N°", chain_idx))
  }) |> patchwork::wrap_plots(nrow = posterior::nchains(draws)) + patchwork::plot_annotation(subtitle = "alpha histograms")
}

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

  lapply(seq_len(dim(draws)[2]), function(chain_idx) {
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
      ggplot2::ggtitle(paste0("Chain N°", chain_idx))
  }) |> patchwork::wrap_plots(nrow = posterior::nchains(draws)) + patchwork::plot_annotation(subtitle = "rho histograms")
}


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
      ggplot2::ggtitle(paste0("Chain N°", chain_idx))
  }) |> patchwork::wrap_plots(nrow = posterior::nchains(draws)) + patchwork::plot_annotation(subtitle = "pi histograms")
}


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
      ggplot2::ggtitle(paste0("Chain N°", chain_idx))
  }) |> patchwork::wrap_plots(nrow = posterior::nchains(draws)) + patchwork::plot_annotation(subtitle = "mean pi histograms")
}

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
    p <- mcmc_hist(chain_draws, pars = true_rho_df$Parameter)

    p +
      geom_vline(
        data = true_rho_df,
        aes(
          xintercept = true,
          colour = Type
        ),
        linetype = "dashed"
      ) + labs(colour = "Parameter type") +
      ggtitle(paste0("Chain N°", chain_idx))
  }) |> patchwork::wrap_plots(nrow = nchains)
}
