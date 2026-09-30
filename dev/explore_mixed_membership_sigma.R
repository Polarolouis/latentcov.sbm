# Compare inference of gibbs_sampling_lbm_cov_poisson_mixed_membership
# using the TRUE (non-diagonal, block-correlated) Sigma vs a misspecified
# diagonal Sigma, on data simulated from the mixed-membership LBM.
#
# Run from the package root with:
#   pkgload::load_all(".")

set.seed(5153)
suppressPackageStartupMessages({ library(posterior); library(tidyverse) })

# ---- Parameters ------------------------------------------------------------
n1 <- 60; n2 <- 70; K <- 3; R <- 2
rho_w <- 0.88; rho_b <- 0.05; Q <- 3; m <- n1 / Q
sigma2_true <- 0.25
alpha_true <- matrix(c(2.5, 0.4, 0.4, 5.0, 1.0, 1.0), K, R)
niter <- 120; niter_metropolis <- 8; warmup <- 40

# ---- Non-diagonal block-correlation among row nodes ------------------------
Sig <- matrix(rho_b, n1, n1)
for (q in seq(Q)) Sig[((q - 1) * m + 1):(q * m), ((q - 1) * m + 1):(q * m)] <- rho_w
diag(Sig) <- 1
stopifnot(min(eigen(Sig, symmetric = TRUE)$values) > 1e-6)
Sigma_true <- Sig
Sigma_diag <- diag(n1)

# ---- Simulate latent positions, memberships, W, alpha and counts -----------
P_true <- t(mvtnorm::rmvnorm(K - 1, mean = rep(0, n1), sigma = sigma2_true * Sigma_true))[, 1:(K - 1)]
probs_true <- ilrInv(P_true)         # soft memberships
Z_true <- probs_true
rho_true <- as.vector(MCMCpack::rdirichlet(1, rep(2, R)))
W_true <- t(sapply(sample.int(R, n2, prob = rho_true, replace = TRUE),
                   function(l) as.integer(seq(R) == l)))
Lambda <- Z_true %*% alpha_true %*% t(W_true)
Y <- matrix(rpois(n1 * n2, Lambda), n1, n2)

# ---- Run the mixed-membership sampler with both Sigma choices --------------
priors <- list(alpha_0 = 1, beta_0 = 1, gammas_0 = rep(2, R), a0 = 1, b0 = 1)
fit_true <- gibbs_sampling_lbm_cov_poisson_mixed_membership(
  Sigma = Sigma_true, Y = Y, init_Z = Z_true, init_W = NULL, K = K, R = R,
  niter = niter, niter_metropolis = niter_metropolis, sigma2_fixed = TRUE,
  priors_hyper_params = priors, tol = 1e-6, P_sampler = sample_P_metropolis_trick_soft)
fit_diag <- gibbs_sampling_lbm_cov_poisson_mixed_membership(
  Sigma = Sigma_diag, Y = Y, init_Z = Z_true, init_W = NULL, K = K, R = R,
  niter = niter, niter_metropolis = niter_metropolis, sigma2_fixed = TRUE,
  priors_hyper_params = priors, tol = 1e-6, P_sampler = sample_P_metropolis_trick_soft)

# ---- Per-iteration membership recovery (label-switch resolved) -------------
z_true_hard <- apply(probs_true, 1, which.max)
per_iter_acc <- function(fit, n, K, z_true, niter, warmup) {
  zp <- as_draws_matrix(subset_draws(fit, variable = "Z_post_probs"))
  its <- (warmup + 1):niter
  sapply(its, function(it) {
    Z <- matrix(0, n, K)
    for (k in seq(K)) for (i in seq(n)) Z[i, k] <- zp[it, (k - 1) * n + i]
    z_hat <- apply(Z, 1, which.max)
    tab <- table(factor(z_true, 1:K), factor(z_hat, 1:K))
    p <- as.integer(clue::solve_LSAP(tab, maximum = TRUE))
    sum(diag(tab[, p])) / n
  })
}
a_true <- per_iter_acc(fit_true, n1, K, z_true_hard, niter, warmup)
a_diag <- per_iter_acc(fit_diag, n1, K, z_true_hard, niter, warmup)

# ---- Posterior decisiveness (max membership probability) -------------------
zp_tm <- as_draws_matrix(subset_draws(fit_true, variable = "Z_post_probs"))
zp_dm <- as_draws_matrix(subset_draws(fit_diag, variable = "Z_post_probs"))
maxZ <- function(zp, n, K) {
  M <- matrix(0, n, K)
  for (k in seq(K)) for (i in seq(n)) M[i, k] <- mean(zp[(warmup + 1):niter, (k - 1) * n + i])
  apply(M / rowSums(M), 1, max)
}
conf_true <- maxZ(zp_tm, n1, K); conf_diag <- maxZ(zp_dm, n1, K)

# ---- Summary ---------------------------------------------------------------
cat("=== Membership recovery (per-iteration, best label permutation) ===\n")
cat(sprintf("  true Sigma: mean %.3f (sd %.3f)\n", mean(a_true), sd(a_true)))
cat(sprintf("  diag Sigma: mean %.3f (sd %.3f)\n\n", mean(a_diag), sd(a_diag)))
cat("=== Posterior membership decisiveness (mean max Z) ===\n")
cat(sprintf("  true Sigma: %.3f\n  diag Sigma: %.3f\n", mean(conf_true), mean(conf_diag)))

# ---- Plot ------------------------------------------------------------------
tibble(
  Sigma = factor(c(rep("True (non-diagonal)", length(a_true)),
                   rep("Diagonal", length(a_diag))),
                 levels = c("True (non-diagonal)", "Diagonal")),
  accuracy = c(a_true, a_diag)
) |>
  ggplot(aes(Sigma, accuracy, fill = Sigma)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.6) +
  geom_hline(yintercept = 1 / K, linetype = 2, colour = "grey40") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Row covariance used in the sampler",
       y = "Per-iteration membership recovery accuracy") +
  theme_minimal()
