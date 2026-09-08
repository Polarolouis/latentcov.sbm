# Benchmark of the C++ Metropolis-Hastings P samplers inside the
# chains_gibbs sampler (gibbs_sampling_lbm_cov_poisson).
#
# Two levels are measured with microbenchmark:
#   1. Standalone P samplers: one full Metropolis pass over all n1 rows
#      (niter_metropolis = 1), i.e. exactly the piece that the Gibbs
#      sampler calls each iteration.  Compares:
#        - sample_P_metropolis_classical      (R random walk)
#        - sample_P_metropolis_classical_cpp  (C++ random walk)
#        - sample_P_metropolis_trick          (R clever proposal)
#        - sample_P_metropolis_trick_cpp      (C++ clever proposal)
#   2. Full per-chain Gibbs sampler
#      gibbs_sampling_lbm_cov_poisson with each P_sampler plugged in.
#      chains_gibbs_sampling_lbm_cov_poisson merely lapplies this
#      function over chains (see R/gibbs-samplers.R), so the per-chain
#      wall-clock is the relevant measure of the C++ gain.
#
# Usage:
#   Rscript dev/bench_cpp_metropolis_gibbs.R
#   or source() it inside an interactive session.

library(microbenchmark)

# --- Load both implementations ----------------------------------------------
# R reference implementation and its helpers
devtools::load_all()

# C++ implementation (compiled on the fly)
Rcpp::sourceCpp(file.path("src", "P_samplers.cpp"), rebuild = TRUE)

# --- Data -------------------------------------------------------------------
# Synthetic latent structure as in tests/testthat helpers
K <- 4
R <- 4
npc <- 30
n1 <- npc * K # number of rows
n2 <- 30 # number of columns
sigma2 <- 1
niter_metropolis <- 1L # single Metropolis step per Gibbs iteration

sigma2_fixed <- TRUE
minibatch <- TRUE

rho_within <- 0.9
rho_between <- -1 / (K * npc - 1)

set.seed(1)
Sigma <- matrix(rho_between, nrow = n1, ncol = n1)
same_cluster <- outer(rep(seq_len(K), each = npc), rep(seq_len(K), each = npc), `==`)
Sigma[same_cluster] <- rho_within
diag(Sigma) <- 1

P0 <- matrix(rnorm(n1 * (K - 1)), nrow = n1)
Z0 <- t(sapply(seq_len(n1), function(j) (seq_len(K) == ((j - 1) %% K) + 1) * 1))
W0 <- t(sapply(seq_len(n2), function(j) (seq_len(R) == ((j - 1) %% R) + 1) * 1))

# Synthetic Poisson network Y | Z, W, alpha_true
set.seed(2)
alpha_true <- matrix(seq(0.5, length.out = K * R, by = 0.1), nrow = K)
lambda <- Z0 %*% alpha_true %*% t(W0)
Y <- matrix(rpois(n1 * n2, lambda), nrow = n1, ncol = n2)

priors_hyper_params <- list(alpha_0 = 1, beta_0 = 1, gammas_0 = rep(2, R), etas_0 = rep(2, K), a0 = 1, b0 = 1)

# --- Sanity check: every implementation runs and returns finite values -----
samplers <- list(
  "R classical" = sample_P_metropolis_classical,
  "R trick" = sample_P_metropolis_trick,
  "C++ classical" = sample_P_metropolis_classical_cpp,
  "C++ trick" = sample_P_metropolis_trick_cpp
)

set.seed(123)
for (nm in names(samplers)) {
  out <- samplers[[nm]](
    P = P0, Z = Z0, Sigma = Sigma, sigma2 = sigma2,
    minibatch = minibatch, niter_metropolis = niter_metropolis
  )
  stopifnot(
    identical(dim(out), dim(P0)),
    all(is.finite(out))
  )
}
cat("All P samplers ran correctly (finite output, same dims).\n")

set.seed(123)
for (nm in names(samplers)) {
  res <- suppressMessages(suppressWarnings(
    gibbs_sampling_lbm_cov_poisson(
      Sigma = Sigma, Y = Y, K = K, R = R,
      init_Z = Z0, init_W = W0,
      niter = 3L, niter_metropolis = niter_metropolis,
      sigma2_fixed = sigma2_fixed, P_sampler = samplers[[nm]],
      priors_hyper_params = priors_hyper_params, verbose = FALSE
    )
  ))
  stopifnot(all(is.finite(as.matrix(res))))
}
cat("Full Gibbs sampler ran correctly with each P sampler.\n")

# --- Benchmark 1: standalone P samplers (one Metropolis pass) ---------------
bm_mh <- microbenchmark(
  "R classical" = sample_P_metropolis_classical(
    P = P0, Z = Z0, Sigma = Sigma, sigma2 = sigma2,
    minibatch = minibatch, niter_metropolis = niter_metropolis
  ),
  "R trick" = sample_P_metropolis_trick(
    P = P0, Z = Z0, Sigma = Sigma, sigma2 = sigma2,
    minibatch = minibatch, niter_metropolis = niter_metropolis
  ),
  "C++ classical" = sample_P_metropolis_classical_cpp(
    P = P0, Z = Z0, Sigma = Sigma, sigma2 = sigma2,
    minibatch = minibatch, niter_metropolis = niter_metropolis
  ),
  "C++ trick" = sample_P_metropolis_trick_cpp(
    P = P0, Z = Z0, Sigma = Sigma, sigma2 = sigma2,
    minibatch = minibatch, niter_metropolis = niter_metropolis
  ),
  times = 10,
  unit = "ms"
)

bm_mh_summary <- summary(bm_mh)
cat("\n--- Benchmark 1: standalone P samplers, one Metropolis pass ---\n")
cat(sprintf("(n1 = %d rows, K = %d, niter_metropolis = %d)\n", n1, K, niter_metropolis))
print(bm_mh_summary)

median_mh <- bm_mh_summary$median
speedup <- median_mh[which(names(median_mh) == "R classical")] / median_mh[which(names(median_mh) == "C++ classical")]
cat(sprintf(
  "Median: R classical = %.2f ms | C++ classical = %.2f ms | C++ classical is %.1fx faster\n",
  median_mh["R classical"], median_mh["C++ classical"], speedup
))
speedup_trick <- median_mh[which(names(median_mh) == "R trick")] / median_mh[which(names(median_mh) == "C++ trick")]
cat(sprintf(
  "Median: R trick = %.2f ms | C++ trick = %.2f ms | C++ trick is %.1fx faster\n",
  median_mh["R trick"], median_mh["C++ trick"], speedup_trick
))

# --- Benchmark 2: full per-chain Gibbs sampler ------------------------------
# chains_gibbs_sampling_lbm_cov_poisson runs this in parallel over
# chains (R/gibbs-samplers.R); the per-chain cost is what matters.
gibbs_call <- function(P_sampler) {
  suppressMessages(suppressWarnings(
    gibbs_sampling_lbm_cov_poisson(
      Sigma = Sigma, Y = Y, K = K, R = R,
      init_Z = Z0, init_W = W0,
      niter = 100L, niter_metropolis = niter_metropolis,
      sigma2_fixed = sigma2_fixed, P_sampler = P_sampler,
      priors_hyper_params = priors_hyper_params, verbose = FALSE
    )
  ))
}

bm_gibbs <- microbenchmark(
  "R classical" = gibbs_call(sample_P_metropolis_classical),
  "R trick" = gibbs_call(sample_P_metropolis_trick),
  "C++ classical" = gibbs_call(sample_P_metropolis_classical_cpp),
  "C++ trick" = gibbs_call(sample_P_metropolis_trick_cpp),
  times = 10,
  unit = "ms"
)

bm_gibbs_summary <- summary(bm_gibbs)
cat("\n--- Benchmark 2: full per-chain Gibbs sampler ---\n")
cat(sprintf("(niter = 10, niter_metropolis = %d, n1 = %d)\n", niter_metropolis, n1))
print(bm_gibbs_summary)

median_gibbs <- bm_gibbs_summary$median
speedup_gibbs <- median_gibbs[which(names(median_gibbs) == "R classical")] / median_gibbs[which(names(median_gibbs) == "C++ classical")]
cat(sprintf(
  "Median over a chain: R classical = %.1f ms | C++ classical = %.1f ms | C++ classical is %.1fx faster\n",
  median_gibbs["R classical"], median_gibbs["C++ classical"], speedup_gibbs
))
speedup_gibbs_trick <- median_gibbs[2] / median_gibbs[4]
cat(sprintf(
  "Median over a chain: R trick = %.1f ms | C++ trick = %.1f ms | C++ trick is %.1fx faster\n",
  median_gibbs[2], median_gibbs[4], speedup_gibbs_trick
))

# --- Plots (optional, requires ggplot2) --------------------------------------
if (requireNamespace("ggplot2", quietly = TRUE)) {
  bm_all <- rbind(
    transform(as.data.frame(bm_mh), level = "P sampler only"),
    transform(as.data.frame(bm_gibbs), level = "Gibbs chain")
  )
  bm_all$expr <- factor(bm_all$expr, levels = names(samplers))
  g <- ggplot2::ggplot(bm_all, ggplot2::aes(x = expr, y = time / 1e6, fill = expr)) +
    ggplot2::geom_boxplot() +
    ggplot2::facet_wrap(~level, scales = "free_y") +
    ggplot2::labs(x = NULL, y = "Time (ms)", title = "C++ vs R Metropolis-Hastings P samplers") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle = 15, hjust = 1))
  gg_path <- file.path("dev", "bench_cpp_metropolis_gibbs.png")
  ggplot2::ggsave(gg_path, g, width = 8, height = 4)
  cat("Boxplot saved to:", gg_path, "\n")
}
