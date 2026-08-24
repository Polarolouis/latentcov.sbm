# Benchmark of the C++ trick Metropolis sampler against the pure R version.
#
# Measures the wall-clock gain of sample_P_metropolis_trick_cpp (C++)
# over sample_P_metropolis_trick (R) for a SINGLE Metropolis iteration
# (niter_metropolis = 1), i.e. one update pass over all individuals.
#
# Usage:
#   Rscript dev/bench_trick_metropolis.R
#   or source() it inside an interactive session.

library(microbenchmark)

# --- Load both implementations ----------------------------------------------
# R reference implementation and its helpers
source(file.path("R", "simulate.R"))
source(file.path("R", "utils.R"))

# C++ implementation (compiled on the fly)
Rcpp::sourceCpp(file.path("src", "P_samplers.cpp"), rebuild = TRUE)

# --- Data -------------------------------------------------------------------
# One-block-per-species latent structure, as in tests/testthat helpers
K <- 4
npc <- 30
n1 <- npc * K
sigma2 <- 1
niter_metropolis <- 1L # single Metropolis iteration
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

# --- Sanity check: the two versions target the same kernel -------------------
set.seed(123)
out_R <- sample_P_metropolis_trick(
  P = P0, Z = Z0, Sigma = Sigma, sigma2 = sigma2,
  minibatch = minibatch, niter_metropolis = niter_metropolis
)
out_Cpp <- sample_P_metropolis_trick_cpp(
  P = P0, Z = Z0, Sigma = Sigma, sigma2 = sigma2,
  minibatch = minibatch, niter_metropolis = niter_metropolis
)
stopifnot(
  identical(dim(out_R), dim(out_Cpp)),
  all(is.finite(out_R)),
  all(is.finite(out_Cpp))
)

# --- Benchmark ---------------------------------------------------------------
bm <- microbenchmark(
  R   = sample_P_metropolis_trick(
    P = P0, Z = Z0, Sigma = Sigma, sigma2 = sigma2,
    minibatch = minibatch, niter_metropolis = niter_metropolis
  ),
  Cpp = sample_P_metropolis_trick_cpp(
    P = P0, Z = Z0, Sigma = Sigma, sigma2 = sigma2,
    minibatch = minibatch, niter_metropolis = niter_metropolis
  ),
  times = 100,
  unit = "ms"
)

bm_summary <- summary(bm)
cat("\nSingle Metropolis iteration (", n1, " individuals, K = ", K, "):\n", sep = "")
print(bm_summary)

median_R <- bm_summary$median[bm_summary$expr == "R"]
median_Cpp <- bm_summary$median[bm_summary$expr == "Cpp"]
cat(sprintf(
  "\nMedian time: R = %.3f ms | C++ = %.3f ms | C++ is %.1fx faster\n",
  median_R, median_Cpp, median_R / median_Cpp
))

# --- Plot (optional, requires ggplot2 via autoplot) ---------------------------
if (requireNamespace("ggplot2", quietly = TRUE)) {
  g <- ggplot2::autoplot(bm)
  gg_path <- file.path("dev", "bench_trick_metropolis.png")
  ggplot2::ggsave(gg_path, g, width = 6, height = 4)
  cat("Boxplot saved to:", gg_path, "\n")
}
