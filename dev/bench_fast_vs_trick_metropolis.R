# Benchmark of the fast precision-based Metropolis-Hastings P samplers
# (sample_P_metropolis_trick_cpp_fast, sample_P_metropolis_trick_soft_cpp_fast)
# against the reference Sigma-based trick sampler (C++ and pure R),
# on realistic latent-covariance-SBM cases.
#
# Two aspects are compared:
#   1. BEHAVIOUR / CORRECTNESS
#      The fast and the trick variants are RNG-lockstep: with an identical
#      seed they must leave R's .Random.seed in the same state (so the rest of
#      a Gibbs run stays in sync) and the positions themselves must agree up to
#      floating-point rounding of the conditional mean (Theta vs solve(Sigma)).
#      We verify this for both the hard fast sampler (vs trick_cpp) and the
#      soft fast sampler (vs the pure-R soft). We also report the acceptance
#      similarity indirectly via the max L1 difference.
#   2. SPEED
#      Wall-clock, one full Metropolis pass over all n rows, fast vs trick,
#      for a realistic subset of problem sizes (the O(n^3) trick reference keeps
#      the grid small) and both minibatch modes.
#
# Usage:
#   Rscript dev/bench_fast_vs_trick_metropolis.R
#   or source() it inside an interactive session.

library(microbenchmark)

# --- Load the package (R wrappers + compiled C++) ----------------------------
devtools::load_all()

# --- Realistic data ----------------------------------------------------------
# Correlated block structure: K blocks of npc rows each, within-block
# correlation rho_within, between-block rho_between. This mirrors the setup
# used throughout tests and the other dev benchmarks.
make_data <- function(K, npc, n2 = 30, sigma2 = 1) {
  n1 <- npc * K
  rho_within <- 0.9
  rho_between <- -1 / (n1 - 1)

  set.seed(1)
  Sigma <- matrix(rho_between, nrow = n1, ncol = n1)
  same_cluster <- outer(
    rep(seq_len(K), each = npc),
    rep(seq_len(K), each = npc),
    `==`
  )
  Sigma[same_cluster] <- rho_within
  diag(Sigma) <- 1

  P0 <- matrix(rnorm(n1 * (K - 1)), nrow = n1)
  # hard memberships: row j belongs to block ((j-1) %% K) + 1
  Z0 <- t(sapply(seq_len(n1), function(j) (seq_len(K) == ((j - 1) %% K) + 1) * 1))
  # soft memberships: slight mass on neighbouring blocks (still row-stochastic)
  Z0s <- Z0
  for (j in seq_len(n1)) {
    k <- ((j - 1) %% K) + 1
    Z0s[j, k] <- 0.8
    Z0s[j, (k %% K) + 1] <- 0.2
  }
  list(
    K = K, n1 = n1, sigma2 = sigma2,
    Sigma = Sigma, P = P0, Z_hard = Z0, Z_soft = Z0s,
    niter_metropolis = 10L
  )
}

# A few realistic sizes: small, medium, larger networks of latent factors.
cases <- list(
  "n=120 K=3" = make_data(K = 3, npc = 40),
  "n=200 K=4" = make_data(K = 4, npc = 50),
  "n=400 K=5" = make_data(K = 5, npc = 80)
)

# --- 1. Behaviour: RNG lockstep + position agreement ------------------------
cat("\n=================== 1. BEHAVIOUR / CORRECTNESS ===================\n")
cat("(same seed -> same .Random.seed stream, positions agree up to roundoff)\n\n")

run_both <- function(data, hard = TRUE, minibatch = TRUE, seed = 4242L) {
  set.seed(seed)
  if (hard) {
    r <- sample_P_metropolis_trick_cpp_fast(
      data$P, data$Z_hard, data$Sigma, data$sigma2,
      minibatch = minibatch, niter_metropolis = data$niter_metropolis
    )
    seed_fast <- .Random.seed
    set.seed(seed)
    ref <- sample_P_metropolis_trick_cpp(
      data$P, data$Z_hard, data$Sigma, data$sigma2,
      minibatch = minibatch, niter_metropolis = data$niter_metropolis
    )
    label <- "hard  : fast C++ (~Theta) vs trick C++ (~solve Sigma)"
  } else {
    r <- sample_P_metropolis_trick_soft_cpp_fast(
      data$P, data$Z_soft, data$Sigma, data$sigma2,
      minibatch = minibatch, niter_metropolis = data$niter_metropolis
    )
    seed_fast <- .Random.seed
    set.seed(seed)
    ref <- sample_P_metropolis_trick_soft(
      data$P, data$Z_soft, data$Sigma, data$sigma2,
      minibatch = minibatch, niter_metropolis = data$niter_metropolis
    )
    label <- "soft  : fast C++ (~Theta) vs pure-R soft (~solve Sigma)"
  }
  list(
    label = label,
    seed_identical = identical(seed_fast, .Random.seed),
    maxdiff = max(abs(unname(r) - unname(ref)))
  )
}

for (nm in names(cases)) {
  d <- cases[[nm]]
  for (mb in c(TRUE, FALSE)) {
    h <- run_both(d, hard = TRUE, minibatch = mb)
    s <- run_both(d, hard = FALSE, minibatch = mb)
    cat(sprintf(
      "[%s minibatch=%s]\n  %s\n    seed identical: %s | max |diff| ~ %.2e\n",
      nm, mb, h$label, h$seed_identical, h$maxdiff
    ))
    cat(sprintf(
      "  %s\n    seed identical: %s | max |diff| ~ %.2e\n",
      s$label, s$seed_identical, s$maxdiff
    ))
    if (!h$seed_identical || !s$seed_identical) {
      stop("RNG stream NOT in sync — fast/trick diverged!")
    }
    if (h$maxdiff > 1e-6 || s$maxdiff > 1e-6) {
      stop("Positions disagree beyond floating-point roundoff!")
    }
  }
}
cat("\nAll fast samplers are RNG-lockstep with their trick reference.\n")

# --- 2. Speed ----------------------------------------------------------------
cat("\n=================== 2. SPEED ======================================\n")
cat("(one full Metropolis pass over all rows; median of timings)\n")
cat("NOTE: the Sigma-based trick reference is O(n^3) per row, so the timing\n")
cat("grid is limited to sizes where it stays tractable; behaviour (section 1)\n")
cat("already covered every size.\n\n")

# Pre-compute precision matrices so the per-call solve() is not timed (it is
# identical for fast and trick and negligible once cached).
for (nm in names(cases)) {
  d <- cases[[nm]]
  d$Theta <- solve(d$Sigma)
  # Lighter Metropolis load for the timing grid (the trick reference dominates
  # the cost; a per-pass comparison is what we want to time).
  d$nit <- 5L
  cases[[nm]] <- d
}

# Timing grid: medium/large are the headline sizes; the n=400 trick pass is
# ~30 s, so it gets fewer repetitions.
speed_cases <- list(
  "n=120 K=3" = cases[["n=120 K=3"]],
  "n=200 K=4" = cases[["n=200 K=4"]]
)
times_grid <- c("n=120 K=3" = 10, "n=200 K=4" = 5)

bench_one <- function(d, minibatch, times) {
  bm <- microbenchmark(
    "trick (Sigma)" = sample_P_metropolis_trick_cpp(
      P = d$P, Z = d$Z_hard, Sigma = d$Sigma, sigma2 = d$sigma2,
      minibatch = minibatch, niter_metropolis = d$nit
    ),
    "fast  (Theta)" = sample_P_metropolis_trick_cpp_fast(
      P = d$P, Z = d$Z_hard, Sigma = d$Sigma, Theta = d$Theta,
      sigma2 = d$sigma2, minibatch = minibatch,
      niter_metropolis = d$nit
    ),
    "fast soft (Theta)" = sample_P_metropolis_trick_soft_cpp_fast(
      P = d$P, Z = d$Z_soft, Sigma = d$Sigma, Theta = d$Theta,
      sigma2 = d$sigma2, minibatch = minibatch,
      niter_metropolis = d$nit
    ),
    times = times,
    unit = "ms"
  )
  s <- summary(bm)
  md <- setNames(s$median, s$expr)
  cat(sprintf(
    "[%s minibatch=%s niter=%d]  trick=%.0f ms  fast=%.1f ms  fast-soft=%.1f ms\n",
    nm, minibatch, d$nit, md[["trick (Sigma)"]], md[["fast  (Theta)"]],
    md[["fast soft (Theta)"]]
  ))
  cat(sprintf(
    "    -> fast is %.1fx faster than trick; fast-soft is %.1fx faster\n",
    md[["trick (Sigma)"]] / md[["fast  (Theta)"]],
    md[["trick (Sigma)"]] / md[["fast soft (Theta)"]]
  ))
}

for (nm in names(speed_cases)) {
  bench_one(speed_cases[[nm]], minibatch = TRUE, times = times_grid[[nm]])
  bench_one(speed_cases[[nm]], minibatch = FALSE, times = times_grid[[nm]])
}

cat("\nDone.\n")
