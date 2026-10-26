# Metropolis-Hastings samplers for P: shape/finiteness invariants and branches.

.P <- P_fx
.sigma2 <- 0.5
.niter <- 5L

# names of the R (dimnames-bearing) samplers only
r_samplers <- list(
  classical_r = function() {
    sample_P_metropolis_classical(
      .P,
      Z_fx,
      Sigma_fx,
      .sigma2,
      minibatch = TRUE,
      niter_metropolis = .niter
    )
  },
  trick_r = function() {
    sample_P_metropolis_trick(
      .P,
      Z_fx,
      Sigma_fx,
      .sigma2,
      minibatch = TRUE,
      niter_metropolis = .niter
    )
  },
  trick_soft = function() {
    sample_P_metropolis_trick_soft(
      .P,
      ilrInv(.P),
      Sigma_fx,
      .sigma2,
      minibatch = TRUE,
      niter_metropolis = .niter
    )
  }
)

cpp_samplers <- list(
  classical_cpp = function() {
    sample_P_metropolis_classical_cpp(
      .P,
      Z_fx,
      Sigma_fx,
      .sigma2,
      minibatch = TRUE,
      niter_metropolis = .niter
    )
  },
  trick_cpp = function() {
    sample_P_metropolis_trick_cpp(
      .P,
      Z_fx,
      Sigma_fx,
      .sigma2,
      minibatch = TRUE,
      niter_metropolis = .niter
    )
  },
  trick_cpp_block = function() {
    sample_P_metropolis_trick_cpp_block(
      .P,
      Z_fx,
      Sigma_fx,
      .sigma2,
      minibatch = TRUE,
      niter_metropolis = .niter,
      block_size = 2L
    )
  },
  trick_cpp_fast = function() {
    sample_P_metropolis_trick_cpp_fast(
      .P,
      Z_fx,
      Sigma_fx,
      .sigma2,
      minibatch = TRUE,
      niter_metropolis = .niter
    )
  },
  trick_soft_cpp_fast = function() {
    sample_P_metropolis_trick_soft_cpp_fast(
      .P,
      Z_fx_soft,
      Sigma_fx,
      .sigma2,
      minibatch = TRUE,
      niter_metropolis = .niter
    )
  }
)

test_that("all MH samplers preserve dims and finiteness", {
  for (spec in c(r_samplers, cpp_samplers)) {
    got <- spec()
    expect_equal(
      dim(got),
      dim(.P),
      label = paste(deparse1(spec), collapse = "")
    )
    expect_false(anyNA(got))
    expect_true(all(is.finite(got)))
  }
})

test_that("R samplers keep the per-row dimnames", {
  for (spec in r_samplers) {
    got <- spec()
    expect_equal(names(dimnames(got)), c("Individual", "Coordinates"))
    expect_equal(dimnames(got)[[1]], paste0("P", seq_len(nrow(.P))))
  }
})

test_that("the sampled positions still map to the simplex through ilrInv", {
  for (spec in c(r_samplers, cpp_samplers)) {
    got <- spec()
    probs <- ilrInv(got)
    expect_equal(
      unname(rowSums(probs)),
      rep(1, nrow(probs)),
      tolerance = 1e-10,
      label = paste(deparse1(spec), collapse = "")
    )
  }
})

test_that("minibatch = FALSE (lexicographic order) runs and preserves shape", {
  out <- sample_P_metropolis_trick(
    .P,
    Z_fx,
    Sigma_fx,
    .sigma2,
    minibatch = FALSE,
    niter_metropolis = 3L
  )
  expect_equal(dim(out), dim(.P))
  out2 <- sample_P_metropolis_classical(
    .P,
    Z_fx,
    Sigma_fx,
    .sigma2,
    minibatch = FALSE,
    niter_metropolis = 3L
  )
  expect_equal(dim(out2), dim(.P))
  out3 <- sample_P_metropolis_trick_cpp(
    .P,
    Z_fx,
    Sigma_fx,
    .sigma2,
    minibatch = FALSE,
    niter_metropolis = 3L
  )
  expect_equal(dim(out3), dim(.P))
  out4 <- sample_P_metropolis_classical_cpp(
    .P,
    Z_fx,
    Sigma_fx,
    .sigma2,
    minibatch = FALSE,
    niter_metropolis = 3L
  )
  expect_equal(dim(out4), dim(.P))
  out5 <- sample_P_metropolis_trick_cpp_fast(
    .P,
    Z_fx,
    Sigma_fx,
    .sigma2,
    minibatch = FALSE,
    niter_metropolis = 3L
  )
  expect_equal(dim(out5), dim(.P))
  out6 <- sample_P_metropolis_trick_soft_cpp_fast(
    .P,
    Z_fx_soft,
    Sigma_fx,
    .sigma2,
    minibatch = FALSE,
    niter_metropolis = 3L
  )
  expect_equal(dim(out6), dim(.P))
})

test_that("the classical sampler reports the acceptance rate when verbose", {
  expect_message(
    sample_P_metropolis_classical(
      .P,
      Z_fx,
      Sigma_fx,
      .sigma2,
      minibatch = TRUE,
      niter_metropolis = 2L,
      verbose = TRUE
    ),
    "Mean accepted rate"
  )
  expect_silent(
    sample_P_metropolis_classical(
      .P,
      Z_fx,
      Sigma_fx,
      .sigma2,
      minibatch = TRUE,
      niter_metropolis = 2L,
      verbose = FALSE
    )
  )
})

test_that("the block sampler handles block sizes that do not divide n", {
  out <- sample_P_metropolis_trick_cpp_block(
    .P,
    Z_fx,
    Sigma_fx,
    .sigma2,
    minibatch = FALSE,
    niter_metropolis = 3L,
    block_size = 4L
  )
  expect_equal(dim(out), dim(.P))
})

test_that("trick samplers barely move a model-consistent P", {
  # A P truly drawn from the model is (almost) a fixed point of the clever
  # proposal: proposals sit at the conditional means, which coincide with the
  # rows themselves up to the (small) shrinkage.
  set.seed(702L)
  P_sim <- simulate_P(K = 3L, Sigma = Sigma_fx, sigma2 = 0.01)
  Z_sim <- .one_hot(apply(ilrInv(P_sim), 1, which.max), 3)
  set.seed(703L)
  out_trick <- sample_P_metropolis_trick(
    P_sim,
    Z_sim,
    Sigma_fx,
    sigma2 = 0.01,
    minibatch = TRUE,
    niter_metropolis = 20L
  )
  set.seed(703L)
  out_cpp <- sample_P_metropolis_trick_cpp(
    P_sim,
    Z_sim,
    Sigma_fx,
    sigma2 = 0.01,
    minibatch = TRUE,
    niter_metropolis = 20L
  )
  expect_equal(unname(out_trick), P_sim, tolerance = 0.2)
  expect_equal(unname(out_cpp), P_sim, tolerance = 0.2)
})

test_that("fast samplers use the same near-fixed-point behaviour", {
  # The fast variant draws from the same full conditional (through Theta), so
  # a model-consistent P is also an approximate fixed point.
  set.seed(712L)
  P_sim <- simulate_P(K = 3L, Sigma = Sigma_fx, sigma2 = 0.01)
  Z_sim <- .one_hot(apply(ilrInv(P_sim), 1, which.max), 3)
  set.seed(713L)
  out_fast <- sample_P_metropolis_trick_cpp_fast(
    P_sim,
    Z_sim,
    Sigma_fx,
    sigma2 = 0.01,
    minibatch = TRUE,
    niter_metropolis = 20L
  )
  expect_equal(unname(out_fast), unname(P_sim), tolerance = 0.25)
})

test_that("the fast samplers share the RNG stream of the Sigma-based trick", {
  # The hard fast sampler must consume exactly the same RNG draws as the
  # Sigma-based `sample_P_metropolis_trick_cpp` (one Rcpp::rnorm per coordinate
  # plus one runif per Metropolis step). With an identical starting seed the
  # two therefore end on the same `.Random.seed` and, because the full
  # conditionals coincide up to floating-point round-off, produce
  # nearly-identical positions.
  Sigma <- Sigma_fx
  run_both <- function(seed) {
    set.seed(seed)
    trick <- latentcov.sbm:::sample_P_metropolis_trick_cpp(
      .P,
      Z_fx,
      Sigma,
      .sigma2,
      minibatch = FALSE,
      niter_metropolis = .niter
    )
    seed_trick <- .Random.seed
    set.seed(seed)
    fast <- sample_P_metropolis_trick_cpp_fast(
      .P,
      Z_fx,
      Sigma,
      .sigma2,
      minibatch = FALSE,
      niter_metropolis = .niter
    )
    list(
      seed_identical = identical(seed_trick, .Random.seed),
      maxdiff = max(abs(unname(trick) - unname(fast)))
    )
  }
  res <- run_both(123L)
  expect_true(res$seed_identical)
  expect_lt(res$maxdiff, 1e-8)
})

test_that("the fast samplers share the RNG stream of the Sigma-based trick and it holds for niter = 50L", {
  # The hard fast sampler must consume exactly the same RNG draws as the
  # Sigma-based `sample_P_metropolis_trick_cpp` (one Rcpp::rnorm per coordinate
  # plus one runif per Metropolis step). With an identical starting seed the
  # two therefore end on the same `.Random.seed` and, because the full
  # conditionals coincide up to floating-point round-off, produce
  # nearly-identical positions.
  Sigma <- Sigma_fx
  run_both <- function(seed) {
    set.seed(seed)
    trick <- latentcov.sbm:::sample_P_metropolis_trick_cpp(
      .P,
      Z_fx,
      Sigma,
      .sigma2,
      minibatch = FALSE,
      niter_metropolis = 50L
    )
    seed_trick <- .Random.seed
    set.seed(seed)
    fast <- sample_P_metropolis_trick_cpp_fast(
      .P,
      Z_fx,
      Sigma,
      .sigma2,
      minibatch = FALSE,
      niter_metropolis = 50L
    )
    list(
      seed_identical = identical(seed_trick, .Random.seed),
      maxdiff = max(abs(unname(trick) - unname(fast)))
    )
  }
  res <- run_both(123L)
  expect_true(res$seed_identical)
  expect_lt(res$maxdiff, 1e-8)
})

test_that("both fast samplers draw per-coordinate normals plus one uniform", {
  # RNG-consumption invariant: over `n * niter` Metropolis steps the fast
  # samplers (hard and soft) must advance R's RNG exactly as `m` standard
  # normals followed by one uniform per step, i.e. in lockstep with the
  # Sigma-based trick samplers.
  Sigma <- Sigma_fx
  m <- ncol(.P)
  niter <- 5L
  n_steps <- nrow(.P) * niter
  expected <- function(seed) {
    set.seed(seed)
    for (i in seq_len(n_steps)) {
      invisible(rnorm(m))
      invisible(runif(1))
    }
    .Random.seed
  }
  set.seed(321L)
  invisible(sample_P_metropolis_trick_cpp_fast(
    .P,
    Z_fx,
    Sigma,
    .sigma2,
    minibatch = FALSE,
    niter_metropolis = niter
  ))
  expect_identical(.Random.seed, expected(321L))
  set.seed(322L)
  invisible(sample_P_metropolis_trick_soft_cpp_fast(
    .P,
    Z_fx_soft,
    Sigma,
    .sigma2,
    minibatch = FALSE,
    niter_metropolis = niter
  ))
  expect_identical(.Random.seed, expected(322L))
})

test_that("the soft R sampler is RNG-lockstep with the fast soft C++ sampler", {
  # `sample_P_metropolis_trick_soft` is a pure-R implementation that consumes
  # R's RNG exactly like `sample_P_metropolis_trick_soft_cpp_fast` (same
  # permutation and per-step draws), so with an identical seed the two must
  # leave R's RNG in the same state and produce positions that agree up to
  # floating-point rounding of the conditional mean.
  Sigma <- Sigma_fx
  run_both <- function(seed) {
    set.seed(seed)
    r <- sample_P_metropolis_trick_soft(
      .P,
      Z_fx_soft,
      Sigma,
      .sigma2,
      minibatch = TRUE,
      niter_metropolis = .niter
    )
    seed_r <- .Random.seed
    set.seed(seed)
    fast <- sample_P_metropolis_trick_soft_cpp_fast(
      .P,
      Z_fx_soft,
      Sigma,
      .sigma2,
      minibatch = TRUE,
      niter_metropolis = .niter
    )
    list(
      seed_identical = identical(seed_r, .Random.seed),
      maxdiff = max(abs(unname(r) - unname(fast)))
    )
  }
  res <- run_both(4242L)
  expect_true(res$seed_identical)
  expect_lt(res$maxdiff, 1e-8)
})
