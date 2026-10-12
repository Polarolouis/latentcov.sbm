# Metropolis-Hastings samplers for P: shape/finiteness invariants and branches.

.P <- P_fx
.sigma2 <- 0.5
.niter <- 5L

# names of the R (dimnames-bearing) samplers only
r_samplers <- list(
  classical_r = function() sample_P_metropolis_classical(.P, Z_fx, Sigma_fx, .sigma2, minibatch = TRUE, niter_metropolis = .niter),
  trick_r = function() sample_P_metropolis_trick(.P, Z_fx, Sigma_fx, .sigma2, minibatch = TRUE, niter_metropolis = .niter),
  trick_soft = function() sample_P_metropolis_trick_soft(.P, ilrInv(.P), Sigma_fx, .sigma2, minibatch = TRUE, niter_metropolis = .niter)
)

cpp_samplers <- list(
  classical_cpp = function() sample_P_metropolis_classical_cpp(.P, Z_fx, Sigma_fx, .sigma2, minibatch = TRUE, niter_metropolis = .niter),
  trick_cpp = function() sample_P_metropolis_trick_cpp(.P, Z_fx, Sigma_fx, .sigma2, minibatch = TRUE, niter_metropolis = .niter),
  trick_cpp_block = function() sample_P_metropolis_trick_cpp_block(.P, Z_fx, Sigma_fx, .sigma2, minibatch = TRUE, niter_metropolis = .niter, block_size = 2L)
)

test_that("all MH samplers preserve dims and finiteness", {
  for (spec in c(r_samplers, cpp_samplers)) {
    got <- spec()
    expect_equal(dim(got), dim(.P), label = paste(deparse1(spec), collapse = ""))
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
    expect_equal(unname(rowSums(probs)), rep(1, nrow(probs)), tolerance = 1e-10, label = paste(deparse1(spec), collapse = ""))
  }
})

test_that("minibatch = FALSE (lexicographic order) runs and preserves shape", {
  out <- sample_P_metropolis_trick(.P, Z_fx, Sigma_fx, .sigma2, minibatch = FALSE, niter_metropolis = 3L)
  expect_equal(dim(out), dim(.P))
  out2 <- sample_P_metropolis_classical(.P, Z_fx, Sigma_fx, .sigma2, minibatch = FALSE, niter_metropolis = 3L)
  expect_equal(dim(out2), dim(.P))
  out3 <- sample_P_metropolis_trick_cpp(.P, Z_fx, Sigma_fx, .sigma2, minibatch = FALSE, niter_metropolis = 3L)
  expect_equal(dim(out3), dim(.P))
  out4 <- sample_P_metropolis_classical_cpp(.P, Z_fx, Sigma_fx, .sigma2, minibatch = FALSE, niter_metropolis = 3L)
  expect_equal(dim(out4), dim(.P))
})

test_that("the classical sampler reports the acceptance rate when verbose", {
  expect_message(
    sample_P_metropolis_classical(.P, Z_fx, Sigma_fx, .sigma2, minibatch = TRUE, niter_metropolis = 2L, verbose = TRUE),
    "Mean accepted rate"
  )
  expect_silent(
    sample_P_metropolis_classical(.P, Z_fx, Sigma_fx, .sigma2, minibatch = TRUE, niter_metropolis = 2L, verbose = FALSE)
  )
})

test_that("the block sampler handles block sizes that do not divide n", {
  out <- sample_P_metropolis_trick_cpp_block(.P, Z_fx, Sigma_fx, .sigma2, minibatch = FALSE, niter_metropolis = 3L, block_size = 4L)
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
  out_trick <- sample_P_metropolis_trick(P_sim, Z_sim, Sigma_fx, sigma2 = 0.01, minibatch = TRUE, niter_metropolis = 20L)
  set.seed(703L)
  out_cpp <- sample_P_metropolis_trick_cpp(P_sim, Z_sim, Sigma_fx, sigma2 = 0.01, minibatch = TRUE, niter_metropolis = 20L)
  expect_equal(unname(out_trick), P_sim, tolerance = 0.2)
  expect_equal(unname(out_cpp), P_sim, tolerance = 0.2)
})
