# --- compositional transforms (ilr, clr, stick-breaking, Psi basis) -----------

test_that("default_Psi_function builds an orthonormal, centered basis", {
  K <- 4L
  Psi <- default_Psi_function(K)
  expect_equal(dim(Psi), c(K - 1, K))
  inprod <- Psi %*% t(Psi)
  expect_equal(inprod, diag(K - 1), tolerance = 1e-12)
  expect_equal(rowSums(Psi), rep(0, K - 1), tolerance = 1e-12)
  expect_equal(default_Psi_function_cpp(K), default_Psi_function(K), tolerance = 1e-12)
})

test_that("clr centers each row of log-counts", {
  x <- matrix(c(1, 2, 3, 10, 20, 30), 2, 3, byrow = TRUE)
  out <- clr(x)
  expect_equal(rowSums(out), rep(0, 2), tolerance = 1e-12)
  expect_equal(out[1, ], log(x[1, ]) - mean(log(x[1, ])), tolerance = 1e-12)
})

test_that("ilr and ilrInv are inverse on the simplex", {
  set.seed(123L)
  probs <- apply(matrix(runif(4 * 3), 4, 3), 1, function(w) w / sum(w))
  probs <- t(probs)
  roundtrip <- ilrInv(ilr(probs))
  expect_equal(roundtrip, probs, tolerance = 1e-12)
  expect_equal(rowSums(roundtrip), rep(1, nrow(probs)), tolerance = 1e-12)
})

test_that("ilr casts vector input to a one-row matrix with a warning", {
  expect_warning(ilr(c(1, 2, 3, 4)))
  expect_equal(ilr(matrix(c(1, 2, 3, 4), nrow = 1)), suppressWarnings(ilr(c(1, 2, 3, 4))))
})

test_that("ilrInv produces row-stochastic probabilities", {
  set.seed(456L)
  z <- matrix(rnorm(6 * 2), 6, 2)
  probs <- ilrInv(z)
  expect_equal(dim(probs), c(6, 3))
  expect_equal(rowSums(probs), rep(1, 6), tolerance = 1e-12)
  expect_true(all(probs > 0))
})

test_that("stick_breaking follows its recursive definition", {
  x <- c(0, 0, 0)
  expect_equal(stick_breaking(x), c(1 / 3, 1 / 4), tolerance = 1e-12)
})

test_that("cat_dist_ilr_given_Pi matches ilrInv at the index", {
  Pi <- c(0.3, -0.7)
  for (Zi in 1:3) {
    expect_equal(cat_dist_ilr_given_Pi(Zi, Pi), ilrInv(matrix(Pi, 1, 2))[1, Zi], tolerance = 1e-12)
  }
})

# --- C++ / R equivalence for the pivot transform ------------------------------

test_that("ilrInvcpp agrees with ilrInv", {
  set.seed(789L)
  z <- matrix(rnorm(4 * 2), 4, 2)
  expect_equal(ilrInvcpp(z), ilrInv(z), tolerance = 1e-10)
})

test_that("ilrInv_cpp exposes the log option", {
  set.seed(101L)
  z <- matrix(rnorm(3 * 2), 3, 2)
  probs <- ilrInv_cpp(z, default_Psi_function(3), log = FALSE)
  logprobs <- ilrInv_cpp(z, default_Psi_function(3), log = TRUE)
  expect_equal(dim(probs), c(3, 3))
  expect_equal(rowSums(probs), rep(1, 3), tolerance = 1e-10)
  expect_equal(logprobs, log(probs), tolerance = 1e-10)
})

# --- inverse-gamma helpers ----------------------------------------------------

test_that("posterior_param_inv_gamma adds observed contributions", {
  P <- matrix(c(1, 0, 0, 1), 2, 2)
  Theta <- diag(2) * 2
  res <- posterior_param_inv_gamma(alpha_0 = 1, beta_0 = 3, P = P, Theta = Theta)
  expect_equal(res$alpha, 1 + nrow(P) / 2)
  # t(P) %*% Theta %*% P = 2 * t(P) P, diag = 2*1 + 2*1 = 4; sum(diag) = 4; *0.5 = 2
  expect_equal(res$beta, 3 + 2)
})

test_that("sample_inv_gamma_rate draws positive values", {
  set.seed(202L)
  draws <- replicate(200, sample_inv_gamma_rate(shape = 2, rate = 1))
  expect_true(all(is.finite(draws)))
  expect_true(all(draws > 0))
})

test_that("sample_inv_gamma_rate respects mean = rate/(shape-1)", {
  set.seed(303L)
  draws <- replicate(20000, sample_inv_gamma_rate(shape = 5, rate = 3))
  expect_equal(mean(draws), 3 / 4, tolerance = 0.1)
})
