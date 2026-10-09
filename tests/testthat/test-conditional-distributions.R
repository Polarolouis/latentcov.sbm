# --- conditional latent position (Pi) -----------------------------------------

test_that("cond_Pi_fast and cond_Pi_given_P_min_i_sigma agree for the same Sigma", {
  set.seed(501L)
  P <- P_fx
  Sigma <- Sigma_fx
  Theta <- solve(Sigma)
  for (i in 1:4) {
    got_fast <- cond_Pi_fast(P, Theta, sigma2 = 0.8, i = i)
    got_full <- cond_Pi_given_P_min_i_sigma(P, Sigma, sigma2 = 0.8, i = i)
    expect_equal(as.vector(got_fast$mean), as.vector(got_full$mean), tolerance = 1e-10)
    expect_equal(as.vector(got_fast$cov), as.vector(got_full$cov), tolerance = 1e-10)
  }
})

test_that("cond_Pi_fast with diagonal Theta recovers an iid normal mean", {
  P <- P_fx
  Theta <- diag(nrow(P)) * 3
  res <- cond_Pi_fast(P, Theta, sigma2 = 1.5, i = 2)
  expect_equal(res$mean, rep(0, ncol(P)))
  expect_equal(res$cov, (1.5 / 3) * diag(ncol(P)))
})

test_that("sample_Pi_given draws from the conditional with the right first moment", {
  set.seed(502L)
  draws <- replicate(20000, sample_Pi_given(P_fx, Sigma_fx, sigma2 = 0.5, i = 1))
  theor <- cond_Pi_given_P_min_i_sigma(P_fx, Sigma_fx, sigma2 = 0.5, i = 1)
  expect_equal(apply(draws, 1, mean), as.vector(theor$mean), tolerance = 0.05)
  expect_equal(nrow(draws), ncol(P_fx))
  expect_false(anyNA(draws))
})

# --- pi | Z -------------------------------------------------------------------

test_that("param_pi_given_Z adds the observed counts to the prior", {
  expect_equal(param_pi_given_Z(etas = c(1, 1, 1), Z = Z_fx), c(3, 3, 3))
  expect_equal(param_pi_given_Z(etas = c(2, 0.5, 1), Z = Z_fx), c(4, 2.5, 3))
})

test_that("sample_pi_given_Z lies on the simplex and matches the Dirichlet mean", {
  set.seed(503L)
  etas_post <- c(5, 3, 2)
  draws <- replicate(20000, sample_pi_given_Z(etas_post))
  expect_equal(apply(draws, 1, mean), etas_post / sum(etas_post), tolerance = 0.05)
  expect_equal(colSums(draws), rep(1, 20000))
  expect_true(all(draws > 0))
})

# --- rho | W ------------------------------------------------------------------

test_that("param_rho_given_W adds the observed column counts to the prior", {
  expect_equal(param_rho_given_W(gammas = c(1, 1), W = W_fx), c(4, 3))
})

test_that("sample_rho_given_W is a positive simplex vector", {
  set.seed(504L)
  draws <- replicate(2000, sample_rho_given_W(c(3, 1)))
  expect_true(all(draws > 0))
  expect_equal(colSums(draws), rep(1, 2000))
  expect_equal(mean(draws[1, ]), 0.75, tolerance = 0.05)
})

# --- multinomial probabilities of Z / W ---------------------------------------

test_that("param_multinom_probs_Z_poisson gives normalized valid probabilities", {
  probs <- param_multinom_probs_Z_poisson(Y_fx, alpha_fx, W_fx, pi_fx, mask_fx)
  expect_equal(dim(probs), c(n1_fx, K_fx))
  expect_equal(rowSums(probs), rep(1, n1_fx), tolerance = 1e-10)
  expect_true(all(probs > 0 & probs < 1))
})

test_that("param_multinom_probs_Z_poisson is NA-proof through the mask", {
  clean <- mask_NA(Y_fx)$YnoNA
  with_na <- param_multinom_probs_Z_poisson(Y_fx, alpha_fx, W_fx, pi_fx, mask_fx)
  without_na <- param_multinom_probs_Z_poisson(clean, alpha_fx, W_fx, pi_fx, mask_fx)
  expect_equal(with_na, without_na, tolerance = 1e-12)
  expect_false(anyNA(with_na))
})

test_that("param_multinom_probs_Z_cov_poisson uses per-row pi from ilrInv(P)", {
  set.seed(505L)
  P_rep <- matrix(rep(ilr(matrix(pi_fx, 1, K_fx)), n1_fx), n1_fx, K_fx - 1, byrow = TRUE)
  cov_probs <- param_multinom_probs_Z_cov_poisson(Y_fx, alpha_fx, W_fx, P_rep, mask_fx)
  plain_probs <- param_multinom_probs_Z_poisson(Y_fx, alpha_fx, W_fx, pi_fx, mask_fx)
  expect_equal(cov_probs, plain_probs, tolerance = 1e-10)
  expect_equal(rowSums(cov_probs), rep(1, n1_fx), tolerance = 1e-10)
})

test_that("param_multinom_probs_W_poisson is normalized and NA-proof", {
  probs <- param_multinom_probs_W_poisson(Y_fx, alpha_fx, Z_fx, rho_fx, mask_fx)
  expect_equal(dim(probs), c(n2_fx, R_fx))
  expect_equal(rowSums(probs), rep(1, n2_fx), tolerance = 1e-10)
  clean <- mask_NA(Y_fx)$YnoNA
  expect_equal(
    param_multinom_probs_W_poisson(clean, alpha_fx, Z_fx, rho_fx, mask_fx),
    probs,
    tolerance = 1e-12
  )
})

test_that("tol clamps extreme membership probabilities", {
  extreme_alpha <- matrix(c(100, 1, 1, 1, 1, 100), K_fx, R_fx)
  probs <- param_multinom_probs_Z_poisson(Y_fx, extreme_alpha, W_fx, pi_fx, mask_fx, tol = 0.05)
  expect_true(all(probs >= 0.05 - 1e-10))
  expect_true(all(probs <= 0.95 + 1e-10))
})

# --- sampling memberships from normalized probs --------------------------------

test_that("sample_Z/W return labels within range and respect degenerate probs", {
  set.seed(506L)
  ones <- .one_hot(c(2, 1, 2, 3, 1, 3), K_fx)
  expect_equal(sample_Z_given_alpha_pi_Y_W(ones), c(2, 1, 2, 3, 1, 3))
  expect_equal(sample_Z_given_alpha_P_Y_W(ones), c(2, 1, 2, 3, 1, 3))
  Wo <- .one_hot(c(1, 2, 1), R_fx)
  expect_equal(sample_W_given_alpha_rho_Y_Z(Wo), c(1, 2, 1))
})

test_that("sample_Z draws from a soft distribution with the right frequency", {
  set.seed(507L)
  probs <- matrix(c(0.8, 0.2, 0.2, 0.8), 2, 2, byrow = TRUE)
  draws <- replicate(2000, sample_Z_given_alpha_pi_Y_W(probs))
  expect_equal(mean(draws[1, ]), 1.2, tolerance = 0.05)
  expect_true(all(draws %in% 1:2))
})

# --- alpha | Y, Z, W ----------------------------------------------------------

test_that("param_alpha_given_Y_Z_W_poisson sums counts per block pair", {
  res <- param_alpha_given_Y_Z_W_poisson(a0 = 0.1, b0 = 0.1, Y = Y_fx, Z = Z_fx, W = W_fx, mask = mask_fx)
  expect_equal(dim(res$shape), c(K_fx, R_fx))
  expect_equal(
    res$shape,
    0.1 + t(Z_fx) %*% mask_observed(Y_fx, mask_fx) %*% W_fx
  )
  expect_equal(res$rate, 0.1 + t(Z_fx) %*% mask_fx %*% W_fx)
})

test_that("NAs do not leak into the alpha sufficient statistics", {
  clean <- mask_NA(Y_fx)$YnoNA
  with_na <- param_alpha_given_Y_Z_W_poisson(0, 0, Y_fx, Z_fx, W_fx, mask_fx)
  without_na <- param_alpha_given_Y_Z_W_poisson(0, 0, clean, Z_fx, W_fx, mask_fx)
  expect_equal(with_na, without_na)
  expect_false(anyNA(with_na$shape))
  expect_false(anyNA(with_na$rate))
})

test_that("with a complete mask, rate equals prior + row/col block size products", {
  full_mask <- matrix(1, n1_fx, n2_fx)
  res <- param_alpha_given_Y_Z_W_poisson(0.5, 2, Y_fx, Z_fx, W_fx, full_mask)
  expect_equal(res$rate, 2 + outer(colSums(Z_fx), colSums(W_fx), `*`))
})

test_that("sample_alpha_given_Y_Z_W_poisson matches the Gamma moments", {
  set.seed(508L)
  shape <- matrix(c(5, 1, 2, 3), 2, 2)
  rate <- matrix(c(2, 1, 1, 6), 2, 2)
  draws <- replicate(5000, sample_alpha_given_Y_Z_W_poisson(shape, rate))
  expect_equal(dim(draws), c(2, 2, 5000))
  expect_true(all(draws > 0))
  expect_equal(mean(draws[1, 1, ]), shape[1, 1] / rate[1, 1], tolerance = 0.15)
  expect_equal(mean(draws[2, 2, ]), shape[2, 2] / rate[2, 2], tolerance = 0.15)
})
