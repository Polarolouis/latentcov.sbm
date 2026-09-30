test_that("Y * mask does not leak NA into alpha parameters (NA * 0 = NA)", {
  Y <- matrix(c(1, NA, 3, NA, 2, 4), 2, 3, byrow = TRUE)
  mask <- (!is.na(Y)) * 1L
  Z <- matrix(c(1, 0, 0, 1), 2, 2)
  W <- matrix(c(1, 0, 0, 1, 0, 0), 3, 2)

  prm <- param_alpha_given_Y_Z_W_poisson(1, 1, Y, Z, W, mask)
  expect_false(anyNA(prm$shape))
  expect_false(anyNA(prm$rate))
  expect_true(all(prm$shape > 0))
  expect_true(all(prm$rate > 0))
})

test_that("multinomial membership probabilities have no NA when Y contains NAs", {
  Y <- matrix(c(1, 2, 3, NA, 2, 4, 5, NA), 2, 4, byrow = TRUE)
  mask <- (!is.na(Y)) * 1L
  K <- 3
  R <- 2
  alpha <- matrix(c(3, 1, 1, 3, 1, 2), nrow = K, ncol = R, byrow = TRUE)
  W <- matrix(c(1, 0, 0, 1), 4, R)
  Z <- matrix(c(1, 0, 0, 0, 1, 0), nrow(Y), K, byrow = TRUE)
  rho <- c(0.5, 0.5)
  pi <- rep(1 / K, K)
  P <- matrix(c(0.3, -0.4, 0.2, 0.1), nrow(Y), K - 1)

  probsZ <- param_multinom_probs_Z_poisson(Y, alpha, W, pi, mask)
  probsW <- param_multinom_probs_W_poisson(Y, alpha, Z, rho, mask)
  probsZc <- param_multinom_probs_Z_cov_poisson(Y, alpha, W, P, mask)
  expect_false(anyNA(probsZ))
  expect_false(anyNA(probsW))
  expect_false(anyNA(probsZc))
})

test_that("sample functions can draw from NA-free probabilities derived from masked Y", {
  Y <- matrix(c(3, NA, 5, NA, 2, 4), 2, 3, byrow = TRUE)
  mask <- (!is.na(Y)) * 1L
  alpha <- matrix(2, 2, 2)
  W <- matrix(c(1, 0, 0, 1, 1, 0), 3, 2, byrow = TRUE)
  pi <- c(0.6, 0.4)
  probs <- param_multinom_probs_Z_poisson(Y, alpha, W, pi, mask)
  expect_false(anyNA(probs))
  set.seed(1)
  Zs <- sample_Z_given_alpha_P_Y_W(probs)
  expect_length(Zs, nrow(probs))
  expect_true(all(Zs %in% seq_len(2)))
})

test_that("gibbs_sampling_lbm_poisson runs on Y with NAs without producing NA alpha", {
  n1 <- 40
  n2 <- 30
  K <- 3
  R <- 4
  set.seed(1234)
  rho <- c(0.4, 0.3, 0.2, 0.1)
  W <- sample.int(R, n2, replace = TRUE, prob = rho)
  mem <- rep(seq(K), each = floor(n1 / K))
  mem <- c(mem, rep(K, n1 - length(mem)))
  Sigma <- matrix(0.05, n1, n1)
  Sigma[outer(mem, mem, `==`)] <- 0.6
  diag(Sigma) <- 1
  P <- simulate_P(K, Sigma, sigma2 = 1)
  Z <- simulate_Z_from_P(P)
  alpha_sharp <- matrix(c(120, 100, 80, 40, 90, 70, 45, 20, 30, 20, 10, 6), K, byrow = TRUE)
  Y_full <- simulate_SBM_alpha_Z_W(alpha_sharp, Z, W)
  Na_idx <- sample.int(n1 * n2, floor(0.4 * n1 * n2))
  Y <- Y_full
  Y[Na_idx] <- NA

  set.seed(1)
  draws <- gibbs_sampling_lbm_poisson(Y, K, R, niter = 20)
  alpha <- posterior::as_draws_df(posterior::subset_draws(draws, variable = "alpha"))
  expect_false(any(is.na(alpha)))
})
