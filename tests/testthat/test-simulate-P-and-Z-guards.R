test_that("simulate_P_and_Z returns the expected components", {
  K <- 3
  n <- 30
  Sigma <- matrix(-0.2, n, n)
  Sigma[outer(rep(1:K, each = n %/% K), rep(1:K, each = n %/% K), `==`)] <- 0.9
  diag(Sigma) <- 1
  res <- simulate_P_and_Z(K = K, Sigma = Sigma, sigma2 = 1)
  expect_type(res, "list")
  expect_true(all(c("Z", "P", "probs") %in% names(res)))
  expect_length(res$Z, n)
  expect_true(all(res$Z %in% seq_len(K)))
  expect_equal(dim(res$P), c(n, K - 1))
  expect_equal(dim(res$probs), c(n, K))
})

test_that("simulate_P_and_Z validates `K`", {
  S <- diag(3)
  expect_error(simulate_P_and_Z(1, S, 1), "K")
  expect_error(simulate_P_and_Z(0L, S, 1), "K")
  expect_error(simulate_P_and_Z(2.5, S, 1), "K")
  expect_error(simulate_P_and_Z("a", S, 1), "K")
  expect_error(simulate_P_and_Z(NA_real_, S, 1), "K")
})

test_that("simulate_P_and_Z validates `sigma2`", {
  S <- diag(3)
  expect_error(simulate_P_and_Z(3, S, 0), "sigma2")
  expect_error(simulate_P_and_Z(3, S, -1), "sigma2")
  expect_error(simulate_P_and_Z(3, S, Inf), "sigma2")
  expect_error(simulate_P_and_Z(3, S, "1"), "sigma2")
})

test_that("simulate_P_and_Z validates the shape of `Sigma`", {
  expect_error(simulate_P_and_Z(3, matrix(1, 3, 4), 1), "square")
  expect_error(simulate_P_and_Z(3, as.numeric(diag(3)), 1), "numeric matrix")
  Sns <- matrix(1, 4, 4)
  Sns[1, 2] <- 0.5
  Sns[2, 1] <- -0.5
  expect_error(simulate_P_and_Z(3, Sns, 1), "symmetric")
  expect_error(simulate_P_and_Z(3, diag(2), 1), "at least `K`")
})

test_that("simulate_P_and_Z rejects matrices that are not positive semi-definite", {
  S <- matrix(-1, 6, 6)
  diag(S) <- 1
  expect_error(simulate_P_and_Z(3, S, 1), "positive semi-definite")
})

test_that("simulate_P_and_Z warns and still samples on (numerically) singular Sigma", {
  S <- matrix(1, 6, 6)
  expect_warning(res <- simulate_P_and_Z(3, S, 1), "singular")
  expect_equal(dim(res$P), c(6, 2))
})

test_that("simulate_P_and_Z warns when Sigma is not a correlation matrix", {
  S <- diag(3) * 2
  expect_warning(simulate_P_and_Z(3, S, 1), "correlation matrix")
})
