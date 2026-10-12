# --- matrix-normal density ----------------------------------------------------

test_that("dmatrixnormal_plain and dmatrixnormal_schur coincide", {
  set.seed(11L)
  X <- matrix(rnorm(6 * 3), 6, 3)
  M <- matrix(rnorm(6 * 3), 6, 3)
  U <- crossprod(matrix(rnorm(36), 6)) + diag(6)
  V <- crossprod(matrix(rnorm(9), 3)) + diag(3)
  expect_equal(
    dmatrixnormal_plain(X, M, U, V),
    dmatrixnormal_schur(X, M, U, V),
    tolerance = 1e-10
  )
  expect_equal(
    dmatrixnormal_plain(X, M, U, V, log_p = TRUE),
    dmatrixnormal_schur(X, M, U, V, log_p = TRUE),
    tolerance = 1e-10
  )
})

test_that("the 1x1 matrix normal is a standard normal density", {
  expect_equal(dmatrixnormal_plain(matrix(0, 1, 1), matrix(0, 1, 1), matrix(1, 1, 1), matrix(1, 1, 1)), dnorm(0))
  expect_equal(dmatrixnormal_plain(matrix(1, 1, 1), matrix(0, 1, 1), diag(1), diag(1)), dnorm(1))
})

test_that("log_p toggles between density and log-density", {
  X <- matrix(c(0.2, -0.4, 1.1, 0.5), 2, 2)
  M <- matrix(0, 2, 2)
  U <- diag(2)
  V <- diag(2)
  expect_equal(
    dmatrixnormal_plain(X, M, U, V, log_p = TRUE),
    log(dmatrixnormal_plain(X, M, U, V, log_p = FALSE))
  )
})

test_that("the normalizing constant scales with det(V)^n and det(U)^p", {
  k <- 0.5
  X <- matrix(c(0.1, 0.7, -0.3, 0.2), 2, 2)
  M <- matrix(0, 2, 2)
  d_scale <- dmatrixnormal_plain(X, M, diag(2), diag(2) * k)
  # X = 0-mean, dev = X ; U = I, V = kI
  expect_equal(
    d_scale,
    (1 / k)^(2) * dmatrixnormal_plain(X / sqrt(k), M, diag(2), diag(2)),
    tolerance = 1e-8
  )
})

test_that("dmatrixnormal_cpp matches the R implementation", {
  set.seed(12L)
  X <- matrix(rnorm(4 * 2), 4, 2)
  M <- matrix(rnorm(4 * 2), 4, 2)
  U <- crossprod(matrix(rnorm(16), 4)) + diag(4)
  V <- crossprod(matrix(rnorm(4), 2)) + diag(2)
  expect_equal(
    dmatrixnormal_cpp(X, M, U, V),
    dmatrixnormal_schur(X, M, U, V),
    tolerance = 1e-10
  )
  expect_equal(
    dmatrixnormal_cpp(X, M, U, V, log_p = TRUE),
    dmatrixnormal_schur(X, M, U, V, log_p = TRUE),
    tolerance = 1e-10
  )
})

# --- matrix-normal sampling ---------------------------------------------------

test_that("rmatrixnormal reaches its moments", {
  set.seed(21L)
  U <- diag(6) * 2
  V <- diag(2) * 1.5
  Y <- replicate(2000, rmatrixnormal(matrix(0, 6, 2), U, V))
  # E[Y] = M and Var(Y_ij) = U_ii * V_jj
  expect_equal(mean(Y), 0, tolerance = 0.1)
  expect_equal(var(as.vector(Y)), U[1, 1] * V[1, 1], tolerance = 0.3)
})

test_that("rmatrixnormal_cpp reproduces rmatrixnormal under the same seed", {
  set.seed(22L)
  a <- rmatrixnormal(matrix(0, 3, 2), diag(3), diag(2))
  set.seed(22L)
  b <- rmatrixnormal_cpp(matrix(0, 3, 2), diag(3), diag(2))
  expect_equal(a, b, tolerance = 1e-12)
})

# --- U inverse helper ---------------------------------------------------------

test_that("compute_Uinv returns the matrix inverse", {
  set.seed(31L)
  U <- crossprod(matrix(rnorm(25), 5)) + diag(5)
  expect_equal(compute_Uinv(U), solve(U), tolerance = 1e-10)
})

test_that("compute_Uinv is symmetric for a symmetric input", {
  expect_equal(compute_Uinv(diag(4) * 3 + 1), solve(diag(4) * 3 + 1), tolerance = 1e-10)
})
