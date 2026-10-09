# --- entropy ------------------------------------------------------------------

test_that("entropy of a uniform 2-category distribution is log(2)", {
  expect_equal(entropy(c(0.5, 0.5)), log(2), tolerance = 1e-12)
  expect_equal(entropy(c(2, 2)), log(2), tolerance = 1e-12)
})

test_that("entropy in log2/log10 scales rescales log entropy", {
  expect_equal(entropy(c(0.5, 0.5), unit = "log2"), 1, tolerance = 1e-12)
  expect_equal(entropy(c(0.5, 0.5), unit = "log10"), log10(2), tolerance = 1e-12)
})

test_that("entropy is maximal for uniform and zero for degenerate", {
  expect_equal(entropy(c(1)), 0)
  expect_equal(entropy(c(1, 0)), 0)
  expect_equal(entropy(c(0.25, 0.25, 0.25, 0.25)), log(4), tolerance = 1e-12)
  H_unif <- entropy(rep(1 / 4, 4))
  H_skew <- entropy(c(0.7, 0.1, 0.1, 0.1))
  expect_lt(H_skew, H_unif)
})

test_that("entropy is invariant to rescaling (frequencies need not sum to 1)", {
  expect_equal(entropy(c(3, 1)), entropy(c(0.75, 0.25)), tolerance = 1e-12)
})

# --- conditional_entropy and mutual_information --------------------------------

test_that("conditional_entropy and mutual_information satisfy the chain rule", {
  tab <- matrix(c(5, 1, 1, 5), 2, 2)
  H_joint <- entropy(tab)
  H_col <- entropy(colSums(tab))
  H_row <- entropy(rowSums(tab))
  expect_equal(conditional_entropy(tab), H_joint - H_col, tolerance = 1e-12)
  expect_equal(mutual_information(tab), H_row + H_col - H_joint, tolerance = 1e-12)
  expect_true(mutual_information(tab) >= 0)
})

test_that("mutual_information of an independent table is (near) zero", {
  tab_indep <- outer(c(0.6, 0.4), c(0.7, 0.2, 0.1)) * 100
  expect_lt(mutual_information(tab_indep), 1e-10)
})

test_that("entropy warns away invalid unit choices via match.arg", {
  expect_error(entropy(c(0.5, 0.5), unit = "log11"))
})

# --- row_normalize_matrix ------------------------------------------------------

test_that("row_normalize_matrix normalizes log-scale rows to sum to 1", {
  x <- matrix(c(log(2), log(4)), 1, 2)
  out <- row_normalize_matrix(x, is_log = TRUE)
  expect_equal(rowSums(out), 1)
  expect_equal(c(out), c(1 / 3, 2 / 3), tolerance = 1e-12)
})

test_that("row_normalize_matrix takes raw (non-log) input when is_log = FALSE", {
  out <- row_normalize_matrix(matrix(c(2, 4), 1, 2), is_log = FALSE)
  expect_equal(c(out), c(1 / 3, 2 / 3), tolerance = 1e-12)
})

test_that("row_normalize_matrix clamps to [tol, 1 - tol] when tol is given", {
  x <- matrix(c(0, 1000), 1, 2)
  out <- row_normalize_matrix(x, is_log = FALSE, tol = 0.1)
  expect_equal(c(out), c(0.1, 0.9))
})

test_that("row_normalize_matrix returns a matrix for several rows", {
  x <- matrix(c(
    log(1), log(3),
    log(2), log(2)
  ), 2, 2, byrow = TRUE)
  out <- row_normalize_matrix(x, is_log = TRUE)
  expect_equal(dim(out), c(2, 2))
  expect_equal(rowSums(out), rep(1, 2), tolerance = 1e-12)
})
