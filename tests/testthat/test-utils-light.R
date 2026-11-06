test_that("extract_MAP_draws can extract Z var as factor", {
  expect_equal(object = extract_MAP_draws(draws, "Z"), true_Z)
})

test_that("extract_MAP_draws can extract W var as factor", {
  expect_equal(object = extract_MAP_draws(draws, "W"), true_W)
})

test_that("extract_MAP_draws can extract Z var", {
  expect_equal(
    object = extract_MAP_draws(draws, "Z", factor = FALSE),
    Z_posterior_probs
  )
})

test_that("extract_MAP_draws can extract W var", {
  expect_equal(
    object = extract_MAP_draws(draws, "W", factor = FALSE),
    W_posterior_probs
  )
})

# Computation functions

test_that(".log is safe", {
  normal_x <- c(seq(0.1, 10, by = 0.1))
  smallest_x <- .Machine[["double.xmin"]]
  expect_equal(.log(normal_x), log(normal_x))
  # the last accurate value
  expect_equal(.log(smallest_x), log(smallest_x))
  # And safety guard
  expect_equal(log(0), -Inf)
  expect_equal(.log(0), log(smallest_x))
  # Safety guards can be disabled
  expect_equal(.log(0, eps = NULL), -Inf)
  expect_equal(.log(0, eps = NULL), log(0))
})
