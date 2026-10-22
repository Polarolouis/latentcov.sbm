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
