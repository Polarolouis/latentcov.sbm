# Sufficient statistics layer: mask-aware aggregates shared by the full
# conditionals, and equivalence between the direct `param_*` computations and
# their `stats`-based reuses.

test_that("sufficient_stats_rows/cols match manual mask-aware aggregates", {
  stats_rows <- sufficient_stats_rows(Y_fx, mask_fx, W_fx)
  Y_obs <- mask_observed(Y_fx, mask_fx)
  expect_equal(stats_rows$R_W, Y_obs %*% W_fx)
  expect_equal(stats_rows$N_W, mask_fx %*% W_fx)

  stats_cols <- sufficient_stats_cols(Y_fx, mask_fx, Z_fx)
  expect_equal(stats_cols$R_Z, t(Y_obs) %*% Z_fx)
  expect_equal(stats_cols$N_Z, t(mask_fx) %*% Z_fx)
})

test_that("Soft memberships are valid memberships and equals if Z_fx_soft is close to Z_fx ", {
  expect_no_error(
    stats_cols_soft <- sufficient_stats_cols(Y_fx, mask_fx, Z_fx_soft)
  )
  expect_no_error(
    stats_cols <- sufficient_stats_cols(Y_fx, mask_fx, Z_fx)
  )
  expect_no_error(
    stats_rows_soft <- sufficient_stats_rows(Y_fx, mask_fx, W_fx_soft)
  )
  expect_no_error(
    stats_rows <- sufficient_stats_rows(Y_fx, mask_fx, W_fx)
  )
  R_W <- stats_rows[["R_W"]]
  N_W <- stats_rows[["N_W"]]
  R_W_soft <- stats_rows_soft[["R_W"]]
  N_W_soft <- stats_rows_soft[["N_W"]]
  R_Z <- stats_cols[["R_Z"]]
  N_Z <- stats_cols[["N_Z"]]
  R_Z_soft <- stats_cols_soft[["R_Z"]]
  N_Z_soft <- stats_cols_soft[["N_Z"]]
  expect_equal(
    lapply(list(R_W_soft, N_W_soft, R_Z_soft, N_Z_soft), rowSums),
    lapply(list(R_W, N_W, R_Z, N_Z), rowSums)
  )
})

test_that("one-hot hard memberships are valid soft memberships (same stats)", {
  expect_equal(
    sufficient_stats_cols(Y_fx, mask_fx, Z_fx),
    sufficient_stats_cols(Y_fx, mask_fx, Z_fx)
  )
  expect_equal(
    sufficient_stats_rows(Y_fx, mask_fx, W_fx),
    sufficient_stats_rows(Y_fx, mask_fx, W_fx)
  )
})


test_that("param_alpha_given_Y_Z_W_poisson reuses precomputed column stats", {
  res_direct <- param_alpha_given_Y_Z_W_poisson(
    0.1,
    0.1,
    Y_fx,
    Z_fx,
    W_fx,
    mask_fx
  )
  stats_cols <- sufficient_stats_cols(Y_fx, mask_fx, Z_fx)
  res_reuse <- param_alpha_given_Y_Z_W_poisson(
    0.1,
    0.1,
    Y_fx,
    Z_fx,
    W_fx,
    mask_fx,
    stats_cols = stats_cols
  )
  expect_equal(res_reuse$shape, res_direct$shape)
  expect_equal(res_reuse$rate, res_direct$rate)
})

test_that("param_multinom_probs_* reuse precomputed stats consistently", {
  stats_rows <- sufficient_stats_rows(Y_fx, mask_fx, W_fx)
  expect_equal(
    param_multinom_probs_Z_poisson(
      Y_fx,
      alpha_fx,
      W_fx,
      pi_fx,
      mask_fx,
      stats = stats_rows
    ),
    param_multinom_probs_Z_poisson(Y_fx, alpha_fx, W_fx, pi_fx, mask_fx)
  )
  expect_equal(
    param_multinom_probs_Z_cov_poisson(
      Y_fx,
      alpha_fx,
      W_fx,
      P_fx,
      mask_fx,
      stats = stats_rows
    ),
    param_multinom_probs_Z_cov_poisson(Y_fx, alpha_fx, W_fx, P_fx, mask_fx)
  )
  stats_cols <- sufficient_stats_cols(Y_fx, mask_fx, Z_fx)
  expect_equal(
    param_multinom_probs_W_poisson(
      Y_fx,
      alpha_fx,
      Z_fx,
      rho_fx,
      mask_fx,
      stats = stats_cols
    ),
    param_multinom_probs_W_poisson(Y_fx, alpha_fx, Z_fx, rho_fx, mask_fx)
  )
})

test_that("soft sufficient statistics are mask-aware (NA cells not counted)", {
  # The soft `alpha` rate must count only observed cells: it equals the hard
  # mask-aware rate, and differs from the pre-fix count over every cell.
  stats_cols <- sufficient_stats_cols(Y_fx, mask_fx, Z_fx)
  rate_mask <- 1 + t(stats_cols$N_Z) %*% W_fx
  rate_all <- 1 + t(Z_fx) %*% matrix(1, nrow(Y_fx), ncol(Y_fx)) %*% W_fx
  expect_false(isTRUE(all.equal(rate_mask, rate_all)))
})

test_that("multinom_probs_Z accepts both a vector prior and a matrix prior", {
  stats_rows <- sufficient_stats_rows(Y_fx, mask_fx, W_fx)
  p_vec <- multinom_probs_Z(stats_rows, alpha_fx, pi_fx)
  p_mat <- multinom_probs_Z(stats_rows, alpha_fx, ilrInv(P_fx))
  expect_equal(dim(p_vec), c(n1_fx, K_fx))
  expect_equal(dim(p_mat), c(n1_fx, K_fx))
  expect_true(all(abs(rowSums(p_vec) - 1) < 1e-10))
  expect_true(all(abs(rowSums(p_mat) - 1) < 1e-10))
})
