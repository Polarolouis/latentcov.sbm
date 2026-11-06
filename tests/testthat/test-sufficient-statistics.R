# Sufficient statistics layer: mask-aware aggregates shared by the full
# conditionals, and equivalence between the direct `param_*` computations and
# their `stats`-based reuses.
emission_distributions <- c("poisson", "bernoulli")

for (emission in emission_distributions) {
  Y <- switch(emission, "bernoulli" = Y_fx_bern, "poisson" = Y_fx)
  test_that(
    glue::glue(
      "{emission} emission: sufficient_stats_rows/cols match manual mask-aware aggregates"
    ),
    {
      stats_rows <- sufficient_stats_rows(
        Y,
        mask_fx,
        W_fx,
        emission = emission
      )

      Y_obs <- mask_observed(Y, mask_fx)

      expected_stats_rows <- switch(
        emission,
        "bernoulli" = list(
          "R_W" = Y_obs %*% W_fx,
          "F_W" = (mask_fx - Y_obs) %*% W_fx
        ),
        "poisson" = list("R_W" = Y_obs %*% W_fx, "N_W" = mask_fx %*% W_fx)
      )
      lapply(seq_along(expected_stats_rows), function(idx) {
        expect_equal(
          stats_rows[[names(expected_stats_rows)[idx]]],
          expected = expected_stats_rows[[idx]]
        )
      })

      stats_cols <- sufficient_stats_cols(Y, mask_fx, Z_fx, emission = emission)
      expected_stats_cols <- switch(
        emission,
        "bernoulli" = list(
          "R_Z" = t(Y_obs) %*% Z_fx,
          "F_Z" = t(mask_fx - Y_obs) %*% Z_fx
        ),
        "poisson" = list("R_Z" = t(Y_obs) %*% Z_fx, "N_Z" = t(mask_fx) %*% Z_fx)
      )
      lapply(seq_along(expected_stats_cols), function(idx) {
        expect_equal(
          stats_cols[[names(expected_stats_cols)[idx]]],
          expected = expected_stats_cols[[idx]]
        )
      })
    }
  )

  test_that(
    glue::glue(
      "{emission} emission: sufficient_stats bind sufficient_stats_rows and sufficient_stats_cols in one list"
    ),
    {
      stats_rows <- sufficient_stats_rows(Y, mask_fx, W_fx, emission = emission)
      stats_cols <- sufficient_stats_cols(Y, mask_fx, Z_fx, emission = emission)

      expect_equal(
        sufficient_stats(Y, mask_fx, Z_fx, W_fx, emission = emission),
        c(stats_cols, stats_rows)
      )
    }
  )

  test_that(
    glue::glue(
      "{emission} emission: soft memberships are valid memberships and equals if Z_fx_soft is close to Z_fx "
    ),
    {
      expect_no_error(
        stats_cols_soft <- sufficient_stats_cols(
          Y,
          mask_fx,
          Z_fx_soft,
          emission = emission
        )
      )
      expect_no_error(
        stats_cols <- sufficient_stats_cols(
          Y,
          mask_fx,
          Z_fx,
          emission = emission
        )
      )
      expect_no_error(
        stats_rows_soft <- sufficient_stats_rows(
          Y,
          mask_fx,
          W_fx_soft,
          emission = emission
        )
      )
      expect_no_error(
        stats_rows <- sufficient_stats_rows(
          Y,
          mask_fx,
          W_fx,
          emission = emission
        )
      )

      expect_equal(
        lapply(c(stats_rows_soft, stats_cols_soft), rowSums),
        lapply(c(stats_rows, stats_cols), rowSums)
      )
    }
  )

  test_that(
    glue::glue(
      "{emission} emission: one-hot hard memberships are valid soft memberships (same stats)"
    ),
    {
      expect_no_error(
        sufficient_stats_cols(Y, mask_fx, Z_fx, emission = emission)
      )
      expect_no_error(
        sufficient_stats_rows(Y, mask_fx, W_fx, emission = emission)
      )
    }
  )

  test_that(
    glue::glue(
      "{emission} emission: param_alpha_given_Y_Z_W_{emission} reuses precomputed column stats"
    ),
    {
      param_alpha_given_Y_Z_W_emission <- switch(
        emission,
        "bernoulli" = param_alpha_given_Y_Z_W_bernoulli,
        "poisson" = param_alpha_given_Y_Z_W_poisson
      )

      res_direct <- param_alpha_given_Y_Z_W_emission(
        0.1,
        0.1,
        Y,
        Z_fx,
        W_fx,
        mask_fx
      )
      stats_cols <- sufficient_stats_cols(Y, mask_fx, Z_fx, emission = emission)
      res_reuse <- param_alpha_given_Y_Z_W_emission(
        0.1,
        0.1,
        Y,
        Z_fx,
        W_fx,
        mask_fx,
        stats_cols = stats_cols
      )

      lapply(seq_along(res_direct), function(idx) {
        expect_equal(res_direct[[idx]], res_reuse[[idx]])
      })
    }
  )

  test_that(
    glue::glue(
      "{emission} emission: param_multinom_probs_{emission} reuse precomputed stats consistently"
    ),
    {
      param_multinom_probs_Z_emission <- switch(
        emission,
        "bernoulli" = param_multinom_probs_Z_bernoulli,
        "poisson" = param_multinom_probs_Z_poisson
      )
      param_multinom_probs_W_emission <- switch(
        emission,
        "bernoulli" = param_multinom_probs_W_bernoulli,
        "poisson" = param_multinom_probs_W_poisson
      )
      stats_rows <- sufficient_stats_rows(Y, mask_fx, W_fx, emission = emission)
      expect_equal(
        param_multinom_probs_Z_emission(
          Y,
          alpha_fx,
          W_fx,
          pi_fx,
          mask_fx,
          stats = stats_rows
        ),
        param_multinom_probs_Z_emission(Y, alpha_fx, W_fx, pi_fx, mask_fx)
      )
      expect_equal(
        param_multinom_probs_Z_emission(
          Y,
          alpha_fx,
          W_fx,
          ilrInv(P_fx),
          mask_fx,
          stats = stats_rows
        ),
        param_multinom_probs_Z_emission(
          Y,
          alpha_fx,
          W_fx,
          ilrInv(P_fx),
          mask_fx
        )
      )
      stats_cols <- sufficient_stats_cols(Y, mask_fx, Z_fx, emission = emission)
      expect_equal(
        param_multinom_probs_W_emission(
          Y,
          alpha_fx,
          Z_fx,
          rho_fx,
          mask_fx,
          stats = stats_cols
        ),
        param_multinom_probs_W_emission(Y, alpha_fx, Z_fx, rho_fx, mask_fx)
      )
    }
  )

  test_that(
    glue::glue(
      "{emission} emission: soft sufficient statistics are mask-aware (NA cells not counted)"
    ),
    {
      # The soft `alpha` rate must count only observed cells: it equals the hard
      # mask-aware rate, and differs from the pre-fix count over every cell.
      stats_cols <- sufficient_stats_cols(Y, mask_fx, Z_fx, emission = emission)

      Y_obs <- mask_observed(Y, mask_fx)

      switch(
        emission,
        "bernoulli" = {
          masked <- t(stats_cols[["F_Z"]]) %*% W_fx
          not_masked <- t(Z_fx) %*% (1 - Y_obs) %*% W_fx
        },
        "poisson" = {
          masked <- 1 + t(stats_cols[["N_Z"]]) %*% W_fx
          not_masked <- 1 + t(Z_fx) %*% matrix(1, nrow(Y), ncol(Y)) %*% W_fx
        }
      )

      expect_false(isTRUE(all.equal(masked, not_masked)))
    }
  )

  test_that(
    glue::glue(
      "{emission} emission: multinom_probs_Z_{emission} accepts both a vector prior and a matrix prior"
    ),
    {
      multinom_probs_Z_emission <- switch(
        emission,
        "bernoulli" = multinom_probs_Z_bernoulli,
        "poisson" = multinom_probs_Z_poisson,
      )
      stats_rows <- sufficient_stats_rows(Y, mask_fx, W_fx, emission = emission)
      p_vec <- multinom_probs_Z_emission(stats_rows, alpha_fx, pi_fx)
      p_mat <- multinom_probs_Z_emission(stats_rows, alpha_fx, ilrInv(P_fx))
      expect_equal(dim(p_vec), c(n1_fx, K_fx))
      expect_equal(dim(p_mat), c(n1_fx, K_fx))
      expect_true(all(abs(rowSums(p_vec) - 1) < 1e-10))
      expect_true(all(abs(rowSums(p_mat) - 1) < 1e-10))
    }
  )

  test_that(
    glue::glue(
      "{emission} emission: multinom_probs_W_{emission} works correctly"
    ),
    {
      multinom_probs_W_emission <- switch(
        emission,
        "bernoulli" = multinom_probs_W_bernoulli,
        "poisson" = multinom_probs_W_poisson,
      )
      stats_cols <- sufficient_stats_cols(Y, mask_fx, Z_fx, emission = emission)
      p_vec <- multinom_probs_W_emission(stats_cols, alpha_fx, rho_fx)
      expect_equal(dim(p_vec), c(n2_fx, R_fx))
      expect_true(all(abs(rowSums(p_vec) - 1) < 1e-10))
    }
  )

  test_that(
    glue::glue(
      "{emission} emission: alpha_params_from_stats_cols_{emission} works correctly"
    ),
    {
      a0 <- 1
      b0 <- 1
      alpha_params_from_stats_cols_emission <- switch(
        emission,
        "bernoulli" = alpha_params_from_stats_cols_bernoulli,
        "poisson" = alpha_params_from_stats_cols_poisson,
      )
      stats_cols <- sufficient_stats_cols(Y, mask_fx, Z_fx, emission = emission)

      expected_alpha_params <- switch(
        emission,
        "bernoulli" = list(
          shape1 = a0 + t(stats_cols[["R_Z"]]) %*% W_fx,
          shape2 = b0 + t(stats_cols[["F_Z"]]) %*% W_fx
        ),
        "poisson" = list(
          shape = a0 + t(stats_cols[["R_Z"]]) %*% W_fx,
          rate = b0 + t(stats_cols[["N_Z"]]) %*% W_fx
        )
      )

      alpha_params <- alpha_params_from_stats_cols_emission(
        a0 = a0,
        b0 = b0,
        W_fx,
        stats_cols
      )
      expect_equal(alpha_params, expected_alpha_params)
    }
  )
}
