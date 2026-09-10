test_that("ilrInvcpp behave the same in C++ and R for orthonormal", {
  expect_equal(object = ilrInvcpp(Pindep), expected = ilrInv(Pindep), tolerance = 1e-6)
  expect_equal(object = ilrInvcpp(Pnonindep), expected = ilrInv(Pnonindep), tolerance = 1e-6)

  expect_equal(object = ilrInvcpp(Pindep, log = TRUE), expected = log(ilrInv(Pindep)), tolerance = 1e-6)
  expect_equal(object = ilrInvcpp(Pnonindep, log = TRUE), expected = log(ilrInv(Pnonindep)), tolerance = 1e-6)
})

test_that("conditional Pi_given_PminI_sigma2 behave the same in C++ and R", {
  res_list <- lapply(seq(n1), function(i) {
    res <- cond_Pi_given_P_min_i_sigma(P = Pnonindep, Sigma = Sigma, sigma2 = sigma2star, i)
  })

  mean_list <- lapply(res_list, "[[", 1)
  cov_list <- lapply(res_list, "[[", 2)

  mean_cpp_list <- lapply(seq(n1), function(i) {
    mean_of_Pi_given_P_min_i_sigma(P = Pnonindep, Sigma = Sigma, sigma2 = sigma2star, i - 1)
  })

  cov_cpp_list <- lapply(seq(n1), function(i) {
    cov_of_Pi_given_P_min_i_sigma(P = Pnonindep, Sigma = Sigma, sigma2 = sigma2star, i - 1)
  })

  expect_equal(mean_cpp_list, mean_list, tolerance = 1e-6)
  expect_equal(cov_cpp_list, cov_list, tolerance = 1e-6)
})
