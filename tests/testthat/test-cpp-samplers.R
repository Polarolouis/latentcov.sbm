test_that("sampler_classical behave coherently", {
  library(microbenchmark)
  bench <- microbenchmark("C++" = {
    P_cpp <- sample_P_metropolis_classical_cpp(P = Pnonindep, Z = Znonindep, Sigma = Sigma, sigma2 = sigma2star, niter_metropolis = 1, minibatch = FALSE)
  }, "R" = {
    P_R <- sample_P_metropolis_classical(P = Pnonindep, Z = Znonindep, Sigma = Sigma, sigma2 = sigma2star, niter_metropolis = 1, minibatch = FALSE)
  }, times = 5)

  repet_cpp <- lapply(seq(10000), function(i) {
    P_cpp <- sample_P_metropolis_classical_cpp(P = Pnonindep[1:10, , drop = FALSE], Z = Znonindep, Sigma = Sigma[1:10, 1:10, drop = FALSE], sigma2 = sigma2star, niter_metropolis = 1, minibatch = FALSE)
    return(data.frame(iter = i, P = P_cpp[1, 1]))
  }) |> do.call(what = "rbind")
  repet_R <- lapply(seq(10000), function(i) {
    P_R <- sample_P_metropolis_classical(P = Pnonindep[1:10, , drop = FALSE], Z = Znonindep, Sigma = Sigma[1:10, 1:10], sigma2 = sigma2star, niter_metropolis = 1, minibatch = FALSE)
    return(data.frame(iter = i, P = P_R[1, 1]))
  }) |> do.call(what = "rbind")
  repet_R |> summarise_at("P", list("mean" = mean, "sd" = sd))
  repet_cpp |> summarise_at("P", list("mean" = mean, "sd" = sd))
  bench_rng <- microbenchmark(
    "C++" = {
      set.seed(3)
      as.vector(arma_normal_test(5))
    },
    "C++ with R RNG" = {
      set.seed(3)
      as.vector(R_normal_test(5))
    },
    "R" = {
      set.seed(3)
      rnorm(5)
    }, times = 1000
  )
  bench_rng
  expect_equal(object = pivot_coord_inv(x = Pindep), expected = pivotCoordInv(Pindep), tolerance = 1e-6)
  expect_equal(object = pivot_coord_inv(x = Pnonindep), expected = pivotCoordInv(Pnonindep), tolerance = 1e-6)

  expect_equal(object = pivot_coord_inv(x = Pindep, log = TRUE), expected = log(pivotCoordInv(Pindep)), tolerance = 1e-6)
  expect_equal(object = pivot_coord_inv(x = Pnonindep, log = TRUE), expected = log(pivotCoordInv(Pnonindep)), tolerance = 1e-6)
})
