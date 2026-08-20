# Common
seed <- 1234L
set.seed(seed)
# Rows
K <- 3
sigma2star <- 1
npc <- 30
n1 <- npc * K
phylo_memberships <- rep(seq(K), each = npc)
rho_within <- 0.9
rho_between <- -1 / (K * npc - 1)

# indep
Pindep <- matrix(rnorm(n1 * (K - 1)), nrow = n1)

# non indep
library(mvtnorm)

Sigma <- matrix(rho_between, nrow = n1, ncol = n1)
same_cluster <- outer(phylo_memberships, phylo_memberships, `==`)
Sigma[same_cluster] <- rho_within
diag(Sigma) <- 1


Pnonindep <- simulate_P_and_Z(K = K, Sigma, sigma2star, 1)[["P"]]


test_that("pivot_coord_inv behave the same in C++ and R", {
  expect_equal(object = pivot_coord_inv(x = Pindep), expected = pivotCoordInv(Pindep), tolerance = 1e-6)
  expect_equal(object = pivot_coord_inv(x = Pnonindep), expected = pivotCoordInv(Pnonindep), tolerance = 1e-6)
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
