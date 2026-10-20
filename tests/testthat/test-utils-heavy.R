test_that("stan_flatten flattens trailing parameter dimensions into bracketed names", {
  x <- array(seq_len(3 * 2 * 2), dim = c(3, 2, 2))
  out <- stan_flatten(x, "beta")
  expect_equal(dim(out), c(3, 4))
  expect_equal(dimnames(out)[[1]], c("1", "2", "3"))
  expect_equal(
    dimnames(out)[[2]],
    c("beta[1,1]", "beta[2,1]", "beta[1,2]", "beta[2,2]")
  )
  expect_equal(out[1, "beta[2,1]"], x[1, 2, 1])

  one <- array(1:3, dim = c(3, 1))
  out1 <- stan_flatten(one, "p")
  expect_equal(dimnames(out1)[[2]], "p")
})

test_that("list_arrays_to_stan binds flattened arrays and strips the _array suffix", {
  niter <- 3
  K <- 3
  R <- 2
  n <- 6
  a_alpha <- array(0, dim = c(niter, K, R))
  a_P <- array(0, dim = c(niter, n, K - 1))
  for (it in 1:niter) {
    a_alpha[it, , ] <- matrix(it * 1:6, 3, 2)
    a_P[it, , ] <- matrix(rep(it, length = n * (K - 1)), n, K - 1)
  }
  out <- list_arrays_to_stan(list(alpha_array = a_alpha, P_array = a_P))
  expect_equal(dim(out), c(niter, K * R + n * (K - 1)))
  expect_equal(
    dimnames(out)[[2]],
    c(
      "alpha[1,1]",
      "alpha[2,1]",
      "alpha[3,1]",
      "alpha[1,2]",
      "alpha[2,2]",
      "alpha[3,2]",
      "P[1,1]",
      "P[2,1]",
      "P[3,1]",
      "P[4,1]",
      "P[5,1]",
      "P[6,1]",
      "P[1,2]",
      "P[2,2]",
      "P[3,2]",
      "P[4,2]",
      "P[5,2]",
      "P[6,2]"
    )
  )
  expect_equal(out[1, "alpha[2,1]"], 2)
  expect_equal(out[1, "P[6,2]"], 1)
})

test_that("list_arrays_to_stan skips all-NA arrays", {
  a <- array(0, dim = c(3, 2))
  out <- list_arrays_to_stan(list(
    ok_array = a,
    allna_array = array(NA_real_, c(3, 2))
  ))
  expect_equal(dimnames(out)[[2]], c("ok[1]", "ok[2]"))
})

test_that("list_arrays_to_stan tolerates a list containing only NA arrays", {
  expect_no_error({
    out <- list_arrays_to_stan(list(zz_array = array(NA_real_, c(3, 2))))
    expect_equal(dim(out), c(3, 0))
  })
})

test_that("list_stan_to_chains_stan stacks chains along a labelled third dimension", {
  m1 <- matrix(
    1:6,
    nrow = 3,
    dimnames = list(Iteration = 1:3, Parameter = letters[1:2])
  )
  m2 <- matrix(
    7:12,
    nrow = 3,
    dimnames = list(Iteration = 1:3, Parameter = letters[1:2])
  )
  out <- list_stan_to_chains_stan(list(m1, m2))
  expect_equal(dim(out), c(3, 2, 2))
  expect_equal(dimnames(out)[[3]], c("1", "2"))
  expect_equal(out[,, 1], m1)
  expect_equal(out[,, 2], m2)
})

# ---- label permutations ------------------------------------------------------

test_that("find_permutation_alphas aligns a permuted connectivity matrix", {
  alpha_ref <- matrix(c(1, 3, 5, 2, 4, 6), 3, 2)
  alpha_perm <- alpha_ref[c(2, 3, 1), c(2, 1)]
  perm <- find_permutation_alphas(alpha_ref, alpha_perm)
  expect_equal(perm$row_perm |> sort(), 1:3)
  expect_equal(perm$col_perm |> sort(), 1:2)
  expect_equal(alpha_perm[perm$row_perm, perm$col_perm], alpha_ref)
})

test_that("find_permutation_alphas returns identity for equal matrices", {
  alpha_ref <- matrix(c(1, 3, 5, 2, 4, 6), 3, 2)
  perm <- find_permutation_alphas(alpha_ref, alpha_ref)
  expect_equal(perm$row_perm, 1:3)
  expect_equal(perm$col_perm, 1:2)
})

test_that("find_permutation_alphas_L2 finds the exact minimal permutation", {
  alpha_ref <- matrix(c(1, 3, 5, 2, 4, 6), 3, 2)
  alpha_perm <- alpha_ref[c(2, 3, 1), c(2, 1)]
  perm <- find_permutation_alphas_L2(alpha_ref, alpha_perm)
  expect_equal(alpha_perm[perm$row_perm, perm$col_perm], alpha_ref)
  expect_equal(perm$row_perm, c(3, 1, 2))
  expect_equal(perm$col_perm, c(2, 1))
})

test_that("find_permutation_alphas_L2 returns identity for equal matrices", {
  alpha_ref <- matrix(c(1, 3, 5, 2, 4, 6), 3, 2)
  perm <- find_permutation_alphas_L2(alpha_ref, alpha_ref)
  expect_equal(perm$row_perm, 1:3)
  expect_equal(perm$col_perm, 1:2)
})

test_that("clue and exhaustive permutations agree on the fixture", {
  alpha_ref <- matrix(c(1, 3, 5, 2, 4, 6), 3, 2)
  alpha_perm <- alpha_ref[c(2, 3, 1), c(2, 1)]
  p_clue <- find_permutation_alphas(alpha_ref, alpha_perm)
  p_l2 <- find_permutation_alphas_L2(alpha_ref, alpha_perm)
  expect_equal(
    alpha_perm[p_clue$row_perm, p_clue$col_perm],
    alpha_perm[p_l2$row_perm, p_l2$col_perm]
  )
})

# ---- identifiability ---------------------------------------------------------

test_that("check_lbm_identifiability accepts an identifiable configuration", {
  net <- matrix(0, nrow = 6, ncol = 5)
  alpha_id <- matrix(c(1, 2, 3, 4, 5, 6), 3, 2)
  pi_id <- c(0.5, 0.3, 0.2)
  rho_id <- c(0.6, 0.4)
  expect_message(
    check_lbm_identifiability(net, alpha_id, pi_id, rho_id, K = 3, R = 2),
    "identifiable"
  )
})

test_that("check_lbm_identifiability aborts when alpha%*%rho has duplicates", {
  net <- matrix(0, nrow = 6, ncol = 5)
  alpha <- matrix(c(1, 1, 3, 2, 2, 6), 3, 2)
  expect_error(
    check_lbm_identifiability(
      net,
      alpha,
      c(0.5, 0.3, 0.2),
      c(0.6, 0.4),
      K = 3,
      R = 2
    ),
    "alpha%\\*%rho"
  )
})

test_that("check_lbm_identifiability aborts when t(pi)%*%alpha has duplicates", {
  net <- matrix(0, nrow = 6, ncol = 5)
  alpha <- matrix(c(1, 2, 3, 1, 2, 3), 3, 2)
  expect_error(
    check_lbm_identifiability(
      net,
      alpha,
      c(0.5, 0.3, 0.2),
      c(0.6, 0.4),
      K = 3,
      R = 2
    ),
    "t\\(pi\\)%\\*%alpha"
  )
})

test_that("check_lbm_identifiability aborts with too few row nodes", {
  net <- matrix(0, nrow = 2, ncol = 5)
  alpha <- matrix(c(1, 2, 3, 4, 5, 6), 3, 2)
  expect_error(
    check_lbm_identifiability(
      net,
      alpha,
      c(0.5, 0.3, 0.2),
      c(0.6, 0.4),
      K = 3,
      R = 2
    ),
    "row nodes"
  )
})

test_that("check_lbm_identifiability aborts with too few column nodes", {
  net <- matrix(0, nrow = 6, ncol = 4)
  alpha <- matrix(c(1, 2, 3, 4, 5, 6), 3, 2)
  expect_error(
    check_lbm_identifiability(
      net,
      alpha,
      c(0.5, 0.3, 0.2),
      c(0.6, 0.4),
      K = 3,
      R = 2
    ),
    "column nodes"
  )
})

# ---- ARI table ---------------------------------------------------------------

test_that("ARI_table returns 1 for a perfectly aligned partition", {
  niter <- 3
  nch <- 2
  n <- 6
  Z_true <- rep(1:3, each = 2)
  var_order <- c(
    paste0("P[", rep(1:n, each = 2), ",", rep(1:2, n), "]"),
    paste0("Z[", 1:n, "]")
  )
  arr <- array(
    seq_len(niter * nch * length(var_order)),
    dim = c(niter, nch, length(var_order)),
    dimnames = list(NULL, NULL, var_order)
  )
  arr[,, grep("^Z\\[", var_order)] <- rep(Z_true, each = niter * nch)
  dr <- posterior::as_draws_array(arr)

  ari <- ARI_table(dr, Z_true, label = "Z")
  expect_equal(dim(ari), c(niter, nch))
  expect_equal(unique(c(ari)), 1)
})

test_that("ARI_table accepts a one-hot true membership", {
  niter <- 2
  nch <- 2
  n <- 6
  Z_true <- rep(1:3, each = 2)
  arr <- array(
    0,
    dim = c(niter, nch, n),
    dimnames = list(NULL, NULL, paste0("Z[", 1:n, "]"))
  )
  arr[,,] <- rep(Z_true, each = niter * nch)
  dr <- posterior::as_draws_array(arr)
  ari <- ARI_table(dr, onehot_encode(Z_true), label = "Z")
  expect_equal(unique(c(ari)), 1)
})

test_that("ARI_table is strictly below 1 for a scrambled partition", {
  niter <- 2
  nch <- 1
  n <- 6
  Z_true <- rep(1:3, each = 2)
  Z_scram <- c(1L, 2L, 1L, 2L, 1L, 2L)
  arr <- array(
    0,
    dim = c(niter, nch, n),
    dimnames = list(NULL, NULL, paste0("Z[", 1:n, "]"))
  )
  arr[,,] <- rep(Z_scram, each = niter * nch)
  dr <- posterior::as_draws_array(arr)
  ari <- ARI_table(dr, Z_true, label = "Z")
  expect_true(all(ari < 1))
})

# ---- compute_meanpi ----------------------------------------------------------

test_that("compute_meanpi appends the mean row proportions", {
  niter <- 3
  nch <- 2
  n <- 6
  K <- 3
  Z1 <- rep(1:3, each = 2)
  var_order <- paste0("Z[", 1:n, "]")
  arr <- array(
    0,
    dim = c(niter, nch, n),
    dimnames = list(NULL, NULL, var_order)
  )
  arr[,,] <- rep(Z1, each = niter * nch)
  dr <- posterior::as_draws_array(arr)

  out <- compute_meanpi(dr, K = 3)
  v <- posterior::variables(out)
  expect_true(all(c("meanpi[1]", "meanpi[2]", "meanpi[3]") %in% v))
  df <- as.data.frame(posterior::as_draws_df(out))
  expect_equal(
    unname(unlist(df[1, c("meanpi[1]", "meanpi[2]", "meanpi[3]")])),
    c(1 / 3, 1 / 3, 1 / 3)
  )
  expect_equal(unname(unlist(df[1, paste0("Z[", 1:n, "]")])), Z1)
})

test_that("compute_meanpi infers the number of groups when K is NULL using P", {
  niter <- 2
  nch <- 1
  K <- 3
  n <- 6
  Z1 <- rep(1:3, each = 2)
  arr <- array(
    0,
    dim = c(niter, nch, n * K),
    dimnames = list(
      NULL,
      NULL,
      c(
        paste0("Z[", seq(n), "]"),
        paste0(
          "P[",
          outer(
            seq(n),
            seq(K - 1),
            FUN = Vectorize(function(x, y) paste0(x, ",", y))
          ),
          "]"
        )
      )
    )
  )
  arr[,, seq(n)] <- rep(Z1, each = niter * nch)
  dr <- posterior::as_draws_array(arr)
  out <- compute_meanpi(dr)
  expect_equal(
    posterior::variables(out)[grep("^meanpi", posterior::variables(out))],
    c("meanpi[1]", "meanpi[2]", "meanpi[3]")
  )
})

test_that("compute_meanpi keeps pre-existing variables", {
  niter <- 2
  nch <- 1
  n <- 6
  Z1 <- rep(1:3, each = 2)
  var_order <- c("sigma2", paste0("Z[", 1:n, "]"))
  arr <- array(
    0,
    dim = c(niter, nch, n + 1),
    dimnames = list(NULL, NULL, var_order)
  )
  arr[,, 1] <- 0.5
  arr[,, -1] <- rep(Z1, each = niter * nch)
  dr <- posterior::as_draws_array(arr)
  out <- compute_meanpi(dr, K = 3)
  expect_true("sigma2" %in% posterior::variables(out))
})

# ---- label switching correction -----------------------------------------------

build_label_switched_draws <- function(niter = 3, nch = 2, vars = "none") {
  K <- 3
  R <- 2
  n <- 6
  alpha1 <- matrix(c(1, 3, 5, 2, 4, 6), 3, 2)
  ra <- c(2L, 3L, 1L)
  ca <- c(2L, 1L)
  alpha2 <- alpha1[ra, ca]
  rho1 <- c(0.6, 0.4)
  rho2 <- rho1[ca]
  P1 <- matrix(seq_len(n * (K - 1)) * 0.1, n, K - 1)
  Psi <- default_Psi_function(K)
  Ta <- round(Psi %*% perm_matrix_from_order(ra) %*% t(Psi), 10)
  P2 <- P1 %*% t(Ta)
  Z1 <- rep(1:3, each = 2)
  Z2 <- ra[Z1]
  W1 <- c(1L, 2L, 1L, 2L, 1L)
  W2 <- ca[W1]

  extra <- character()
  extra1 <- extra2 <- NULL
  if (vars == "meanpi") {
    meanpi1 <- c(0.3, 0.4, 0.3)
    extra <- paste0("meanpi[", 1:K, "]")
    extra1 <- meanpi1
    extra2 <- meanpi1[ra]
  } else if (vars == "pi") {
    pi1 <- c(0.3, 0.4, 0.3)
    extra <- paste0("pi[", 1:K, "]")
    extra1 <- pi1
    extra2 <- pi1[ra]
  }

  var_order <- c(
    paste0("alpha[", rep(1:K, R), ",", rep(1:R, each = K), "]"),
    "rho[1]",
    "rho[2]",
    paste0("P[", rep(1:n, K - 1), ",", rep(1:(K - 1), each = n), "]"),
    extra,
    paste0("Z[", 1:n, "]"),
    paste0("W[", 1:5, "]")
  )
  arr <- array(
    0,
    dim = c(niter, nch, length(var_order)),
    dimnames = list(NULL, NULL, var_order)
  )
  for (it in 1:niter) {
    for (ch in 1:nch) {
      arr[it, ch, ] <- if (ch == 1) {
        c(alpha1, rho1, c(P1), extra1, Z1, W1)
      } else {
        c(alpha2, rho2, c(P2), extra2, Z2, W2)
      }
    }
  }
  posterior::as_draws_array(arr)
}

chains_sampled_equal <- function(draws) {
  df <- as.data.frame(posterior::as_draws_df(draws))
  sel <- !grepl("^\\.", names(df))
  c1 <- as.matrix(df[df$.chain == 1, sel])
  c2 <- as.matrix(df[df$.chain == 2, sel])
  expect_equal(dim(c1), dim(c2))
  expect_lt(max(abs(c1 - c2)), 1e-9)
}

test_that("delabel_switch_stan aligns a permuted second chain onto the first", {
  ddr <- build_label_switched_draws(vars = "meanpi")
  out <- delabel_switch_stan(ddr, K = 3, R = 2)
  chains_sampled_equal(out)
  df <- as.data.frame(posterior::as_draws_df(out))
  expect_equal(unique(df[["meanpi[1]"]]), 0.3)
  expect_equal(unique(df[["meanpi[2]"]]), 0.4)
  expect_equal(unique(df[["meanpi[3]"]]), 0.3)
})

test_that("delabel_switch_stan honours an explicit alpha_ref", {
  dr <- build_label_switched_draws()
  alpha_ref <- matrix(c(1, 3, 5, 2, 4, 6), 3, 2)
  out <- delabel_switch_stan(dr, K = 3, R = 2, alpha_ref = alpha_ref)
  df <- as.data.frame(posterior::as_draws_df(out))
  sel <- !grepl("^\\.", names(df))
  expect_equal(
    unname(as.matrix(df[df$.chain == 1, sel])),
    unname(as.matrix(df[df$.chain == 2, sel])),
    tolerance = 1e-9
  )
})

test_that("delabel_switch_stan requires at least two chains or an alpha_ref", {
  dr <- build_label_switched_draws(nch = 1)
  expect_error(delabel_switch_stan(dr, K = 3, R = 2), "at least two chains")
})

test_that("delabel_switch_stan flattens chain labels onto the second permutation", {
  dr <- build_label_switched_draws(niter = 3, nch = 3)
  out <- delabel_switch_stan(dr, K = 3, R = 2)
  df <- as.data.frame(posterior::as_draws_df(out))
  sel <- !grepl("^\\.", names(df))
  c1 <- as.matrix(df[df$.chain == 1, sel])
  c3 <- as.matrix(df[df$.chain == 3, sel])
  expect_lt(max(abs(c1 - c3)), 1e-9)
})

# ---- lbm_results_to_stan_draws -----------------------------------------------

mk_chain_results <- function(alph, rho, P, Z, W, niter = 4) {
  n <- nrow(P)
  K <- nrow(alph)
  R <- ncol(alph)
  list(
    alpha_array = array(rep(c(alph), each = niter), dim = c(niter, K, R)),
    rho_array = array(rep(rho, each = niter), dim = c(niter, R)),
    P_array = array(rep(c(P), each = niter), dim = c(niter, n, K - 1)),
    Z_array = array(rep(Z, each = niter), dim = c(niter, n)),
    W_array = array(rep(W, each = niter), dim = c(niter, length(W)))
  )
}

lbm_fixture_chains <- function(niter = 4) {
  K <- 3
  R <- 2
  n <- 6
  alpha1 <- matrix(c(1, 3, 5, 2, 4, 6), 3, 2)
  ra <- c(2L, 3L, 1L)
  ca <- c(2L, 1L)
  rho1 <- c(0.6, 0.4)
  P1 <- matrix(seq_len(n * (K - 1)) * 0.1, n, K - 1)
  Psi <- default_Psi_function(K)
  Ta <- round(Psi %*% perm_matrix_from_order(ra) %*% t(Psi), 10)
  Z1 <- rep(1:3, each = 2)
  W1 <- c(1L, 2L, 1L, 2L, 1L)
  c1 <- mk_chain_results(alpha1, rho1, P1, Z1, W1, niter)
  c2 <- mk_chain_results(
    alpha1[ra, ca],
    rho1[ca],
    P1 %*% t(Ta),
    ra[Z1],
    ca[W1],
    niter
  )
  list(c1, c2)
}

test_that("lbm_results_to_stan_draws reports a single provided chain", {
  chains <- lbm_fixture_chains()
  expect_message(
    lbm_results_to_stan_draws(chains[[1]], K = 3, R = 2),
    "Only one chain provided"
  )
})

test_that("lbm_results_to_stan_draws returns a draws_array with one chain for a single chain", {
  chains <- lbm_fixture_chains()
  out <- suppressMessages(lbm_results_to_stan_draws(chains[[1]], K = 3, R = 2))
  expect_s3_class(out, "draws_array")
  expect_equal(dim(out), c(4, 1, 31))
})

test_that("lbm_results_to_stan_draws delabels several chains onto each other", {
  chains <- lbm_fixture_chains()
  out <- suppressMessages(lbm_results_to_stan_draws(chains, K = 3, R = 2))
  expect_equal(dim(out), c(4, 2, 31))
  df <- as.data.frame(posterior::as_draws_df(out))
  sel <- !grepl("^\\.", names(df))
  expect_equal(
    unname(as.matrix(df[df$.chain == 1, sel])),
    unname(as.matrix(df[df$.chain == 2, sel])),
    tolerance = 1e-9
  )
})

test_that("lbm_results_to_stan_draws can skip the delabeling step", {
  chains <- lbm_fixture_chains()
  out <- lbm_results_to_stan_draws(chains, K = 3, R = 2, delabel_switch = FALSE)
  df <- as.data.frame(posterior::as_draws_df(out))
  expect_equal(unique(df[df$.chain == 2, "alpha[1,1]"]), 4)
  expect_equal(unique(df[df$.chain == 1, "alpha[1,1]"]), 1)
})

test_that("lbm_results_to_stan_draws applies burn-in and thinning", {
  chains <- lbm_fixture_chains(niter = 20)
  out <- suppressMessages(lbm_results_to_stan_draws(
    chains,
    K = 3,
    R = 2,
    apply_burnin_thinning = TRUE
  ))
  expect_equal(dim(out)[1], 1)
})

# ---- build_covariance_matrix -------------------------------------------------

test_that("build_covariance_matrix computes the per-chain covariance of the mean positions", {
  niter <- 8
  nch <- 2
  n <- 6
  K <- 3
  P1 <- matrix(seq_len(n * (K - 1)) * 0.1, n, K - 1)
  var_order <- c(
    paste0("P[", rep(1:n, K - 1), ",", rep(1:(K - 1), each = n), "]"),
    "sigma2"
  )
  arr <- array(
    0,
    dim = c(niter, nch, length(var_order)),
    dimnames = list(NULL, NULL, var_order)
  )
  for (it in 1:niter) {
    for (ch in 1:nch) {
      arr[it, ch, ] <- c(c(P1) * it, 1)
    }
  }
  dr <- posterior::as_draws_array(arr)

  lifecycle::expect_deprecated(
    build_covariance_matrix(dr, n = n, K = K)
  )
  out <- suppressWarnings(build_covariance_matrix(dr, n = n, K = K))
  expect_equal(dim(out), c(n, n, nch))

  # manual expectation: per chain, the coordinate-space covariance of the
  # iteration-wise mean positions, across the n individuals
  Pcols <- paste0("P[", rep(1:n, K - 1), ",", rep(1:(K - 1), each = n), "]")
  df <- as.data.frame(posterior::as_draws_df(dr))
  for (ch in 1:nch) {
    M <- as.matrix(df[df$.chain == ch, Pcols]) # iterations x P vars
    Pbar <- matrix(colMeans(M), nrow = n) # (individual, coordinate)
    expect_equal(unname(out[,, ch]), unname(cov(t(Pbar))), tolerance = 1e-9)
  }

  corr <- suppressWarnings(build_covariance_matrix(
    dr,
    n = n,
    K = K,
    corr = TRUE
  ))
  expect_equal(dim(corr), c(n, n, nch))
})
