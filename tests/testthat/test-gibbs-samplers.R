# Gibbs samplers: structure, deterministic `known_*` branches, NA handling.

.ny <- rbind(
  c(0, 1, 2, 0, 1), c(3, 4, 0, 2, 1), c(5, 0, 1, 3, 4),
  c(0, 0, 1, 2, 1), c(2, 1, 0, 0, 3), c(1, 3, 2, 1, 0)
)
.ny_na <- .ny
.ny_na[2, 3] <- NA
.ny_na[5, 1] <- NA

.z_init <- .one_hot(rep(c(1, 2, 3), each = 2), 3)
.w_init <- .one_hot(c(1, 2, 1, 2, 1), 2)
.hp <- list(etas_0 = rep(2, 3), gammas_0 = rep(2, 2), a0 = 1, b0 = 1)

# ---- classical Poisson LBM ---------------------------------------------------

test_that("classical Gibbs returns a finite draws_array of the right shape", {
  set.seed(10)
  out <- gibbs_sampling_lbm_poisson(.ny, K = 3, R = 2, niter = 5, priors_hyper_params = .hp, init_Z = .z_init, init_W = .w_init)
  expect_s3_class(out, "draws_array")
  df <- as.data.frame(posterior::as_draws_df(out))
  expect_equal(nrow(df), 5)
  expect_false(anyNA(as.data.frame(posterior::as_draws_df(out))))
  v <- posterior::variables(out)
  expect_true(all(c("rho[1]", "rho[2]", "pi[1]", "pi[2]", "pi[3]", "alpha[1,1]") %in% v))
})

test_that("classical Gibbs handles NAs without leaking them", {
  set.seed(11)
  out <- gibbs_sampling_lbm_poisson(.ny_na, K = 3, R = 2, niter = 5, priors_hyper_params = .hp, init_Z = .z_init, init_W = .w_init)
  expect_false(anyNA(as.data.frame(posterior::as_draws_df(out))))
})

test_that("known_pi pins the pi trajectory", {
  kp <- c(0.5, 0.3, 0.2)
  set.seed(12)
  out <- gibbs_sampling_lbm_poisson(.ny, K = 3, R = 2, niter = 5, priors_hyper_params = .hp, init_Z = .z_init, init_W = .w_init, known_pi = kp)
  df <- as.data.frame(posterior::as_draws_df(out))
  for (k in 1:3) {
    expect_equal(as.vector(df[[paste0("pi[", k, "]")]]), rep(kp[k], 5))
  }
})

test_that("known_alpha pins the alpha trajectory", {
  set.seed(13)
  out <- gibbs_sampling_lbm_poisson(.ny, K = 3, R = 2, niter = 5, priors_hyper_params = .hp, init_Z = .z_init, init_W = .w_init, known_alpha = alpha_fx)
  df <- as.data.frame(posterior::as_draws_df(out))
  for (k in 1:3) {
    for (r in 1:2) {
      expect_equal(as.vector(df[[paste0("alpha[", k, ",", r, "]")]]), rep(alpha_fx[k, r], 5))
    }
  }
})

test_that("known_Z and known_W matrices pin the membership trajectories", {
  set.seed(14)
  out <- gibbs_sampling_lbm_poisson(.ny,
    K = 3, R = 2, niter = 5, priors_hyper_params = .hp,
    init_Z = .z_init, init_W = .w_init, known_Z = .z_init, known_W = .w_init
  )
  df <- as.data.frame(posterior::as_draws_df(out))
  z_lab <- apply(.z_init, 1, which.max)
  w_lab <- apply(.w_init, 1, which.max)
  for (i in seq_len(nrow(.ny))) {
    expect_equal(as.vector(df[[paste0("Z[", i, "]")]]), rep(z_lab[i], 5))
  }
  for (j in seq_len(ncol(.ny))) {
    expect_equal(as.vector(df[[paste0("W[", j, "]")]]), rep(w_lab[j], 5))
  }
})

test_that("force_order stores pi and rho sorted decreasing each iteration", {
  set.seed(15)
  out <- gibbs_sampling_lbm_poisson(.ny, K = 3, R = 2, niter = 5, priors_hyper_params = .hp, init_Z = .z_init, init_W = .w_init, force_order = TRUE)
  df <- as.data.frame(posterior::as_draws_df(out))
  pi_rows <- df[, c("pi[1]", "pi[2]", "pi[3]")]
  rho_rows <- df[, c("rho[1]", "rho[2]")]
  expect_true(all(apply(pi_rows, 1, function(z) all(diff(z) <= 1e-12))))
  expect_true(all(apply(rho_rows, 1, function(z) all(diff(z) <= 1e-12))))
})

test_that("the classical Gibbs sampler logs every 10 iterations when verbose", {
  expect_message(
    gibbs_sampling_lbm_poisson(.ny, K = 3, R = 2, niter = 10, priors_hyper_params = .hp, init_Z = .z_init, init_W = .w_init, verbose = TRUE),
    "Iter : 10"
  )
  expect_silent(
    gibbs_sampling_lbm_poisson(.ny, K = 3, R = 2, niter = 3, priors_hyper_params = .hp, init_Z = .z_init, init_W = .w_init, verbose = FALSE)
  )
})

# ---- latent covariance Poisson LBM ------------------------------------------

test_that("cov Gibbs keeps sigma2 fixed at 1 by default", {
  set.seed(20)
  out <- gibbs_sampling_lbm_cov_poisson(Sigma = Sigma_fx, Y = .ny, K = 3, R = 2, niter = 5, init_Z = .z_init, init_W = .w_init)
  expect_false(anyNA(as.data.frame(posterior::as_draws_df(out))))
  df <- as.data.frame(posterior::as_draws_df(out))
  expect_equal(as.vector(df[["sigma2"]]), rep(1, 5))
})

test_that("a numeric sigma2_fixed pins sigma2 and logs it", {
  set.seed(21)
  expect_message(
    out <- gibbs_sampling_lbm_cov_poisson(
      Sigma = Sigma_fx, Y = .ny, K = 3, R = 2, niter = 5,
      init_Z = .z_init, init_W = .w_init, sigma2_fixed = 0.5
    ),
    "Using sigma2=0.5"
  )
  df <- as.data.frame(posterior::as_draws_df(out))
  expect_equal(as.vector(df[["sigma2"]]), rep(0.5, 5))
})

test_that("sigma2_fixed = FALSE samples strictly positive sigma2", {
  set.seed(22)
  out <- gibbs_sampling_lbm_cov_poisson(
    Sigma = Sigma_fx, Y = .ny, K = 3, R = 2, niter = 15,
    init_Z = .z_init, init_W = .w_init, sigma2_fixed = FALSE
  )
  df <- as.data.frame(posterior::as_draws_df(out))
  s2 <- as.vector(df[["sigma2"]])
  expect_true(all(s2 > 0) && all(is.finite(s2)))
})

test_that("an identity P_sampler pins P and known_* pin their trajectories", {
  set.seed(23)
  out <- gibbs_sampling_lbm_cov_poisson(
    Sigma = Sigma_fx, Y = .ny, K = 3, R = 2, niter = 5,
    init_Z = .z_init, init_W = .w_init,
    known_alpha = alpha_fx, known_P = P_fx, known_Z = .z_init, known_W = .w_init,
    P_sampler = function(P, ...) P
  )
  df <- as.data.frame(posterior::as_draws_df(out))
  for (i in seq_len(nrow(.ny))) {
    for (crd in 1:2) {
      expect_equal(as.vector(df[[paste0("P[", i, ",", crd, "]")]]), rep(unname(P_fx[i, crd]), 5))
    }
  }
  for (k in 1:3) {
    for (r in 1:2) {
      expect_equal(as.vector(df[[paste0("alpha[", k, ",", r, "]")]]), rep(alpha_fx[k, r], 5))
    }
  }
  for (i in seq_len(nrow(.ny))) {
    expect_equal(as.vector(df[[paste0("Z[", i, "]")]]), rep(apply(.z_init, 1, which.max)[i], 5))
  }
})

test_that("cov Gibbs handles NAs and reports verbose iterations", {
  set.seed(24)
  out <- gibbs_sampling_lbm_cov_poisson(Sigma = Sigma_fx, Y = .ny_na, K = 3, R = 2, niter = 5, init_Z = .z_init, init_W = .w_init)
  expect_false(anyNA(as.data.frame(posterior::as_draws_df(out))))
  expect_message(
    gibbs_sampling_lbm_cov_poisson(Sigma = Sigma_fx, Y = .ny, K = 3, R = 2, niter = 10, init_Z = .z_init, init_W = .w_init, verbose = TRUE),
    "Iter : 10"
  )
})

# ---- mixed-membership latent covariance Poisson LBM --------------------------

test_that("mixed-membership Gibbs returns soft membership draws summing to 1", {
  set.seed(30)
  zsoft <- ilrInv(P_fx)
  out <- gibbs_sampling_lbm_cov_poisson_mixed_membership(
    Sigma = Sigma_fx, Y = .ny, init_Z = zsoft, init_W = .w_init,
    K = 3, R = 2, niter = 5, known_alpha = alpha_fx, known_P = P_fx, known_W = .w_init
  )
  expect_false(anyNA(as.data.frame(posterior::as_draws_df(out))))
  v <- posterior::variables(out)
  expect_true(any(grepl("^Z_post_probs", v)))
  df <- as.data.frame(posterior::as_draws_df(out))
  for (i in seq_len(nrow(.ny))) {
    cols <- paste0("Z_post_probs[", i, ",", 1:3, "]")
    expect_equal(unname(rowSums(df[, cols])), rep(1, 5), tolerance = 1e-10)
  }
})

test_that("mixed-membership GG handles NAs (soft alpha | Y, Z, W statistics)", {
  set.seed(31)
  zsoft <- ilrInv(P_fx)
  out <- gibbs_sampling_lbm_cov_poisson_mixed_membership(
    Sigma = Sigma_fx, Y = .ny_na, init_Z = zsoft, init_W = .w_init,
    K = 3, R = 2, niter = 5, known_alpha = alpha_fx, known_P = P_fx, known_W = .w_init
  )
  expect_false(anyNA(as.data.frame(posterior::as_draws_df(out))))
})

test_that("mixed-membership Gibbs samples alpha and W when unknown", {
  set.seed(32)
  zsoft <- ilrInv(P_fx)
  out <- gibbs_sampling_lbm_cov_poisson_mixed_membership(
    Sigma = Sigma_fx, Y = .ny, init_Z = zsoft, init_W = .w_init,
    K = 3, R = 2, niter = 5, known_P = P_fx
  )
  expect_false(anyNA(as.data.frame(posterior::as_draws_df(out))))
})

test_that("mixed-membership Gibbs samples alpha and W from NA-containing data", {
  set.seed(34)
  zsoft <- ilrInv(P_fx)
  out <- gibbs_sampling_lbm_cov_poisson_mixed_membership(
    Sigma = Sigma_fx, Y = .ny_na, init_Z = zsoft, init_W = .w_init,
    K = 3, R = 2, niter = 5, known_P = P_fx
  )
  expect_false(anyNA(as.data.frame(posterior::as_draws_df(out))))
})

test_that("mixed-membership Gibbs samples sigma2 when not fixed", {
  set.seed(33)
  zsoft <- ilrInv(P_fx)
  out <- gibbs_sampling_lbm_cov_poisson_mixed_membership(
    Sigma = Sigma_fx, Y = .ny, init_Z = zsoft, init_W = .w_init,
    K = 3, R = 2, niter = 10, known_alpha = alpha_fx, known_P = P_fx, known_W = .w_init, sigma2_fixed = FALSE
  )
  df <- as.data.frame(posterior::as_draws_df(out))
  expect_true(all(as.vector(df[["sigma2"]]) > 0))
})

test_that("known_membership_to_labels turns a one-hot matrix into labels", {
  memb <- rbind(c(0, 1, 0), c(1, 0, 0), c(0, 0, 1))
  expect_equal(known_membership_to_labels(memb), c(2L, 1L, 3L))
})

test_that("known_membership_to_labels passes a label vector through", {
  expect_equal(known_membership_to_labels(c(2, 1, 3)), c(2L, 1L, 3L))
})

# ---- coverage: initial proportion draws, force_order, auto_save, chains ----

test_that("classical Gibbs initializes NULL proportions sorted by force_order", {
  set.seed(40)
  out <- gibbs_sampling_lbm_poisson(.ny,
    K = 3, R = 2, niter = 5, priors_hyper_params = .hp,
    init_Z = NULL, init_W = NULL, force_order = TRUE
  )
  df <- as.data.frame(posterior::as_draws_df(out))
  expect_true(all(apply(df[, c("pi[1]", "pi[2]", "pi[3]")], 1, function(z) all(diff(z) <= 1e-12))))
  expect_true(all(apply(df[, c("rho[1]", "rho[2]")], 1, function(z) all(diff(z) <= 1e-12))))
})

test_that("classical Gibbs auto-saves when an iteration lands on a multiple of 100", {
  set.seed(41)
  path <- tempfile(fileext = ".Rds")
  out <- gibbs_sampling_lbm_poisson(.ny,
    K = 3, R = 2, niter = 100, priors_hyper_params = .hp,
    init_Z = .z_init, init_W = .w_init, auto_save_path = path
  )
  expect_true(file.exists(path))
  expect_false(anyNA(as.data.frame(posterior::as_draws_df(out))))
})

test_that("cov Gibbs draws initial proportions when init_W and init_Z are NULL", {
  set.seed(42)
  out <- gibbs_sampling_lbm_cov_poisson(Sigma = Sigma_fx, Y = .ny, K = 3, R = 2, niter = 5, init_Z = NULL, init_W = NULL)
  expect_false(anyNA(as.data.frame(posterior::as_draws_df(out))))
})

test_that("cov Gibbs auto-saves every 100 iterations", {
  set.seed(43)
  path <- tempfile(fileext = ".Rds")
  out <- gibbs_sampling_lbm_cov_poisson(
    Sigma = Sigma_fx, Y = .ny, K = 3, R = 2, niter = 100,
    init_Z = .z_init, init_W = .w_init, auto_save_path = path
  )
  expect_true(file.exists(path))
})

test_that("mixed-membership Gibbs logs its state when verbose", {
  set.seed(44)
  zsoft <- ilrInv(P_fx)
  expect_message(
    gibbs_sampling_lbm_cov_poisson_mixed_membership(
      Sigma = Sigma_fx, Y = .ny, init_Z = zsoft, init_W = .w_init,
      K = 3, R = 2, niter = 5, known_alpha = alpha_fx, known_P = P_fx, verbose = TRUE
    ),
    "state of the chain"
  )
})

test_that("mixed-membership Gibbs draws proportions when init_Z and init_W are NULL", {
  set.seed(45)
  out <- gibbs_sampling_lbm_cov_poisson_mixed_membership(
    Sigma = Sigma_fx, Y = .ny, init_Z = NULL, init_W = NULL,
    K = 3, R = 2, niter = 5, known_alpha = alpha_fx
  )
  expect_false(anyNA(as.data.frame(posterior::as_draws_df(out))))
})

test_that("mixed-membership Gibbs accepts a numeric sigma2_fixed and logs it", {
  set.seed(46)
  zsoft <- ilrInv(P_fx)
  expect_message(
    out <- gibbs_sampling_lbm_cov_poisson_mixed_membership(
      Sigma = Sigma_fx, Y = .ny, init_Z = zsoft, init_W = .w_init,
      K = 3, R = 2, niter = 5, known_alpha = alpha_fx, known_P = P_fx, sigma2_fixed = 0.5
    ),
    "Using sigma2=0.5"
  )
  df <- as.data.frame(posterior::as_draws_df(out))
  expect_equal(as.vector(df[["sigma2"]]), rep(0.5, 5))
})

test_that("mixed-membership Gibbs samples P when known_P is NULL", {
  set.seed(47)
  zsoft <- ilrInv(P_fx)
  out <- gibbs_sampling_lbm_cov_poisson_mixed_membership(
    Sigma = Sigma_fx, Y = .ny, init_Z = zsoft, init_W = .w_init,
    K = 3, R = 2, niter = 5, known_alpha = alpha_fx
  )
  expect_false(anyNA(as.data.frame(posterior::as_draws_df(out))))
})

test_that("mixed-membership Gibbs auto-saves every 100 iterations", {
  set.seed(48)
  zsoft <- ilrInv(P_fx)
  path <- tempfile(fileext = ".Rds")
  out <- gibbs_sampling_lbm_cov_poisson_mixed_membership(
    Sigma = Sigma_fx, Y = .ny, init_Z = zsoft, init_W = .w_init,
    K = 3, R = 2, niter = 100, known_alpha = alpha_fx, known_P = P_fx, auto_save_path = path
  )
  expect_true(file.exists(path))
})

test_that("the three chains wrappers run several chains and log when verbose", {
  set.seed(50)
  expect_message(
    res <- chains_gibbs_sampling_lbm_poisson(
      nchains = 2, Y = .ny, K = 3, R = 2, niter = 3,
      priors_hyper_params = .hp, init_Z = .z_init, init_W = .w_init, verbose = TRUE
    ),
    "Chains will be saved"
  )
  expect_equal(posterior::nchains(res), 2)
  set.seed(51)
  expect_message(
    res2 <- chains_gibbs_sampling_lbm_cov_poisson(
      nchains = 2, Sigma = Sigma_fx, Y = .ny, K = 3, R = 2, niter = 3,
      init_Z = .z_init, init_W = .w_init, verbose = TRUE
    ),
    "Chains will be saved"
  )
  expect_equal(posterior::nchains(res2), 2)
  set.seed(52)
  expect_message(
    res3 <- chains_gibbs_sampling_lbm_cov_poisson_mixed_membership(
      nchains = 2, Sigma = Sigma_fx, Y = .ny,
      init_Z = ilrInv(P_fx), init_W = .w_init, K = 3, R = 2, niter = 3, known_alpha = alpha_fx, verbose = TRUE
    ),
    "Chains will be saved"
  )
  expect_equal(posterior::nchains(res3), 2)
  expect_false(anyNA(as.data.frame(posterior::as_draws_df(res3))))
})
