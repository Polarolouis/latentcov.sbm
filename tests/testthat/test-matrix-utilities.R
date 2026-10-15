# --- one-hot encodings --------------------------------------------------------

test_that("onehot_encode builds row-stochastic indicator matrices", {
  enc <- onehot_encode(c(1, 2, 1), K = 3)
  expect_equal(dim(enc), c(3, 3))
  expect_equal(enc[, 1], c(1, 0, 1))
  expect_equal(enc[, 2], c(0, 1, 0))
  expect_equal(rowSums(enc), rep(1, 3))
})

test_that("onehot_encode infers K from the unique labels by default", {
  expect_equal(dim(onehot_encode(c(1, 2, 3, 2, 1))), c(5, 3))
  expect_equal(ncol(onehot_encode(factor(c("a", "b", "a")))), 2)
})

test_that(".one_hot matches onehot_encode positions", {
  expect_equal(.one_hot(c(2, 1, 2), 3), onehot_encode(c(2, 1, 2), K = 3))
  expect_equal(.one_hot(c(3, 3), 4)[, 3], c(1, 1))
})

test_that(".rev_one_hot recovers the labels", {
  X <- .one_hot(c(2, 1, 3), 3)
  expect_equal(.rev_one_hot(X), c(2, 1, 3))
})

# --- permutation matrix -------------------------------------------------------

test_that("perm_matrix_from_order builds P with x[order] = P %*% x", {
  x <- c(5, 9, 7)
  order <- c(2, 1, 3)
  P <- perm_matrix_from_order(order)
  expect_equal(as.vector(P %*% x), x[order])
  expect_equal(rowSums(P), rep(1, 3))
  expect_equal(colSums(P), rep(1, 3))
})

test_that("perm_matrix_from_order handles the identity", {
  expect_equal(perm_matrix_from_order(1:4), diag(4))
})

# --- mse ----------------------------------------------------------------------

test_that("mse returns the sum of squared differences", {
  expect_equal(mse(c(1, 2), c(3, 6)), 20)
  expect_equal(mse(matrix(0, 2, 2), matrix(0, 2, 2)), 0)
})

# --- NA masking ---------------------------------------------------------------

test_that("mask_NA splits Y into YnoNA and a 0/1 mask", {
  res <- mask_NA(Y_fx, replace_value = 7)
  expect_equal(res$mask, mask_fx)
  expect_equal(res$YnoNA[is.na(Y_fx)], rep(7, sum(is.na(Y_fx))))
  expect_equal(res$YnoNA[!is.na(Y_fx)], Y_fx[!is.na(Y_fx)])
})

test_that("mask_observed zeroes out the masked (missing) entries, no NA leak", {
  out <- mask_observed(Y_fx, mask_fx)
  expect_false(anyNA(out))
  expect_true(all(out[which(mask_fx == 0)] == 0))
  expect_equal(out[which(mask_fx == 1)], Y_fx[!is.na(Y_fx)])
})

test_that("mask_observed is robust whether Y already has 0-substituted NAs or not", {
  raw <- mask_observed(Y_fx, mask_fx)
  clean <- mask_NA(Y_fx)$YnoNA
  expect_equal(mask_observed(clean, mask_fx), raw)
})

# --- auto_save ----------------------------------------------------------------

test_that("auto_save writes an .Rds and returns its path invisibly", {
  obj <- list(a = 1, b = 2)
  path <- tempfile(fileext = ".Rds")
  expect_message(auto_save(obj, path, verbose = TRUE), "Saving")
  expect_true(file.exists(path))
  expect_equal(readRDS(path), obj)
  expect_invisible(auto_save(obj, tempfile(fileext = ".Rds"), verbose = FALSE))
})

# --- .onLoad ------------------------------------------------------------------

test_that(".onLoad registers the verbose option", {
  op <- options(latentcov.sbm.verbose = NULL)
  on.exit(options(op), add = TRUE)

  .onLoad(libname = NULL, pkgname = NULL)
  expect_false(getOption("latentcov.sbm.verbose", default = TRUE))
})

test_that(".onLoad does not override a user-set option", {
  op <- options(latentcov.sbm.verbose = TRUE)
  on.exit(options(op), add = TRUE)

  .onLoad(libname = NULL, pkgname = NULL)
  expect_true(getOption("latentcov.sbm.verbose", default = FALSE))
})
