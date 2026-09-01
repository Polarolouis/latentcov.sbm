#' Simulate P and Z matrices
#'
#' Simulates the P and Z matrices for the Latent Covariance Stochastic Block Model
#' @param K number of clusters
#' @param Sigma covariance matrix
#' @param sigma2 variance parameter
#'
#' @return A list containing Z, P, and probs
#' @importFrom stats rmultinom
#' @export
simulate_P_and_Z <- function(K, Sigma, sigma2) {
  # ---- Guard rails ---------------------------------------------------------

  if (!is.numeric(K) || length(K) != 1L || is.na(K) || K != round(K) || K < 2L) {
    cli::cli_abort(c(
      "x" = "`K` must be a single integer larger or equal to 2 (the number of row clusters).",
      "i" = "You provided {.val {K}}."
    ))
  }

  if (!is.numeric(sigma2) || length(sigma2) != 1L || is.na(sigma2) || !is.finite(sigma2) || sigma2 <= 0) {
    cli::cli_abort(c(
      "x" = "`sigma2` must be a single strictly positive scalar (the variance of the latent coordinates).",
      "i" = "You provided {.val {sigma2}}."
    ))
  }

  if (!is.matrix(Sigma) || !is.numeric(Sigma)) {
    cli::cli_abort(c(
      "x" = "`Sigma` must be a numeric matrix.",
      "i" = "You provided an object of class {.cls {class(Sigma)}}."
    ))
  }

  n1 <- nrow(Sigma)
  if (n1 != ncol(Sigma)) {
    cli::cli_abort(c(
      "x" = "`Sigma` must be a square matrix.",
      "i" = "It has {n1} row(s) and {.val {ncol(Sigma)}} column(s)."
    ))
  }

  if (n1 < K) {
    cli::cli_abort(c(
      "x" = "`Sigma` must have at least `K` rows/columns.",
      "i" = "It has {n1} and `K` = {K}."
    ))
  }

  if (!isSymmetric(Sigma)) {
    cli::cli_abort(c(
      "x" = "`Sigma` must be symmetric.",
      "i" = "The largest |Sigma - t(Sigma)| entry is {.val {max(abs(Sigma - t(Sigma)))}}."
    ))
  }

  diag_dev <- max(abs(diag(Sigma) - 1))
  if (diag_dev > 1e-4) {
    cli::cli_warn(c(
      "!" = "`Sigma` does not look like a correlation matrix: its diagonal deviates from 1.",
      "i" = "max |diag(Sigma) - 1| = {.val {diag_dev}}.",
      "i" = "Since the rows of `P` have variance {.val sigma2 * Sigma[i, i]}, a diagonal different from 1 breaks the marginal standardisation of the latent coordinates."
    ))
  }

  e <- eigen(Sigma, symmetric = TRUE)
  lambda <- e$values
  lambda_scale <- max(1, max(abs(lambda)))
  tol_min <- 1e-8 * lambda_scale

  if (min(lambda) < -tol_min) {
    cli::cli_abort(c(
      "x" = "`Sigma` is not positive semi-definite, hence cannot be a correlation matrix.",
      "i" = "Its smallest eigenvalue is {.val {min(lambda)}}.",
      "i" = "For a block-structured matrix (within-block {.val rho_w}, between-block {.val rho_b}) the minimal feasible {.code rho_b} is the unique negative root of {.code 1 + rho_b * sum(m_q / ((rho_w - rho_b) * m_q + (1 - rho_w))) = 0} ({.code rho_b = -(1 + (m - 1) * rho_w) / (m * (Q - 1))} for equal blocks).",
      "i" = "Project `Sigma` onto the correlation-matrix cone first (e.g. a Higham nearest-correlation-matrix algorithm such as the one provided in the experiments repository)."
    ))
  }

  if (min(lambda) < tol_min) {
    cli::cli_warn(c(
      "!" = "`Sigma` is (numerically) singular: its smallest eigenvalue is close to 0 ({.val {min(lambda)}}).",
      "i" = "An eigen (SVD-like) based sampler is used, robust to positive semi-definite matrices.",
      "i" = "For better-conditioned draws consider making `Sigma` positive definite (small ridge on the diagonal or spectral clipping)."
    ))
  }

  # ---- Sampling ------------------------------------------------------------

  if (min(lambda) >= tol_min) {
    # Positive definite: use the classic (chol based) matrix-normal draw.
    P <- rmatrixnormal(M = matrix(0, n1, K - 1), U = Sigma, V = sigma2 * diag(1, K - 1))
  } else {
    # (Numerically) singular: eigen-based draw, vec(P) ~ N(0, sigma2 * I %x% Sigma).
    A <- e$vectors %*% diag(sqrt(pmax(lambda, 0)), n1, n1) %*% t(e$vectors) # A A^T = Sigma
    X <- matrix(stats::rnorm(n1 * (K - 1)), n1, K - 1)
    P <- sqrt(sigma2) * (A %*% X)
  }

  probs <- pivotCoordInv(P)

  Z <- sapply(seq_len(nrow(probs)), function(i) {
    sample.int(n = K, size = 1, replace = TRUE, prob = probs[i, ])
  })

  return(list(Z = Z, P = P, probs = probs))
}
