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
  P <- t(mvtnorm::rmvnorm(n = K - 1, mean = rep(0, nrow(Sigma)), sigma = sigma2 * Sigma))

  probs <- pivotCoordInv(P)

  Z <- sapply(seq_len(nrow(probs)), function(i) {
    sample.int(n = K, size = 1, replace = TRUE, prob = probs[i, ])
  })

  return(list(Z = Z, P = P, probs = probs))
}
