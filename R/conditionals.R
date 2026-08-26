# Tool constants
TOL <- 1e-6

#' Conditional distribution for Pi (fast version)
#'
#' Computes the conditional distribution for a single row of P given the rest of P, Theta, and sigma2
#' @param P a matrix of size \eqn{n_1 \times K-1} specifying latent positions
#' @param Theta precision matrix (inverse of Sigma)
#' @param sigma2 variance parameter
#' @param i index of the row to compute the conditional for
#' @return A list with mean and covariance of the conditional distribution
cond_Pi_fast <- function(P, Theta, sigma2, i) {
    n <- nrow(P)
    K1 <- ncol(P)

    idx <- setdiff(1:n, i)

    Theta_ii <- Theta[i, i]
    Theta_i_rest <- Theta[i, idx]
    P_rest <- P[idx, , drop = FALSE]

    mean_i <- -(1 / Theta_ii) * (Theta_i_rest %*% P_rest)
    cov_i <- (sigma2 / Theta_ii) * diag(K1)

    list(mean = as.vector(mean_i), cov = cov_i)
}

#' Conditional distribution for Pi given P minus i and sigma
#'
#' Computes the conditional distribution for a single row of P given the rest of P, Sigma, and sigma2
#' @param P a matrix of size \eqn{n_1 \times K-1} specifying latent positions
#' @param Sigma covariance matrix
#' @param sigma2 variance parameter
#' @param i index of the row to compute the conditional for
#' @return A list with mean and covariance of the conditional distribution
cond_Pi_given_P_min_i_sigma <- function(P, Sigma, sigma2, i) {
    n <- nrow(P)
    K_minus_1 <- ncol(P)
    minus_i <- setdiff(1:n, i)

    Si <- Sigma[i, minus_i, drop = FALSE] %*% solve(Sigma[minus_i, minus_i, drop = FALSE])

    mean_i <- Si %*% P[minus_i, , drop = FALSE]

    cov_i <- sigma2 * as.numeric(Sigma[i, i] - Si %*% Sigma[-i, i]) * diag(1, nrow = K_minus_1)

    return(list(mean = mean_i, cov = cov_i))
}

#' Sample from the conditional distribution of Pi
#'
#' Samples a single row of P from its conditional distribution given the rest of P, Sigma, and sigma2
#' @param P a matrix of size \eqn{n_1 \times K-1} specifying latent positions
#' @param Sigma covariance matrix
#' @param sigma2 variance parameter
#' @param i index of the row to sample
#' @return Sampled row of P
sample_Pi_given <- function(P, Sigma, sigma2, i) {
    res <- cond_Pi_given_P_min_i_sigma(P, Sigma, sigma2, i)

    # Draw from multivariate normal
    Pi_sample <- MASS::mvrnorm(
        n = 1,
        mu = res$mean,
        Sigma = res$cov
    )

    return(Pi_sample)
}

#' Compute categorical distribution in ilr coordinates
#'
#' Computes the categorical distribution for a given row index and latent position
#' @param Zi row index
#' @param Pi latent position
#' @return probability of the categorical distribution
cat_dist_ilr_given_Pi <- function(Zi, Pi) {
    probs <- pivotCoordInv(matrix(Pi, nrow = 1))
    return(probs[Zi])
}

#' Compute posterior parameters for inverse gamma distribution
#'
#' Computes the posterior parameters for the inverse gamma distribution of sigma2
#' @param alpha_0 shape parameter of the prior inverse gamma
#' @param beta_0 rate parameter of the prior inverse gamma
#' @param P a matrix of size \eqn{n_1 \times K-1} specifying latent positions
#' @param Theta precision matrix (inverse of Sigma)
#' @return A list with alpha and beta parameters of the posterior inverse gamma
posterior_param_inv_gamma <- param_sigma2_given_P <- function(alpha_0, beta_0, P, Theta) {
    return(list(alpha = alpha_0 + (nrow(P) / 2), beta = beta_0 + 0.5 * sum(diag(t(P) %*% Theta %*% P))))
}

#' Sample from the Inverse-Gamma
#'
#' Explicit inverse-gamma sampler using the identity:
#' if X ~ IG(shape = a, rate = b) then 1 / X ~ Gamma(shape = a, rate = b)
#' @importFrom stats rgamma
#' @param shape the shape parameter of the inverse gamma
#' @param rate the rate parameter of the inverse gamma
#' @return a sample drawn from the inverse gamma distribution
sample_inv_gamma_rate <- sample_sigma2_given_P <- function(shape, rate) {
    return(1 / rgamma(n = 1, shape = shape, rate = rate))
}

## CLASSICAL LBM Z and pi


# pi | Z

#' Compute posterior parameters for the Dirichlet of row groups proportions
#'
#' Computes the posterior parameters of the Dirichlet distribution of
#' \eqn{\pi \mid Z}
#'
#' @param etas a vector of size K, the prior parameters of the Dirichlet
#' @param Z a matrix of size \eqn{n_1 \times K} with a single 1 per line
#' indicating the membership of row node \eqn{i}, \eqn{Z_{i,k} = 1} if \eqn{i} is in group \eqn{k} 0 else
#'
#' @return a vector of size K with the updated Dirichlet parameters
param_pi_given_Z <- function(etas, Z) {
    etas + colSums(Z)
}

#' Sample from the posterior \eqn{\pi \mid Z}
#'
#' A function to sample from the posterior distribution \eqn{\pi \mid Z}
#' which is a Dirichlet with the posterior parameters
#'
#' @param etas_post a vector of size K containing the Dirichlet parameters from which to sample
#'
#' @seealso [param_pi_given_Z()] for the computations of the posterior parameters
sample_pi_given_Z <- function(etas_post) {
    as.vector(MCMCpack::rdirichlet(n = 1, alpha = etas_post))
}

# Z | alpha, pi, Y, W
# Multinomial prob to normalize

#' Compute multinomial probabilities of Z in the classical Poisson LBM
#'
#' Computes the normalized conditional probabilities
#' \eqn{Z_i \mid \alpha, W, Y, \pi} of each row-node block membership
#' under a Poisson latent block model.
#'
#' @param Y Non-negative integer matrix of observed counts.
#' @param alpha a matrix (\eqn{K \times R}) of connectivity coefficients
#' @param W a matrix of size \eqn{n_2 \times R} with a single 1 per line
#' indicating the membership of column node \eqn{j}
#' @param pi a vector of size K containing the row block proportions
#' @param tol Numeric; a small number used to clamp probabilities to
#' \code{[tol, 1-tol]}. Default to \link{.Machine}$double.eps based tolerance
#'
#' @return a matrix (\eqn{n_1 \times K}) of normalized membership probabilities
param_multinom_probs_Z_poisson <- function(Y, alpha, W, pi, tol = TOL) {
    R_W <- Y %*% W
    N_W <- diag(colSums(W))
    unormalized_log_probs <- matrix(1, nrow = nrow(Y)) %*% log(pi) + R_W %*% log(t(alpha)) - matrix(1, nrow = nrow(Y), ncol = ncol(W)) %*% N_W %*% t(alpha)

    return(row_normalize_matrix(unormalized_log_probs, tol = tol))
}

#' Sample row-block memberships Z
#'
#' Draws the block label of each row node independently from its
#' categorical distribution.
#'
#' @param probs a matrix (\eqn{n_1 \times K}) of per-row membership probabilities
#'
#' @return an integer vector of length \eqn{n_1} with one sampled group label per row node
sample_Z_given_alpha_pi_Y_W <- function(probs) {
    sapply(seq_len(nrow(probs)), function(i) {
        sample.int(n = ncol(probs), size = 1, replace = TRUE, prob = probs[i, ])
    })
}

## END OF CLASSICAL LBM

# rho | W

#' Compute the posterior parameters for the Dirichlet of row groups
#'
#' @param gammas a vector of size R (the number of column groups),
#' the prior of the Dirichlet
#' @param W a matrix of size \eqn{n_2 \times R} with a single 1 per line
#' indicating the membership of column node \eqn{j}, \eqn{W_{j,r} = 1} if \eqn{j} is in group \eqn{r} 0 else
#'
#' @return a vector of size R with the updated Dirichlet parameters
param_rho_given_W <- function(gammas, W) {
    gammas + colSums(W)
}

#' Sample from the posterior \eqn{\rho\mid W}
#'
#' A function to sample from the posterior distribution \eqn{\rho\mid W} which is a Dirichlet of posterior parameters
#'
#' @param gammas_post a vector of size R containing the Dirichlet parameter from which to sample
#'
#' @importFrom MCMCpack rdirichlet
#'
#' @seealso [param_rho_given_W()] for the computations of the posterior parameters
sample_rho_given_W <- function(gammas_post) {
    as.vector(MCMCpack::rdirichlet(n = 1, alpha = gammas_post))
}

# Z | alpha, P, Y, W
# Multinomial prob to normalize

#' Compute multinomial probabilities of Z in the latent covariance Poisson LBM
#'
#' Computes the normalized conditional probabilities
#' \eqn{Z_i \mid \alpha, W, Y, P} of each row-node block membership
#' under a Poisson latent block model with latent correlated positions.
#' The probabilities are derived from the latent positions through their
#' inverse pivot coordinates.
#'
#' @param Y Non-negative integer matrix of observed counts.
#' @param alpha a matrix (\eqn{K \times R}) of connectivity coefficients
#' @param W a matrix of size \eqn{n_2 \times R} with a single 1 per line
#' indicating the membership of column node \eqn{j}
#' @param P a matrix of size \eqn{n_1 \times K-1} specifying latent positions
#' @param tol Numeric; a small number used to clamp probabilities to
#' \code{[tol, 1-tol]}. Default to \link{.Machine}$double.eps based tolerance
#'
#' @return a matrix (\eqn{n_1 \times K}) of normalized membership probabilities
param_multinom_probs_Z_cov_poisson <- function(Y, alpha, W, P, tol = TOL) {
    R_W <- Y %*% W
    N_W <- diag(colSums(W))
    unormalized_log_probs <- log(pivotCoordInv(P)) + R_W %*% log(t(alpha)) - matrix(1, nrow = nrow(Y), ncol = ncol(W)) %*% N_W %*% t(alpha)

    return(row_normalize_matrix(unormalized_log_probs, tol = tol))
}

#' Sample row-block memberships Z
#'
#' Draws the block label of each row node independently from its
#' categorical distribution.
#'
#' @param probs a matrix (\eqn{n_1 \times K}) of per-row membership probabilities
#'
#' @return an integer vector of length \eqn{n_1} with one sampled group label per row node
sample_Z_given_alpha_P_Y_W <- function(probs) {
    sapply(seq_len(nrow(probs)), function(i) {
        sample.int(n = ncol(probs), size = 1, replace = TRUE, prob = probs[i, ])
    })
}

# W | alpha, rho, Y, Z

#' Compute multinomial probabilities of W in the Poisson LBM
#'
#' Computes the normalized conditional probabilities
#' \eqn{W_j \mid \alpha, Z, Y, \rho} of each column-node block membership
#' under a Poisson latent block model.
#'
#' @param Y Non-negative integer matrix of observed counts.
#' @param alpha a matrix (\eqn{K \times R}) of connectivity coefficients
#' @param Z a matrix of size \eqn{n_1 \times K} with a single 1 per line
#' indicating the membership of row node \eqn{i}
#' @param rho a vector of size R containing the column block proportions
#' @param tol Numeric; a small number used to clamp probabilities to
#' \code{[tol, 1-tol]}. Default to \link{.Machine}$double.eps based tolerance
#'
#' @return a matrix (\eqn{n_2 \times R}) of normalized membership probabilities
param_multinom_probs_W_poisson <- function(Y, alpha, Z, rho, tol = TOL) {
    R_Z <- t(Y) %*% Z
    N_Z <- diag(colSums(Z))
    unormalized_log_probs <- matrix(1, nrow = ncol(Y)) %*% log(rho) + R_Z %*% log(alpha) - matrix(1, nrow = ncol(Y), ncol = ncol(Z)) %*% N_Z %*% alpha

    return(row_normalize_matrix(unormalized_log_probs, tol = tol))
}

#' Sample column-block memberships W
#'
#' Draws the block label of each column node independently from its
#' categorical distribution.
#'
#' @param probs a matrix (\eqn{n_2 \times R}) of per-column membership probabilities
#'
#' @return an integer vector of length \eqn{n_2} with one sampled group label per column node
sample_W_given_alpha_rho_Y_Z <- function(probs) {
    sapply(seq_len(nrow(probs)), function(i) {
        sample.int(n = ncol(probs), size = 1, replace = TRUE, prob = probs[i, ])
    })
}


# alpha | Y,Z,W

#' Compute posterior parameters for the Gamma connectivity coefficients
#'
#' Computes the shape and rate matrices of the full conditional
#' \eqn{\alpha \mid Y, Z, W} whose entries follow independent Gamma distributions.
#'
#' @param a0 prior shape parameter of the Gamma
#' @param b0 prior rate parameter of the Gamma
#' @param Y Non-negative integer matrix of observed counts.
#' @param Z a matrix of size \eqn{n_1 \times K} with a single 1 per line
#' indicating the membership of row node \eqn{i}
#' @param W a matrix of size \eqn{n_2 \times R} with a single 1 per line
#' indicating the membership of column node \eqn{j}
#'
#' @return a list with `shape` and `rate` matrices (\eqn{K \times R}) of the posterior Gamma parameters
param_alpha_given_Y_Z_W_poisson <- function(a0, b0, Y, Z, W) {
    return(list(shape = a0 + t(Z) %*% Y %*% W, rate = b0 + t(Z) %*% matrix(1, nrow = nrow(Z), ncol = nrow(W)) %*% W))
}

#' Sample the connectivity matrix alpha
#'
#' Draws each entry of the connectivity matrix independently from its
#' Gamma full-conditional distribution.
#'
#' @param shape a matrix (\eqn{K \times R}) containing the shape parameter from which to sample
#' @param rate a matrix (\eqn{K \times R}) containing the rate parameter from which to sample
#'
#' @return a matrix (\eqn{K \times R}) of sampled connectivity coefficients
sample_alpha_given_Y_Z_W_poisson <- function(shape, rate) {
    matrix(rgamma(length(shape), shape = shape, rate = rate), nrow = nrow(shape))
}
