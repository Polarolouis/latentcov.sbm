estimate_and_select_model <- function(
  Y,
  Sigma = NULL,
  emission = "poisson",
  niter = 1000L,
  priors_hyper_params = default_priors(emission)
) {}

#' Compute the exact ICL for the integrated alpha parameter
#'
#' @inheritParams icl_logL_rho
#'
#' @return the integrated term of the exact ICL for the emission
#' parameters
#'
#' @keywords internal
icl_logL_alpha <- function(
  emission = "poisson",
  K,
  R,
  priors,
  Y,
  Z,
  W,
  mask = (!is.na(Y)) * 1L
) {
  check_emission_distribution(emission)

  suff_emission_stats <- sufficient_stats_emission(Y, Z, W, emission, mask)

  switch(
    emission,
    "poisson" = {
      a0 <- priors[["a0"]]
      b0 <- priors[["b0"]]
      S <- suff_emission_stats[["S"]]
      E <- suff_emission_stats[["E"]]
      logfactY <- suff_emission_stats[["logfactY"]]
      constant_term <- K * R * (a0 * log(b0) - lgamma(a0))
      stats_term <- sum(lgamma(a0 + S) - (a0 + S) * log(b0 + E))

      return(constant_term + stats_term - logfactY)
    },
    "bernoulli" = {
      a0 <- priors[["a0"]]
      b0 <- priors[["b0"]]
      S <- suff_emission_stats[["S"]]
      F <- suff_emission_stats[["F"]]
      constant_term <- -K * R * lbeta(a0, b0)
      stats_term <- sum(lbeta(a0 + S, b0 + F))

      return(constant_term + stats_term)
    }
  )
}

#' Compute the exact ICL for the integrated rho parameter
#'
#' @param R the number of column groups
#' @param priors the priors list
#' @param W a matrix of size \eqn{n_2 \times R} indicating hard or
#' soft-memberhips for the column nodes
#'
#' @return the integrated term of the exact ICL for the column
#' parameters
#'
#' @keywords internal
icl_logL_rho <- function(
  R,
  priors,
  W
) {
  gammas <- priors[["gammas_0"]]
  N_W <- colSums(W)

  return(mvlbeta(N_W + gammas) - mvlbeta(gammas))
}


icl_logL_P_approx <- function(Pstar, Sigma, sigma2) {}


#' Fixed-point algorithm to compute \eqn{P^{*}} the maximum a
#' posteriori
#'
#' @param Pinit the inital value of P from which to start the
#' estimation
#' @inheritParams log_prior_plus_log_Z_given_P
#' @param tolerance the tolerance to which we stop the algorithm.
#' Defaults to 1e-3
#' @param max_iterations the maximum iterations of the
#' optimization. Defaults to 1000.
#'
#' @return estimated \eqn{P^{*}}
relaxed_fixed_point_MAP_P <- function(
  Pinit,
  Sigma,
  sigma2,
  Z,
  basis = default_Psi_function(ncol(Pinit) + 1),
  transformation = ilrInv,
  tolerance = 1e-3,
  max_iterations = 10000L,
  relax_coeff = 0.5
) {
  Pold <- Pinit
  Pold2 <- Pold
  probs_old <- transformation(Pold, basis = basis)
  Pnew <- relax_coeff *
    (sigma2 * Sigma %*% (Z - probs_old) %*% t(basis)) +
    (1 - relax_coeff) * Pold
  iter <- 1
  obj_vec <- c(log_prior_plus_log_Z_given_P(
    Pnew,
    Z,
    Sigma,
    sigma2,
    basis,
    transformation
  ))
  difference_vec <- c(norm(Pnew - Pold, type = "F"))
  grad_vec <- c(
    sum(log_prior_plus_log_Z_given_P_grad(
      Pnew,
      Z,
      Sigma,
      sigma2,
      basis,
      transformation
    ))^2
  )
  while (
    iter < max_iterations &&
      norm(Pnew - Pold, type = "F") > tolerance &&
      all(Pnew != Pold2)
  ) {
    Pold2 <- Pold
    Pold <- Pnew
    probs_old <- transformation(Pold, basis = basis)
    Pnew <- relax_coeff *
      (sigma2 * Sigma %*% (Z - probs_old) %*% t(basis)) +
      (1 - relax_coeff) * Pold
    obj_vec <- c(
      obj_vec,
      log_prior_plus_log_Z_given_P(
        Pnew,
        Z,
        Sigma,
        sigma2,
        basis,
        transformation
      )
    )
    difference_vec <- c(difference_vec, norm(Pnew - Pold, type = "F"))
    grad_vec <- c(
      grad_vec,
      sum(log_prior_plus_log_Z_given_P_grad(
        Pnew,
        Z,
        Sigma,
        sigma2,
        basis,
        transformation
      ))^2
    )
    iter <- iter + 1
  }
  return(list(P = Pnew, obj = obj_vec, grad = grad_vec, diff = difference_vec))
}

#' Unfinished attempt at Newton implementation
newton_MAP_P <- function(
  Pinit,
  Sigma,
  sigma2,
  Z,
  basis = default_Psi_function(ncol(Pinit) + 1),
  transformation = ilrInv,
  tolerance = 1e-3,
  max_iterations = 10000L
) {
  Pold <- Pinit
  Pold2 <- Pold
  # The Newton pace
  Pnew <- Pold -
    matrix(
      solve(log_prior_plus_log_Z_given_P_hessian_vec_row(
        Pold,
        Sigma,
        sigma2,
        basis,
        transformation
      )) %*%
        matrix(
          t(log_prior_plus_log_Z_given_P_grad(
            Pold,
            Z,
            Sigma,
            sigma2,
            basis,
            transformation
          )),
          ncol = 1
        ),
      byrow = TRUE,
      ncol = ncol(Pold)
    )
  iter <- 1
  obj_vec <- c(log_prior_plus_log_Z_given_P(
    Pnew,
    Z,
    Sigma,
    sigma2,
    basis,
    transformation
  ))
  difference_vec <- c(norm(Pnew - Pold, type = "F"))
  grad_vec <- c(
    sum(log_prior_plus_log_Z_given_P_grad(
      Pnew,
      Z,
      Sigma,
      sigma2,
      basis,
      transformation
    ))^2
  )
  while (
    iter < max_iterations &&
      norm(Pnew - Pold, type = "F") > tolerance &&
      all(Pnew != Pold2)
  ) {
    Pold2 <- Pold
    Pold <- Pnew
    Pnew <- Pold -
      matrix(
        solve(log_prior_plus_log_Z_given_P_hessian_vec_row(
          Pold,
          Sigma,
          sigma2,
          basis,
          transformation
        )) %*%
          matrix(
            t(log_prior_plus_log_Z_given_P_grad(
              Pold,
              Z,
              Sigma,
              sigma2,
              basis,
              transformation
            )),
            ncol = 1
          ),
        ncol = ncol(Pold)
      )
    obj_vec <- c(
      obj_vec,
      log_prior_plus_log_Z_given_P(
        Pnew,
        Z,
        Sigma,
        sigma2,
        basis,
        transformation
      )
    )
    difference_vec <- c(difference_vec, norm(Pnew - Pold, type = "F"))
    grad_vec <- c(
      grad_vec,
      sum(log_prior_plus_log_Z_given_P_grad(
        Pnew,
        Z,
        Sigma,
        sigma2,
        basis,
        transformation
      ))^2
    )
    iter <- iter + 1
  }
  return(list(P = Pnew, obj = obj_vec, grad = grad_vec, diff = difference_vec))
}
