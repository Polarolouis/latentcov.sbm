#' Calculate entropy of a probability distribution
#'
#' @param freqs A vector of frequencies or probabilities
#' @param unit The logarithm base to use for entropy calculation. Options are "log", "log2", or "log10"
#' @return The entropy value
#' @export
entropy <- function(freqs, unit = c("log", "log2", "log10")) {
  unit <- match.arg(unit)

  freqs <- freqs / sum(freqs) # just to make sure ...

  H <- -sum(ifelse(freqs > 0, freqs * log(freqs), 0))

  if (unit == "log2") H <- H / log(2) # change from log to log2 scale
  if (unit == "log10") H <- H / log(10) # change from log to log10 scale

  return(H)
}

#' Calculate conditional entropy
#'
#' @param tab A contingency table (matrix)
#' @param unit The logarithm base to use for entropy calculation. Options are "log", "log2", or "log10"
#' @return The conditional entropy value
#' @export
conditional_entropy <- function(tab, unit = c("log", "log2", "log10")) {
  entropy(tab, unit = unit) - entropy(colSums(tab), unit = unit)
}

#' Calculate mutual information
#'
#' @param tab A contingency table (matrix)
#' @param unit The logarithm base to use for entropy calculation. Options are "log", "log2", or "log10"
#' @return The mutual information value
#' @export
mutual_information <- function(tab, unit = c("log", "log2", "log10")) {
  entropy(rowSums(tab), unit = unit) - conditional_entropy(tab, unit = unit)
}

#' Normalize each row of a matrix
#'
#' This function normalizes the rows of a matrix so that each row sums to 1.
#' It can handle both log-scale and regular-scale matrices.
#'
#' @param mat The matrix to normalize
#' @param is_log A boolean indicating if the provided matrix is the log of the unnormalized one. Default is TRUE.
#' @param tol The tolerance around which values too close to 1 or 0 are clamped to 1-tol and tol. Default is NULL, meaning no clamping happens.
#' @return The normalized matrix
#' @export
row_normalize_matrix <- function(mat, is_log = TRUE, tol = NULL) {
  if (!is_log) {
    mat <- log(mat)
  }
  mi <- apply(mat, 1, max)
  mat_exp <- exp(mat - mi)
  normalized_mat <- mat_exp / apply(mat_exp, 1, sum)

  if (!is.null(tol)) {
    normalized_mat <- pmin(pmax(normalized_mat, tol), 1 - tol)
  }

  return(normalized_mat)
}

#' Build a covariance matrix from posterior samples
#'
#' @param posterior_array A posterior array in the form Iteration x Chain x Parameter.
#' @param n The number of individuals, the rows of the outputted covariance matrix.
#' @param K The number of blocks, in the latent continuous space there are K-1 columns.
#' @return A covariance matrix
#' @export
build_covariance_matrix <- function(posterior_array, n, K) {
  # Function implementation would go here
}


#' Flatten a posterior array into a stan-style parameter matrix
#'
#' Reshapes a multi-dimensional posterior array (whose first dimension is
#' the iterations) into a matrix of dimensions \eqn{\text{niter} \times
#' \text{nparams}}, building lookalike parameter names such as
#' \code{"name[2,3]"} from the trailing dimensions.
#'
#' @param x A multi-dimensional array with iterations in the first dimension.
#' @param name The base (prefix) parameter name.
#'
#' @return A matrix with dimensions named \code{Iteration} and \code{Parameter}.
#'
#' @seealso [list_arrays_to_stan()], [list_stan_to_chains_stan()]
#' @keywords internal
stan_flatten <- function(x, name) {
  d <- dim(x)
  niter <- d[1]
  idx <- expand.grid(lapply(d[-1], seq_len))
  if (nrow(idx) > 1) {
    param_names <- apply(idx, 1, function(i) {
      paste0(name, paste0("[", paste0(i, collapse = ","), "]", collapse = ""))
    })
  } else {
    param_names <- name
  }
  out <- matrix(x, nrow = niter)
  dimnames(out) <- list(
    Iteration = seq_len(niter),
    Parameter = param_names
  )
  out
}

#' Convert a list of posterior arrays into a stan-style parameter matrix
#'
#' Combines several posterior arrays (typically the outputs of the model
#' samplers such as [gibbs_sampling_lbm_poisson()]) into a single matrix
#' whose columns are the flattened parameters. The trailing \code{"_array"}
#' suffix of each element name is stripped before building the parameter
#' names.
#'
#' @param array_list A named list of posterior arrays sharing the same
#'   number of iterations.
#'
#' @return A matrix with dimensions named \code{Iteration} and
#'   \code{Parameter}.
#'
#' @seealso [gibbs_sampling_lbm_poisson()], [gibbs_sampling_lbm_cov_poisson()],
#'   [lbm_results_to_stan_draws()]
#' @keywords internal
list_arrays_to_stan <- function(array_list) {
  out <- do.call("cbind", sapply(seq_along(array_list), function(idx) {
    param_name <- names(array_list)[idx]
    param_name <- substr(param_name, 1, nchar(param_name) - 6)
    if (!all(is.na(array_list[[idx]]))) {
      return(stan_flatten(array_list[[idx]], param_name))
    } else {
      return(NULL)
    }
  }))
  names(dimnames(out)) <- c("Iteration", "Parameter")
  return(out)
}

#' Stack a list of stan-style matrices into a chains array
#'
#' Takes a list of \eqn{\text{iteration} \times \text{parameter}} matrices
#' (one per MCMC chain) and stacks them along a new third dimension
#' labelled \code{Chain}.
#'
#' @param list_stan A list of matrices with identical dimensions, each
#'   produced by [list_arrays_to_stan()].
#'
#' @return A 3-dimensional array of dimensions
#'   \eqn{\text{iteration} \times \text{parameter} \times \text{chain}}.
#'
#' @seealso [list_arrays_to_stan()], [lbm_results_to_stan_draws()]
#' @keywords internal
list_stan_to_chains_stan <- function(list_stan) {
  chains_stan_array <- simplify2array(c(list_stan))
  names(dimnames(chains_stan_array))[3] <- "Chain"
  dimnames(chains_stan_array)[["Chain"]] <- seq(dim(chains_stan_array)[3])
  return(chains_stan_array)
}

#' Find row and column label permutations between two connectivity matrices
#'
#' Estimates the permutations of the row and column block labels that align
#' \code{alpha_to_align} onto the reference matrix \code{alpha_ref}. Rows are
#' matched first with an assignment problem solved on the pairwise Euclidean
#' distances ([clue::solve_LSAP]); the columns are then matched after applying
#' the row permutation.
#'
#' @param alpha_ref The reference connectivity matrix (\eqn{K \times R}).
#' @param alpha_to_align The connectivity matrix to permute.
#'
#' @return A list with two elements: `row_perm` and `col_perm`, the integer
#'   permutation vectors to apply to the rows and columns of
#'   \code{alpha_to_align} respectively.
#'
#' @seealso [find_permutation_alphas_L2()] for an exact (but exponential)
#'   alternative, [delabel_switch_stan()] for its use in label-switching
#'   correction.
#' @importFrom stats dist
find_permutation_alphas <- function(alpha_ref, alpha_to_align) {
  ## Match rows
  row_cost <- as.matrix(dist(rbind(alpha_ref, alpha_to_align)))[
    seq_len(nrow(alpha_ref)),
    nrow(alpha_ref) + seq_len(nrow(alpha_to_align))
  ]

  row_perm <- clue::solve_LSAP(row_cost)

  alpha2 <- alpha_to_align[row_perm, , drop = FALSE]

  ## Match columns
  col_cost <- as.matrix(dist(rbind(t(alpha_ref), t(alpha2))))[
    seq_len(ncol(alpha_ref)),
    ncol(alpha_ref) + seq_len(ncol(alpha2))
  ]

  col_perm <- clue::solve_LSAP(col_cost)

  return(list(
    row_perm = as.vector(row_perm),
    col_perm = as.vector(col_perm)
  ))
}

#' Find label permutations by exhaustive L2 minimization
#'
#' Exactly minimises the sum of squared differences between \code{alpha_ref}
#' and \code{alpha_to_align} over all row and column permutations. Exhaustive,
#' hence only suited for small matrices.
#'
#' @inheritParams find_permutation_alphas
#'
#' @return A list with two elements: `row_perm` and `col_perm`, the integer
#'   permutation vectors to apply to the rows and columns of
#'   \code{alpha_to_align} respectively.
#'
#' @seealso [find_permutation_alphas()] for the fast (but approximate)
#'   alternative, [delabel_switch_stan()] for its use in label-switching
#'   correction.
find_permutation_alphas_L2 <- function(alpha_ref, alpha_to_align) {
  all_row_perms <- gtools::permutations(nrow(alpha_to_align), nrow(alpha_to_align))
  all_col_perms <- gtools::permutations(ncol(alpha_to_align), ncol(alpha_to_align))
  all_perms_loss <- outer(seq_len(nrow(all_row_perms)), seq_len(nrow(all_col_perms)), FUN = Vectorize(function(i_row, j_col) {
    sum((alpha_ref - alpha_to_align[all_row_perms[i_row, ], all_col_perms[j_col, ]])^2)
  }))
  best_perm <- which(all_perms_loss == min(all_perms_loss), arr.ind = TRUE)
  row_perm <- all_row_perms[best_perm[1], ]
  col_perm <- all_col_perms[best_perm[2], ]

  return(list(
    row_perm = as.vector(row_perm),
    col_perm = as.vector(col_perm)
  ))
}

#' Default Psi contrast matrix
#'
#' Builds the default \eqn{(K-1) \times K} contrast (basis) matrix used by the
#' pivot coordinate (isometric log-ratio style) transformation of the block
#' proportions.
#'
#' @param K Number of groups (blocks).
#'
#' @return A matrix of size \eqn{(K-1)\times K}.
default_Psi_function <- function(K) {
  t(sapply(seq(K - 1), function(j) {
    sqrt((K - j) / (K - j + 1)) * c(rep(0, j - 1), -1, rep(1 / (K - j), K - j))
  }))
}

#' Build a permutation matrix from an ordering
#'
#' @param order An integer vector giving the new order of a set of elements
#'   (a permutation of \code{seq_along(order)}).
#'
#' @return A permutation matrix \eqn{P} such that \eqn{x[order] = Px}.
#' @keywords internal
perm_matrix_from_order <- function(order) {
  sapply(order, function(i) {
    column <- rep(0, length(order))
    column[i] <- 1
    column
  })
}

#' Correct label switching across MCMC chains
#'
#' Relabels the block indices of every MCMC chain (except the first, used as
#' reference) so that all chains share a common labelling, alleviating the
#' label-switching problem. Permutations are estimated on the posterior means
#' of the connectivity matrix `alpha`, then propagated to `rho`, `pi`, `P`,
#' `Z` and `W` parameters.
#'
#' @param draws A [posterior::draws_array] object with dimensions
#'   iteration \eqn{\times} chain \eqn{\times} parameter.
#' @param K Number of row groups.
#' @param R Number of column groups.
#' @param Psi_function A function building the contrast matrix used to
#'   rotate the `P` coordinates. Defaults to [default_Psi_function()].
#' @param find_permutations A function estimating the row and column
#'   permutations between two alpha matrices, e.g.
#'   [find_permutation_alphas()] or [find_permutation_alphas_L2()].
#'
#' @return The relabelled draws array, with the same structure as `draws`.
#'
#' @seealso [lbm_results_to_stan_draws()], [find_permutation_alphas()]
#' @importFrom utils head tail
delabel_switch_stan <- function(draws, K, R, Psi_function = default_Psi_function, find_permutations = find_permutation_alphas_L2) {
  stopifnot("There must be at least two chains" = dim(draws)[2] > 1)

  var_idx_alphas <- which(startsWith(dimnames(draws)[[3]], "alpha"))
  start_alpha_var <- head(var_idx_alphas, 1)
  end_alpha_var <- tail(var_idx_alphas, 1)

  mean_alphas_array <- apply(draws[, , var_idx_alphas], 2:3, mean)

  alpha_matrices <- lapply(seq_len(nrow(mean_alphas_array)), function(row) matrix(mean_alphas_array[row, ], nrow = K, ncol = R))
  alpha_ref <- alpha_matrices[[1]]
  alpha_matrices <- alpha_matrices[-1]

  permutations_list <- lapply(alpha_matrices, find_permutations, alpha_ref = alpha_ref)

  # Apply permutations to relabel all chains
  draws_delabeled <- draws
  Psi <- Psi_function(K)

  # Process each chain starting from chain 2 (chain 1 is the reference)
  for (chain_idx in 2:dim(draws)[2]) {
    perm <- permutations_list[[chain_idx - 1]]
    row_perm <- perm$row_perm
    col_perm <- perm$col_perm

    # Apply row permutation to alpha
    for (iter in seq(posterior::niterations(draws))) {
      alpha_matrix <- matrix(draws[iter, chain_idx, start_alpha_var:end_alpha_var], nrow = K, ncol = R)
      alpha_matrix <- alpha_matrix[row_perm, col_perm, drop = FALSE]
      draws_delabeled[iter, chain_idx, start_alpha_var:end_alpha_var] <- as.vector(alpha_matrix)
    }

    # Apply column permutation to other block-indexed parameters

    ## rho

    rho_idx <- which(startsWith(dimnames(draws_delabeled)[[3]], "rho"))
    draws_delabeled[, chain_idx, rho_idx] <- draws[, chain_idx, rho_idx][, , col_perm, drop = FALSE]

    ##  pi (if they exists)
    pi_idx <- which(startsWith(dimnames(draws_delabeled)[[3]], "pi"))
    if (length(pi_idx) > 0) {
      draws_delabeled[, chain_idx, pi_idx] <- draws[, chain_idx, pi_idx][, , row_perm, drop = FALSE]
    }

    ## P
    P_idx <- which(startsWith(dimnames(draws_delabeled)[[3]], "P"))

    Perm_mat <- perm_matrix_from_order(row_perm)
    T_c <- round(Psi %*% Perm_mat %*% t(Psi), digits = 10)
    nind <- length(P_idx) / (K - 1)
    for (it in seq_len(posterior::niterations(draws_delabeled))) {
      # reconstruction de la matrice P
      Pmat <- matrix(
        draws[it, chain_idx, P_idx],
        nrow = nind,
        ncol = K - 1,
        byrow = FALSE
      )

      # permutation
      Pmat <- Pmat %*% t(T_c)

      # remise dans le draws_array
      draws_delabeled[it, chain_idx, P_idx] <- c(Pmat)
    }

    # Z
    Z_idx <- which(startsWith(dimnames(draws_delabeled)[[3]], "Z"))
    chain_Z <- draws[, chain_idx, Z_idx]
    draws_delabeled[, chain_idx, Z_idx] <- array(
      row_perm[chain_Z],
      dim = dim(chain_Z)
    )

    # W
    W_idx <- which(startsWith(dimnames(draws_delabeled)[[3]], "W"))
    chain_W <- draws[, chain_idx, W_idx]
    draws_delabeled[, chain_idx, W_idx] <- array(
      col_perm[chain_W],
      dim = dim(chain_W)
    )
  }

  return(draws_delabeled)
}

#' @describeIn delabel_switch_stan Alias kept for backward compatibility.
#' @keywords internal
delabel_switch_stan_per_iteration <- function(draws, K, R, Psi_function = default_Psi_function, find_permutations = find_permutation_alphas_L2) {
  stopifnot("There must be at least two chains" = dim(draws)[2] > 1)

  var_idx_alphas <- which(startsWith(dimnames(draws)[[3]], "alpha"))
  start_alpha_var <- head(var_idx_alphas, 1)
  end_alpha_var <- tail(var_idx_alphas, 1)

  mean_alphas_array <- apply(draws[, , var_idx_alphas], 2:3, mean)

  alpha_matrices <- lapply(seq_len(nrow(mean_alphas_array)), function(row) matrix(mean_alphas_array[row, ], nrow = K, ncol = R))
  alpha_ref <- alpha_matrices[[1]]
  alpha_matrices <- alpha_matrices[-1]

  permutations_list <- lapply(alpha_matrices, find_permutations, alpha_ref = alpha_ref)

  # Apply permutations to relabel all chains
  draws_delabeled <- draws
  Psi <- Psi_function(K)

  # Process each chain starting from chain 2 (chain 1 is the reference)
  for (chain_idx in 2:dim(draws)[2]) {
    perm <- permutations_list[[chain_idx - 1]]
    row_perm <- perm$row_perm
    col_perm <- perm$col_perm

    # Apply row permutation to alpha
    for (iter in seq(posterior::niterations(draws))) {
      alpha_matrix <- matrix(draws[iter, chain_idx, start_alpha_var:end_alpha_var], nrow = K, ncol = R)
      alpha_matrix <- alpha_matrix[row_perm, col_perm, drop = FALSE]
      draws_delabeled[iter, chain_idx, start_alpha_var:end_alpha_var] <- as.vector(alpha_matrix)
    }

    # Apply column permutation to other block-indexed parameters

    ## rho

    rho_idx <- which(startsWith(dimnames(draws_delabeled)[[3]], "rho"))
    draws_delabeled[, chain_idx, rho_idx] <- draws[, chain_idx, rho_idx][, , col_perm, drop = FALSE]

    ##  pi (if they exists)
    pi_idx <- which(startsWith(dimnames(draws_delabeled)[[3]], "pi"))
    if (length(pi_idx) > 0) {
      draws_delabeled[, chain_idx, pi_idx] <- draws[, chain_idx, pi_idx][, , row_perm, drop = FALSE]
    }

    ## P
    P_idx <- which(startsWith(dimnames(draws_delabeled)[[3]], "P"))

    Perm_mat <- perm_matrix_from_order(row_perm)
    T_c <- round(Psi %*% Perm_mat %*% t(Psi), digits = 10)
    nind <- length(P_idx) / (K - 1)
    for (it in seq_len(posterior::niterations(draws_delabeled))) {
      # reconstruction de la matrice P
      Pmat <- matrix(
        draws[it, chain_idx, P_idx],
        nrow = nind,
        ncol = K - 1,
        byrow = FALSE
      )

      # permutation
      Pmat <- Pmat %*% t(T_c)

      # remise dans le draws_array
      draws_delabeled[it, chain_idx, P_idx] <- c(Pmat)
    }

    # Z
    Z_idx <- which(startsWith(dimnames(draws_delabeled)[[3]], "Z"))
    chain_Z <- draws[, chain_idx, Z_idx]
    draws_delabeled[, chain_idx, Z_idx] <- array(
      row_perm[chain_Z],
      dim = dim(chain_Z)
    )

    # W
    W_idx <- which(startsWith(dimnames(draws_delabeled)[[3]], "W"))
    chain_W <- draws[, chain_idx, W_idx]
    draws_delabeled[, chain_idx, W_idx] <- array(
      col_perm[chain_W],
      dim = dim(chain_W)
    )
  }

  return(draws_delabeled)
}


#' Convert LBM Gibbs sampling results into stan-style posterior draws
#'
#' Gathers the outputs of one or several runs of the model Gibbs samplers,
#' flattens them into a [posterior::draws_array] object and, when several
#' chains are provided, applies the [delabel_switch_stan()] label-switching
#' correction. An optional burn-in/thinning can be applied.
#'
#' @param lbm_results The output of [gibbs_sampling_lbm_poisson()] or
#'   [gibbs_sampling_lbm_cov_poisson()]. Either a single named list (one
#'   chain) or a list of such lists (several chains).
#' @param K Number of row groups.
#' @param R Number of column groups.
#' @param Psi_function A function building the contrast matrix. Defaults to
#'   [default_Psi_function()].
#' @param find_permutation A function estimating the permutations between
#'   two alpha matrices. Defaults to [find_permutation_alphas()].
#' @param apply_burnin_thinning Logical; if `TRUE` the first half of the
#'   iterations are discarded as burn-in and the remaining ones are thinned
#'   by a factor of 10. Defaults to `FALSE`.
#'
#' @return A [posterior::draws_array] with dimensions
#'   iteration \eqn{\times} chain \eqn{\times} parameter.
#'
#' @seealso [list_arrays_to_stan()], [list_stan_to_chains_stan()],
#'   [delabel_switch_stan()]
lbm_results_to_stan_draws <- function(lbm_results, K, R, Psi_function = default_Psi_function, find_permutation = find_permutation_alphas, apply_burnin_thinning = FALSE) {
  if (!is.null(names(lbm_results))) {
    message("Only one chain provided !")
    multiple_lbm_results <- list(lbm_results)
  } else {
    multiple_lbm_results <- lbm_results
  }

  list_stan <- lapply(multiple_lbm_results, list_arrays_to_stan)

  stan_results <- list_stan_to_chains_stan(list_stan)

  draws <- posterior::as_draws_array(aperm(stan_results, c(1, 3, 2)))
  if (is.null(names(lbm_results))) {
    message("Delabel switching")
    draws <- delabel_switch_stan(draws, K = K, R = R, Psi_function, find_permutation)
  }
  if (apply_burnin_thinning) {
    burnin <- floor(dim(draws)[1] / 2)
    thinning <- 10L

    draws <- draws[seq(burnin + 1, dim(draws)[1], by = thinning), , ]
  }
  return(draws)
}

#' Compute the Adjusted Rand Index across chains and iterations
#'
#' Computes the Adjusted Rand Index (ARI) between the block labels sampled at
#' each iteration and chain and the true partition.
#'
#' @param draws A [posterior::draws_array] containing sampled block labels.
#' @param true The true partition. Either a vector of labels or a one-hot
#'   indicator matrix.
#' @param label The name of the label variable to extract from `draws`.
#'   Defaults to `"Z"`.
#'
#' @return An array (iteration \eqn{\times} chain) of ARI values.
#'
#' @seealso [onehot_encode()]
ARI_table <- function(draws, true, label = "Z") {
  label_draws <- posterior::subset_draws(draws, variable = label)

  if (is.matrix(true)) {
    true <- .rev_one_hot(true)
  }

  aritable <- t(sapply(seq(posterior::niterations(label_draws)), function(iter) {
    sapply(seq(posterior::nchains(label_draws)), function(chain_idx) {
      aricode::ARI(as.vector(label_draws[iter, chain_idx, ]), true)
    })
  }))

  aritable <- array(aritable, dim(aritable), dimnames = list("Iteration" = seq(posterior::niterations(draws)), "Chains" = seq(posterior::nchains(draws))))
  return(aritable)
}

#' One-hot encode a vector of labels
#'
#' @param fact A vector of group labels.
#' @param K The number of groups. Defaults to the unique values of `fact`.
#'
#' @return A matrix (\eqn{n \times K}) with a single 1 per row indicating the
#'   membership of each individual.
#'
#' @seealso [ARI_table()], `.one_hot`
onehot_encode <- function(fact, K = unique(fact)) {
  t(sapply(fact, function(label) as.integer(seq(K) == label)))
}

#' Check the identifiability of a Poisson LBM configuration
#'
#' Verifies, in the sense of Keribin et al., that a Poisson latent block
#' model configuration is identifiable: uniqueness of `alpha %*% rho` and
#' `t(pi) %*% alpha`, and enough row and column nodes. Aborts with an error
#' message if the configuration is not identifiable, otherwise prints a
#' success message.
#'
#' @param netMat The observed network (count) matrix.
#' @param alpha The connectivity matrix (\eqn{K \times R}).
#' @param pi The row block proportions (vector of size \eqn{K}).
#' @param rho The column block proportions (vector of size \eqn{R}).
#' @param K Number of row groups.
#' @param R Number of column groups.
#'
#' @return Invisibly `NULL`; used for its side effects (messages).
check_lbm_identifiability <- function(netMat, alpha, pi, rho, K, R) {
  # From Keribin et al
  taus <- as.vector(alpha %*% rho)
  if (any(duplicated(taus))) {
    cli::cli_abort(c("x" = "All elements of `alpha%*%rho` should be uniques !", "i" = "{cli::qty(length(unique(taus)))}The only unique value{?s} {?no/is/are} {.val {unique(taus)}}"))
  }
  sigmas <- as.vector(t(pi) %*% alpha)
  if (any(duplicated(sigmas))) {
    cli::cli_abort(c("x" = "All elements of `t(pi)%*%alpha` should be uniques !", "i" = "{cli::qty(length(unique(sigmas)))}The only unique value{?s} {?no/is/are} {.val {unique(sigmas)}}"))
  }

  if (nrow(netMat) < 2 * R - 1) {
    cli::cli_abort(c("x" = "There are not enough row nodes for the network to be identifiable !", "i" = "There should be at least {.val {2*R-1}} but there is only {.val {nrow(netMat)}}"))
  }

  if (ncol(netMat) < 2 * K - 1) {
    cli::cli_abort(c("x" = "There are not enough column nodes for the network to be identifiable !", "i" = "There should be at least {.val {2*K-1}} but there is only {.val {ncol(netMat)}}"))
  }

  cli::cli_alert_success("This configuration is identifiable in the sense of Keribin et al. !")
}

#' One-hot encode a vector of labels (internal)
#'
#' @param x A vector of group labels (integers from 1 to `Q`).
#' @param Q The number of groups.
#'
#' @return A matrix (\eqn{\text{length}(x) \times Q}) with a single 1 per row.
#' @keywords internal
.one_hot <- function(x, Q) {
  O <- matrix(0, length(x), Q)
  O[cbind(seq.int(length(x)), x)] <- 1
  return(O)
}

#' Convert a one-hot matrix back to labels (internal)
#'
#' @param X A one-hot indicator matrix.
#'
#' @return A vector of labels (the column index of the 1 of each row).
#' @keywords internal
.rev_one_hot <- function(X) {
  return(as.vector(max.col(X)))
}
