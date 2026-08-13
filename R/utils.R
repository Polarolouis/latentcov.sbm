# Entropy, stolen from entropy package
entropy <- function(freqs, unit = c("log", "log2", "log10")) {
    unit <- match.arg(unit)

    freqs <- freqs / sum(freqs) # just to make sure ...

    H <- -sum(ifelse(freqs > 0, freqs * log(freqs), 0))

    if (unit == "log2") H <- H / log(2) # change from log to log2 scale
    if (unit == "log10") H <- H / log(10) # change from log to log10 scale

    return(H)
}

conditional_entropy <- function(tab, unit = c("log", "log2", "log10")) {
    entropy(tab, unit = unit) - entropy(colSums(tab), unit = unit)
}

mutual_information <- function(tab, unit = c("log", "log2", "log10")) {
    entropy(rowSums(tab), unit = unit) - conditional_entropy(tab, unit = unit)
}

#' A function to normalize each row of a given matrix
#'
#' @param mat the matrix for which the rows must be normalized
#' @param is_log a boolean indicating if the provided matrix is the
#' log of the unnormalized one. Default to TRUE.
#' @param tol the tolerance around which the values too close to 1 or
#' 0 are clamped to 1-tol and tol. Default to NULL, meaning no
#' clamping happens
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

#' A function
#' @param posterior_array A posterior array in the form Iteration x Chain x Parameter.
#' @param n the number of individuals, the rows of the outputted covariance matrix.
#' @param K the number of blocks, in the latent continuous space there are K-1 columns.
build_covariance_matrix <- function(posterior_array, n, K) {}


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

list_stan_to_chains_stan <- function(list_stan) {
    chains_stan_array <- simplify2array(c(list_stan))
    names(dimnames(chains_stan_array))[3] <- "Chain"
    dimnames(chains_stan_array)[["Chain"]] <- seq(dim(chains_stan_array)[3])
    return(chains_stan_array)
}

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

default_Psi_function <- function(K) {
    t(sapply(seq(K - 1), function(j) {
        sqrt((K - j) / (K - j + 1)) * c(rep(0, j - 1), -1, rep(1 / (K - j), K - j))
    }))
}

perm_matrix_from_order <- function(order) {
    sapply(order, function(i) {
        column <- rep(0, length(order))
        column[i] <- 1
        column
    })
}

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

onehot_encode <- function(fact, K = unique(fact)) {
    t(sapply(fact, function(label) as.integer(seq(K) == label)))
}

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
