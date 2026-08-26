#include "P_samplers.h"
#include "mvnorm.h"

// [[Rcpp::depends(RcppArmadillo)]]

//' Inverse pivot coordinate transformation (C++)
//'
//' Performs the inverse pivot coordinate transformation for the latent
//' covariance stochastic block model, mapping the \eqn{K-1} latent
//' coordinates of a row back to the \eqn{K}-simplex of membership
//' probabilities.
//' @param x Input matrix (\eqn{n \times K-1}) of latent coordinates
//' @param norm Normalization type ("orthogonal" or "orthonormal")
//' @param log if \code{TRUE}, the log-probabilities are returned
//' @return A matrix (\eqn{n \times K}) of (log) simplex probabilities
//'
//' @export
// [[Rcpp::export]]
arma::mat pivot_coord_inv(arma::mat &x, std::string norm = "orthonormal",
                          bool log = false) {
  // Mirror precisely the R implementation: x <- -x and then operate on that
  arma::mat xneg = -x;
  arma::mat xback;
  arma::mat y(x.n_rows, x.n_cols + 1, arma::fill::zeros);
  int D = x.n_cols + 1;
  double first_fill = 1.0;

  if (norm != "orthogonal" && norm != "orthonormal") {
    Rcpp::stop("Norm %s not implemented !", norm);
  }

  if (norm == "orthonormal") {
    first_fill = -std::sqrt((double)(D - 1) / (double)D);
  } else if (norm == "orthogonal") {
    first_fill = 1.0;
  }

  y.col(0) = first_fill * xneg.col(0);

  for (int i = 1; i < (int)y.n_cols; ++i) {
    for (int j = 0; j < i; ++j) {
      unsigned int ull_i = static_cast<unsigned int>(i);
      unsigned int ull_j = static_cast<unsigned int>(j);
      double denom = 1.0;
      if (norm == "orthonormal") {
        denom = std::sqrt((double)(D - ull_j) * (double)(D - ull_j - 1.0));
      }
      y.col(ull_i) += xneg.col(ull_j) / denom;
    }
  }

  for (int i = 1; i < (int)y.n_cols - 1; ++i) {
    unsigned int ull_i = static_cast<unsigned int>(i);
    double multip = 1.0;
    if (norm == "orthonormal") {
      multip = std::sqrt((double)(D - i - 1) / (double)(D - i));
    }
    y.col(ull_i) -= xneg.col(ull_i) * multip;
  }

  arma::vec max_rows = arma::max(y, 1);
  arma::mat yexp = arma::exp(y - arma::repmat(max_rows, 1, y.n_cols));
  xback = yexp / arma::repmat(arma::sum(yexp, 1), 1, y.n_cols);

  if (log) {
    xback = arma::log(xback);
  }

  return xback;
}

using namespace arma;
//' A Metropolis-Hastings sampler for P with a classical random walk (C++)
//'
//' Proposes each row \eqn{P_i} with an isotropic Gaussian random walk of
//' variance \code{rho} and accepts with the full Metropolis-Hastings ratio
//' combining the multinomial prior and the \eqn{N(\mu_i, \Sigma_i)}
//' conditional.
//' @param P a matrix of size \eqn{n_1 \times K-1} specifying latent
//'   position of the nodes probabilities of membership in the latent space
//' @param Z a matrix of size \eqn{n_1\times K} with a single 1 per line
//'   indicating the membership of row node \eqn{i}, \eqn{Z_{i,k} = 1} if
//'   \eqn{i} is in group \eqn{k} 0 else
//' @param Sigma a covariance matrix describing the covariance between row
//'   nodes \eqn{(n_1\times n_1)}
//' @param sigma2 a variance parameter indicating the variance between the
//'   K-1 columns of P
//' @param minibatch a boolean indicating wether to update the rows in the
//'   lexical order or to sample at each iteration. Default to TRUE
//' @param niter_metropolis an integer specifying the number of metropolis
//'   iterations to perform. Defaults to 50.
//' @param rho a double indicating the variance of the gaussian proposal
//'   distribution.
//' @return an updated matrix of latent positions, of the same size as
//'   \code{P}
//' @seealso [sample_P_metropolis_classical()],
//' [sample_P_metropolis_trick()],
//' [sample_P_metropolis_trick_cpp()]
//' @export
// [[Rcpp::export]]
arma::mat sample_P_metropolis_classical_cpp(arma::mat &P, arma::mat &Z,
                                            arma::mat &Sigma, double sigma2,
                                            bool minibatch = true,
                                            int niter_metropolis = 50,
                                            double rho = 1.0) {
  Rcpp::RNGScope scope;
  mat running_P = P;
  int n = running_P.n_rows;

  arma::rowvec rho_values = {1.0, 0.1, 10};

  arma::uvec row_order;

  if (minibatch) {
    row_order = arma::randperm(n);
  } else {
    row_order = arma::regspace<uvec>(0, n - 1);
  }

  // Loop over the individuals
  int accepted_count = 0;
  arma::uvec::iterator ind_end = row_order.end();
  for (arma::uvec::iterator ind_id = row_order.begin(); ind_id != ind_end;
       ++ind_id) {
    int indiv_idx = *ind_id;
    for (int iter_metro = 0; iter_metro < niter_metropolis; ++iter_metro) {
      // Sample the variance of the normal
      double rho_iter =
          rho * Rcpp::RcppArmadillo::sample(rho_values, 1, false)(0);

      arma::rowvec Pi_candidate = running_P.row(indiv_idx);
      arma::rowvec noise =
          arma::randn<arma::rowvec>(P.n_cols) * std::sqrt(rho_iter);

      Pi_candidate += noise;

      arma::rowvec Pi_old = running_P.row(indiv_idx);

      double log_u = log(Rcpp::runif(1)(0));

      arma::uvec Zi_vec = find(Z.row(indiv_idx) == 1, 1, "first");
      uint Zi = Zi_vec(0);

      if (DEBUG_SAMPLE) {
        Rcpp::Rcout << "Z_" << indiv_idx + 1 << " = " << Zi + 1 << std::endl;
      }

      arma::rowvec ind_mu =
          mean_of_Pi_given_P_min_i_sigma(running_P, Sigma, sigma2, indiv_idx);
      arma::mat ind_cov =
          cov_of_Pi_given_P_min_i_sigma(running_P, Sigma, sigma2, indiv_idx);

      double log_accept =
          pivot_coord_inv(Pi_candidate, "orthonormal", true)(Zi) -
          pivot_coord_inv(Pi_old, "orthonormal", true)(Zi);

      if (DEBUG_SAMPLE) {
        Rcpp::Rcout << "(multinom) accept_prob = " << exp(log_accept)
                    << std::endl;
      }

      log_accept += dmvnorm(Pi_candidate, ind_mu, ind_cov, true)(0) -
                    dmvnorm(Pi_old, ind_mu, ind_cov, true)(0);
      if (DEBUG_SAMPLE) {
        Rcpp::Rcout << "(multinom+mvnorm) accept_prob = " << exp(log_accept)
                    << " | u = " << exp(log_u) << std::endl;
      }

      if (log_u < log_accept) {
        if (DEBUG_SAMPLE) {
          Rcpp::Rcout << "Accepted !" << std::endl << "Old Pi = ";
          Pi_old.print(Rcpp::Rcout);
          Rcpp::Rcout << "New Pi = ";
          Pi_candidate.print(Rcpp::Rcout);
        }

        running_P.row(indiv_idx) = Pi_candidate;
        accepted_count++;
      }
    }
  }

  if (DEBUG_SAMPLE) {
    Rcpp::Rcout << "Mean accepted rate = "
                << (float)accepted_count / (float)(niter_metropolis * n)
                << std::endl;
  }

  return running_P;
};

//' A Metropolis-Hastings sampler for P with the clever proposition (C++)
//'
//' Proposes each row \eqn{P_i} from its full conditional
//' \eqn{N(\mu_i, \Sigma_i)} (as the R \code{sample_P_metropolis_trick}):
//' the Gaussian terms cancel in the acceptance ratio, leaving only the
//' multinomial prior ratio. This is the C++ version of
//' \code{sample_P_metropolis_trick}.
//' @param rho unused double for compatibility only
//' @inheritParams sample_P_metropolis_trick
//' @return an updated matrix of latent positions, of the same size
//' as \code{P}
//' @seealso [sample_P_metropolis_classical()],
//' [sample_P_metropolis_classical_cpp()],
//' [sample_P_metropolis_trick()]
//' @export
// [[Rcpp::export]]
arma::mat sample_P_metropolis_trick_cpp(arma::mat &P, arma::mat &Z,
                                        arma::mat &Sigma, double sigma2,
                                        bool minibatch = true,
                                        int niter_metropolis = 50,
                                        double rho = 1.0) {
  Rcpp::RNGScope scope;
  (void)rho;
  mat running_P = P;
  int n = running_P.n_rows;

  arma::uvec row_order;

  if (minibatch) {
    row_order = arma::randperm(n);
  } else {
    row_order = arma::regspace<uvec>(0, n - 1);
  }

  // Loop over the individuals
  int accepted_count = 0;
  arma::uvec::iterator ind_end = row_order.end();
  for (arma::uvec::iterator ind_id = row_order.begin(); ind_id != ind_end;
       ++ind_id) {
    int indiv_idx = *ind_id;
    for (int iter_metro = 0; iter_metro < niter_metropolis; ++iter_metro) {
      // Conditional distribution of Pi given P_minus_i
      arma::rowvec ind_mu =
          mean_of_Pi_given_P_min_i_sigma(running_P, Sigma, sigma2, indiv_idx);
      arma::mat ind_cov =
          cov_of_Pi_given_P_min_i_sigma(running_P, Sigma, sigma2, indiv_idx);

      // Proposal drawn from the conditional Gaussian (sample_Pi_given): the
      // mvnorm terms cancel in the acceptance ratio, leaving only the
      // multinomial prior ratio.
      arma::rowvec Pi_old = running_P.row(indiv_idx);
      arma::rowvec Pi_candidate = rmvnorm(1, ind_mu, ind_cov).row(0);

      double log_u = log(Rcpp::runif(1)(0));

      arma::uvec Zi_vec = find(Z.row(indiv_idx) == 1, 1, "first");
      uint Zi = Zi_vec(0);

      double log_accept =
          pivot_coord_inv(Pi_candidate, "orthonormal", true)(Zi) -
          pivot_coord_inv(Pi_old, "orthonormal", true)(Zi);

      if (DEBUG_SAMPLE) {
        Rcpp::Rcout << "(multinom) accept_prob = " << exp(log_accept)
                    << " | u = " << exp(log_u) << std::endl;
      }

      if (log_u < log_accept) {
        if (DEBUG_SAMPLE) {
          Rcpp::Rcout << "Accepted !" << std::endl << "Old Pi = ";
          Pi_old.print(Rcpp::Rcout);
          Rcpp::Rcout << "New Pi = ";
          Pi_candidate.print(Rcpp::Rcout);
        }

        running_P.row(indiv_idx) = Pi_candidate;
        accepted_count++;
      }
    }
  }

  if (DEBUG_SAMPLE) {
    Rcpp::Rcout << "Mean accepted rate = "
                << (float)accepted_count / (float)(niter_metropolis * n)
                << std::endl;
  }

  return running_P;
};

//' @inheritParams sample_P_metropolis_trick
//' @param block_size the size of blocks to update simultaneously (default: 1
// for individual updates) ' @return an updated matrix of latent positions, of
// the same size ' as \code{P} ' @seealso [sample_P_metropolis_classical()], '
//[sample_P_metropolis_classical_cpp()], ' [sample_P_metropolis_trick()] '
//@export
// [[Rcpp::export]]
arma::mat sample_P_metropolis_trick_cpp_block(arma::mat &P, arma::mat &Z,
                                              arma::mat &Sigma, double sigma2,
                                              bool minibatch = true,
                                              int niter_metropolis = 50,
                                              double rho = 1.0,
                                              int block_size = 1) {
  Rcpp::RNGScope scope;
  (void)rho;
  mat running_P = P;
  int n = running_P.n_rows;

  arma::uvec row_order;

  if (minibatch) {
    row_order = arma::randperm(n);
  } else {
    row_order = arma::regspace<uvec>(0, n - 1);
  }

  // Loop over the individuals in blocks
  int accepted_count = 0;
  int total_proposals = 0;

  // Process individuals in blocks
  arma::uvec::iterator ind_end = row_order.end();
  arma::uvec::iterator ind_id = row_order.begin();

  while (ind_id != ind_end) {
    // Determine the current block size (may be smaller for the last block)
    int current_block_size = std::min(block_size, (int)(ind_end - ind_id));

    // Collect indices for the current block
    arma::uvec block_indices(current_block_size);
    for (int i = 0; i < current_block_size; ++i) {
      block_indices(i) = *(ind_id + i);
    }

    // For each iteration of Metropolis-Hastings within the block
    for (int iter_metro = 0; iter_metro < niter_metropolis; ++iter_metro) {
      // Store old values for all individuals in the block
      arma::mat old_P_block = running_P.rows(block_indices);

      // Generate proposals for all individuals in the block
      arma::mat new_P_block(old_P_block.n_rows, old_P_block.n_cols);

      // For each individual in the block, compute conditional distribution and
      // propose
      for (int i = 0; i < current_block_size; ++i) {
        int indiv_idx = block_indices(i);

        // Conditional distribution of Pi given P_minus_i
        arma::rowvec ind_mu =
            mean_of_Pi_given_P_min_i_sigma(running_P, Sigma, sigma2, indiv_idx);
        arma::mat ind_cov =
            cov_of_Pi_given_P_min_i_sigma(running_P, Sigma, sigma2, indiv_idx);

        // Proposal drawn from the conditional Gaussian
        arma::rowvec Pi_old = running_P.row(indiv_idx);
        arma::rowvec Pi_candidate = rmvnorm(1, ind_mu, ind_cov).row(0);
        new_P_block.row(i) = Pi_candidate;
      }

      // Compute acceptance probability for the entire block
      double log_accept_total = 0.0;
      for (int i = 0; i < current_block_size; ++i) {
        int indiv_idx = block_indices(i);
        arma::rowvec Pi_old = old_P_block.row(i);
        arma::rowvec Pi_candidate = new_P_block.row(i);

        arma::uvec Zi_vec = find(Z.row(indiv_idx) == 1, 1, "first");
        uint Zi = Zi_vec(0);

        double log_accept_single =
            pivot_coord_inv(Pi_candidate, "orthonormal", true)(Zi) -
            pivot_coord_inv(Pi_old, "orthonormal", true)(Zi);

        log_accept_total += log_accept_single;
      }

      double log_u = log(Rcpp::runif(1)(0));
      total_proposals += current_block_size;

      if (log_u < log_accept_total) {
        // Accept the entire block
        running_P.rows(block_indices) = new_P_block;
        accepted_count += current_block_size;
      }
    }

    // Move to the next block
    std::advance(ind_id, current_block_size);
  }

  if (DEBUG_SAMPLE) {
    Rcpp::Rcout << "Mean accepted rate = "
                << (float)accepted_count / (float)total_proposals << std::endl;
  }

  return running_P;
};

//' Mean of Pi given P minus i and Sigma (C++)
//'
//' Computes the mean of the conditional distribution of a single row of P
//' given the rest of P, Sigma, and sigma2. This is the C++ version of
//' \code{cond_Pi_given_P_min_i_sigma}.
//' @param P a matrix of size \eqn{n_1 \times K-1} specifying latent
//'   positions
//' @param Sigma a covariance matrix between the row nodes
//' @param sigma2 variance parameter indicating the variance between the
//'   K-1 columns of P
//' @param i 0-based index of the row to compute the conditional for
//' @return the conditional mean, a row vector of length \eqn{K-1}
// [[Rcpp::export]]
arma::rowvec mean_of_Pi_given_P_min_i_sigma(arma::mat &P, arma::mat &Sigma,
                                            double sigma2, int i) {
  int n = P.n_rows;

  // Indices minus_i = setdiff(1:n, i)  (0-based: tous sauf i)
  arma::uvec minus_i = arma::regspace<arma::uvec>(0, n - 1);
  minus_i.shed_row(i);

  // Si = Sigma[i, minus_i] %*% solve(Sigma[minus_i, minus_i])
  arma::rowvec Sigma_i_minus = Sigma.row(i);
  Sigma_i_minus.shed_col(i); // Sigma[i, minus_i]
  arma::mat Sigma_minus_minus = Sigma.submat(minus_i, minus_i);
  // Resolve the linear system directly instead of computing the explicit
  // inverse: more stable numerically and closer to R's solve().
  arma::vec Si_t = arma::solve(Sigma_minus_minus, Sigma_i_minus.t());
  arma::rowvec Si = Si_t.t();

  // mean_i = Si %*% P[minus_i, ]
  arma::rowvec mean_i = Si * P.rows(minus_i);
  return mean_i;
}

//' Covariance of Pi given P minus i and Sigma (C++)
//'
//' Computes the covariance of the conditional distribution of a single row
//' of P given the rest of P, Sigma, and sigma2. This is the C++ version of
//' \code{cond_Pi_given_P_min_i_sigma}.
//' @param P a matrix of size \eqn{n_1 \times K-1} specifying latent
//'   positions
//' @param Sigma a covariance matrix between the row nodes
//' @param sigma2 variance parameter indicating the variance between the
//'   K-1 columns of P
//' @param i 0-based index of the row to compute the conditional for
//' @return the conditional covariance, a matrix of size
//'   \eqn{K-1 \times K-1}
// [[Rcpp::export]]
arma::mat cov_of_Pi_given_P_min_i_sigma(arma::mat &P, arma::mat &Sigma,
                                        double sigma2, int i) {
  int n = P.n_rows;
  int K_minus_1 = P.n_cols;

  // Indices minus_i = setdiff(1:n, i)  (0-based: tous sauf i)
  arma::uvec minus_i = arma::regspace<arma::uvec>(0, n - 1);
  minus_i.shed_row(i);

  // Si = Sigma[i, minus_i] %*% solve(Sigma[minus_i, minus_i])
  arma::rowvec Sigma_i_minus = Sigma.row(i);
  Sigma_i_minus.shed_col(i); // Sigma[i, minus_i]
  arma::mat Sigma_minus_minus = Sigma.submat(minus_i, minus_i);
  // Resolve the linear system directly instead of computing the explicit
  // inverse: more stable numerically and closer to R's solve().
  arma::vec Si_t = arma::solve(Sigma_minus_minus, Sigma_i_minus.t());
  arma::rowvec Si = Si_t.t();

  // cov_i = sigma2 * (Sigma[i,i] - Si %*% Sigma[-i, i]) * I
  arma::vec Sigma_minus_i = Sigma.col(i);
  Sigma_minus_i.shed_row(i); // Sigma[-i, i]
  double scalar = Sigma(i, i) - arma::as_scalar(Si * Sigma_minus_i);
  arma::mat cov_i =
      (sigma2 * scalar) * arma::eye<arma::mat>(K_minus_1, K_minus_1);

  return cov_i;
}

//' Mean of P_b given P minus b and Sigma (C++)
//'
//' Computes the mean of the conditional distribution of a block of P
//' given the rest of P, Sigma, and sigma2. This is the C++ version of
//' \code{cond_Pi_given_P_min_i_sigma}.
//' @param P a matrix of size \eqn{n_1 \times K-1} specifying latent
//'   positions
//' @param Sigma a covariance matrix between the row nodes
//' @param sigma2 variance parameter indicating the variance between the
//'   K-1 columns of P
//' @param indiv_indicees 0-based indices of the rows to compute the conditional
// for ' @return the conditional mean, a row vector of length \eqn{K-1}
// [[Rcpp::export]]
arma::mat block_mean_of_Pi_given_P_min_i_sigma(arma::mat &P, arma::mat &Sigma,
                                               double sigma2,
                                               arma::uvec indiv_indices) {
  int n = P.n_rows;

  // Indices minus_i = setdiff(1:n, i)  (0-based: tous sauf i)
  arma::uvec minus_i = arma::regspace<arma::uvec>(0, n - 1);
  minus_i.shed_rows(indiv_indices);

  // Si = Sigma[i, minus_i] %*% solve(Sigma[minus_i, minus_i])
  arma::mat Sigma_i_minus = Sigma.rows(indiv_indices);
  Sigma_i_minus.shed_cols(indiv_indices); // Sigma[i, minus_i]
  arma::mat Sigma_minus_minus = Sigma.submat(minus_i, minus_i);
  // Resolve the linear system directly instead of computing the explicit
  // inverse: more stable numerically and closer to R's solve().
  arma::mat Si_t = arma::solve(Sigma_minus_minus, Sigma_i_minus.t());
  arma::mat Si = Si_t.t();

  // mean_i = Si %*% P[minus_i, ]
  arma::mat mean_i = Si * P.rows(minus_i);
  return mean_i;
}
