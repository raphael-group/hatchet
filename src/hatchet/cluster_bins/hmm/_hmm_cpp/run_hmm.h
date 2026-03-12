#pragma once

#include <cstdint>
#include <vector>

/**
 * Return value of run_hmm_cpp.  All vectors are flat C-order arrays.
 */
struct RunHMMResult {
    std::vector<double> rdr_means;       // (K, M)
    std::vector<double> rdr_vars;        // (K, M)
    std::vector<double> baf_means;       // (K, M)
    std::vector<double> baf_taus;        // (M,)
    std::vector<double> log_startprobs;  // (K, 2)
    std::vector<double> posts;           // (N, K, 2)
    std::vector<double> lls0;            // (N, K)
    std::vector<double> lls1;            // (N, K)
    std::vector<double> elbo_trace;      // (n_iters+1)  — starts at -inf
    double loglik;
    int    n_iters_done;
};

/**
 * Full C++ EM loop for the 2-mixture BAF+RDR HMM.
 *
 * Runs n_iter Baum-Welch iterations in C++, eliminating per-iteration
 * Python/JAX dispatch overhead.  Parameters are copied from the initial
 * arrays and updated in place internally; results are returned in a
 * RunHMMResult struct.
 *
 * Array shapes (all C-contiguous float64 unless noted):
 *   X_rdrs        (N, M) — RDR observations
 *   X_alphas      (N, M) — A-allele counts (for loglik + BAF M-step + tau MLE)
 *   X_betas       (N, M) — B-allele counts (for loglik + BAF M-step + tau MLE)
 *   X_totals      (N, M) — total allele counts (for loglik)
 *   seg_lengths   (S,) int64 — number of bins per segment
 *   log_switchprobs (N,) — log phase-switch probability per bin
 *   log_stayprobs   (N,) — log phase-stay   probability per bin
 *   log_transmat  (K, K) — log cluster transition matrix (fixed)
 *   rdr_means0    (K, M) — initial RDR means
 *   rdr_vars0     (K, M) — initial RDR variances
 *   baf_means0    (K, M) — initial BAF means
 *   baf_taus0     (M,)   — initial BB dispersion
 */
RunHMMResult run_hmm_cpp(
    int K, int N, int M, int S,
    const double*  X_rdrs,
    const double*  X_alphas,
    const double*  X_betas,
    const double*  X_totals,
    const int64_t* seg_lengths,
    const double*  log_switchprobs,
    const double*  log_stayprobs,
    const double*  log_transmat,
    const double*  rdr_means0,
    const double*  rdr_vars0,
    const double*  baf_means0,
    const double*  baf_taus0,
    int    n_iter,
    double min_covar,
    double tol_ll,
    double tol,
    int    tau_iters,
    double min_tau,
    double max_tau,
    double baf_eps,
    double ig_alpha = 10.0,
    const double* ig_beta = nullptr
);
