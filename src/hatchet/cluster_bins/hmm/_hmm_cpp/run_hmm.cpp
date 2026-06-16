/**
 * Full C++ EM loop for the 2-mixture BAF+RDR HMM.
 *
 * Calls, in each iteration:
 *   1. compute_loglik_cpp    — OpenMP-parallelised loglik kernel
 *   2. forward_backward_cpp  — OpenMP-parallelised fwd-bwd over segments
 *   3. update_start_probs    — segment-start posterior aggregation
 *   4. update_rdr_params_cpp — closed-form Gaussian M-step
 *   5. update_baf_tau_cpp    — (first tau_iters iters only) Brent tau MLE
 *   6. update_baf_means_cpp  — Brent BAF M-step per (k, m), uses updated tau
 *
 * Note: mhBAF folding (flipping BAF means > 0.5) is applied post-decoding in
 * cluster_bins.py, not here, to preserve EM monotonicity.
 */

#include "run_hmm.h"
#include "loglik.h"
#include "fwd_bwd.h"
#include "m_steps.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif


// ---- Internal helpers -------------------------------------------------------

// Transpose (N, M) -> (M, N).
static void transpose_nm_to_mn(const double* src, double* dst, int N, int M)
{
    for (int n = 0; n < N; ++n)
        for (int m = 0; m < M; ++m)
            dst[(long)m * N + n] = src[(long)n * M + m];
}


// Build posts_kn2 (K, N, 2) from posts (N, K, 2).
// Tiled to keep K active rows (K*TILE*2 doubles ≈ 12KB) in L1/L2 cache.
static void build_posts_kn2(const double* posts, double* posts_kn2, int N, int K)
{
    static const int TILE = 32;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int n0 = 0; n0 < N; n0 += TILE) {
        int nEnd = std::min(n0 + TILE, N);
        for (int k = 0; k < K; ++k)
            for (int n = n0; n < nEnd; ++n) {
                posts_kn2[(long)k * N * 2 + (long)n * 2 + 0]
                    = posts[(long)n * K * 2 + (long)k * 2 + 0];
                posts_kn2[(long)k * N * 2 + (long)n * 2 + 1]
                    = posts[(long)n * K * 2 + (long)k * 2 + 1];
            }
    }
}


// Update log_startprobs from posteriors at segment-start bins.
static void update_start_probs(
    const double*  posts,         // (N, K, 2)
    const int64_t* seg_starts,    // (S,)
    int S, int K,
    double tol,
    double* log_startprobs)       // (K, 2) — written in place
{
    std::vector<double> gamma0(K * 2, 0.0);
    for (int s = 0; s < S; ++s) {
        long n = seg_starts[s];
        for (int k = 0; k < K; ++k) {
            gamma0[(long)k * 2 + 0] +=
                std::max(posts[n * K * 2 + (long)k * 2 + 0], tol);
            gamma0[(long)k * 2 + 1] +=
                std::max(posts[n * K * 2 + (long)k * 2 + 1], tol);
        }
    }
    double total = 0.0;
    for (int i = 0; i < K * 2; ++i) total += gamma0[i];
    for (int i = 0; i < K * 2; ++i)
        log_startprobs[i] = std::log(gamma0[i] / total);
}


// ---- Main EM loop -----------------------------------------------------------

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
    bool   share_tau,
    double baf_eps,
    double ig_alpha,
    const double* ig_beta,
    int baf_k_start)
{
    // ---- Precompute seg_starts ----
    std::vector<int64_t> seg_starts(S);
    seg_starts[0] = 0;
    for (int s = 1; s < S; ++s)
        seg_starts[s] = seg_starts[s - 1] + seg_lengths[s - 1];

    // ---- Working arrays (reused each iteration) ----
    std::vector<double> lls0(N * K);
    std::vector<double> lls1(N * K);
    std::vector<double> posts(N * K * 2);
    std::vector<double> posts_nk(N * K);        // marginal over h (N, K)
    std::vector<double> posts_kn2(K * N * 2);   // transposed for BAF M-step

    // ---- Copy initial params (updated in place during EM) ----
    std::vector<double> rdr_means(rdr_means0, rdr_means0 + K * M);
    std::vector<double> rdr_vars(rdr_vars0,   rdr_vars0  + K * M);
    std::vector<double> baf_means(baf_means0, baf_means0 + K * M);
    std::vector<double> baf_taus(baf_taus0,   baf_taus0  + (long)K * M);

    // ---- Uniform log start probs: log(1 / (2K)) ----
    std::vector<double> log_startprobs(K * 2, std::log(1.0 / (2 * K)));

    // ---- Pre-transpose X_alphas/X_betas (N,M) -> (M,N) for BAF M-step ----
    std::vector<double> alphas_mn(M * N);
    std::vector<double> betas_mn(M * N);
    transpose_nm_to_mn(X_alphas, alphas_mn.data(), N, M);
    transpose_nm_to_mn(X_betas,  betas_mn.data(),  N, M);

    // ---- Parameter traces (row 0 = init params) ----
    std::vector<double> trace_rdr_means(rdr_means);
    std::vector<double> trace_rdr_vars(rdr_vars);
    std::vector<double> trace_baf_means(baf_means);
    std::vector<double> trace_baf_taus(baf_taus);

    // ---- EM loop ----
    std::vector<double> elbo_trace;
    elbo_trace.push_back(-std::numeric_limits<double>::infinity());

    double loglik = -std::numeric_limits<double>::infinity();
    int    n_done = n_iter;

    for (int it = 0; it < n_iter; ++it) {

        // E-step: log-likelihood kernel
        compute_loglik_cpp(
            X_rdrs, X_alphas, X_betas, X_totals,
            rdr_means.data(), rdr_vars.data(), baf_means.data(), baf_taus.data(),
            lls0.data(), lls1.data(), N, K, M, share_tau);

        // E-step: forward-backward
        loglik = forward_backward_cpp(
            lls0.data(), lls1.data(),
            seg_starts.data(), seg_lengths,
            log_startprobs.data(),
            log_switchprobs, log_stayprobs,
            log_transmat,
            posts.data(), N, K, S);

        // Compute posts_nk = sum(posts, axis=-1), i.e. (N, K)
        for (int n = 0; n < N; ++n)
            for (int k = 0; k < K; ++k)
                posts_nk[(long)n * K + k] =
                    posts[(long)n * K * 2 + (long)k * 2 + 0]
                  + posts[(long)n * K * 2 + (long)k * 2 + 1];

        // Penalized ELBO (IG log-prior on RDR variance)
        double loglik_penalized = loglik;
        if (ig_alpha > 0 && ig_beta != nullptr) {
            double ig_log_prior = 0.0;
            for (int k = 0; k < K; ++k)
                for (int m = 0; m < M; ++m) {
                    double v = rdr_vars[(long)k * M + m];
                    ig_log_prior += -(ig_alpha + 1.0) * std::log(v) - ig_beta[m] / v;
                }
            loglik_penalized += ig_log_prior;
        }

        // Convergence check (delta normalized to mean per-bin units, so tol_ll
        // is independent of dataset size; matches scikit-learn's mean lower-bound)
        double delta_ll = loglik_penalized - elbo_trace.back();
        elbo_trace.push_back(loglik_penalized);
        if (std::abs(delta_ll) / N < tol_ll) {
            n_done = it + 1;
            break;
        }

        // M-step: start probabilities
        update_start_probs(
            posts.data(), seg_starts.data(), S, K, tol,
            log_startprobs.data());

        // M-step: RDR (closed-form, with IG prior if provided)
        update_rdr_params_cpp(
            X_rdrs, posts_nk.data(),
            rdr_means.data(), rdr_vars.data(),
            N, K, M, min_covar, ig_alpha, ig_beta);

        // M-step: BAF tau (first tau_iters iterations only, before BAF means)
        if (it < tau_iters) {
            update_baf_tau_cpp(
                X_alphas, X_betas, posts.data(), baf_means.data(),
                baf_taus.data(), N, K, M, min_tau, max_tau, share_tau);
        }

        // M-step: BAF means (Brent per (k,m), uses possibly updated tau)
        build_posts_kn2(posts.data(), posts_kn2.data(), N, K);
        update_baf_means_cpp(
            alphas_mn.data(), betas_mn.data(),
            posts_kn2.data(), baf_taus.data(),
            baf_means.data(), N, K, M, baf_eps, baf_k_start);

        trace_rdr_means.insert(trace_rdr_means.end(), rdr_means.begin(), rdr_means.end());
        trace_rdr_vars.insert(trace_rdr_vars.end(), rdr_vars.begin(), rdr_vars.end());
        trace_baf_means.insert(trace_baf_means.end(), baf_means.begin(), baf_means.end());
        trace_baf_taus.insert(trace_baf_taus.end(), baf_taus.begin(), baf_taus.end());
    }

    // Final loglik with fitted params (needed for Viterbi decoding)
    compute_loglik_cpp(
        X_rdrs, X_alphas, X_betas, X_totals,
        rdr_means.data(), rdr_vars.data(), baf_means.data(), baf_taus.data(),
        lls0.data(), lls1.data(), N, K, M, share_tau);

    RunHMMResult result;
    result.rdr_means      = std::move(rdr_means);
    result.rdr_vars       = std::move(rdr_vars);
    result.baf_means      = std::move(baf_means);
    result.baf_taus       = std::move(baf_taus);
    result.log_startprobs = std::move(log_startprobs);
    result.posts          = std::move(posts);
    result.lls0           = std::move(lls0);
    result.lls1           = std::move(lls1);
    result.elbo_trace     = std::move(elbo_trace);
    result.trace_rdr_means  = std::move(trace_rdr_means);
    result.trace_rdr_vars   = std::move(trace_rdr_vars);
    result.trace_baf_means  = std::move(trace_baf_means);
    result.trace_baf_taus   = std::move(trace_baf_taus);
    result.loglik         = result.elbo_trace.back();
    result.data_loglik    = loglik;
    result.n_iters_done   = n_done;
    return result;
}
