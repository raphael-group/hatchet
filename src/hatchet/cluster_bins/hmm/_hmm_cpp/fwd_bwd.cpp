/**
 * Forward-backward algorithm for the 2-mixture HMM.
 *
 * Segments parallelised with OpenMP. Log-space normalisation via logsumexp.
 */

#include "fwd_bwd.h"

#include <cmath>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif


static inline double logaddexp(double a, double b) {
    if (a > b) return a + std::log1p(std::exp(b - a));
    return       b + std::log1p(std::exp(a - b));
}

static inline double logsumexp_k(const double* x, int K) {
    double m = x[0];
    for (int i = 1; i < K; ++i) if (x[i] > m) m = x[i];
    double s = 0.0;
    for (int i = 0; i < K; ++i) s += std::exp(x[i] - m);
    return m + std::log(s);
}


/**
 * Forward-backward for one segment.
 *
 * lls0_seg  (T, K),  lls1_seg  (T, K)
 * log_startprobs (K, 2)
 * log_sw_seg (T,),   log_st_seg (T,)
 * log_transmat (K, K)
 * posts_seg  (T, K, 2) — written in place
 *
 * Returns segment log-likelihood.
 */
static double fwd_bwd_seg_cpp(
    const double* lls0_seg,
    const double* lls1_seg,
    const double* log_startprobs,
    const double* log_sw_seg,
    const double* log_st_seg,
    const double* log_transmat,
    double*       posts_seg,
    int T, int K)
{
    std::vector<double> fwd(T * K * 2);
    std::vector<double> bwd(T * K * 2, 0.0);
    std::vector<double> log_c(T);
    std::vector<double> tmp(K);

    // ---- Forward t=0 ----
    for (int k = 0; k < K; ++k) {
        fwd[k*2 + 0] = lls0_seg[k] + log_startprobs[k*2 + 0];
        fwd[k*2 + 1] = lls1_seg[k] + log_startprobs[k*2 + 1];
    }
    {
        double m = fwd[0];
        for (int k = 0; k < K; ++k)
            for (int h = 0; h < 2; ++h)
                if (fwd[k*2 + h] > m) m = fwd[k*2 + h];
        double s = 0.0;
        for (int k = 0; k < K; ++k)
            for (int h = 0; h < 2; ++h)
                s += std::exp(fwd[k*2 + h] - m);
        log_c[0] = m + std::log(s);
        for (int k = 0; k < K; ++k) {
            fwd[k*2 + 0] -= log_c[0];
            fwd[k*2 + 1] -= log_c[0];
        }
    }

    // ---- Forward t=1..T-1 ----
    for (int t = 1; t < T; ++t) {
        double pswitch = log_sw_seg[t];
        double pstay   = log_st_seg[t];
        const double* prev = &fwd[(t-1) * K * 2];
        double*       curr = &fwd[t     * K * 2];

        for (int k = 0; k < K; ++k) {
            // h=0: stay from h=0, or switch from h=1
            for (int j = 0; j < K; ++j)
                tmp[j] = prev[j*2 + 0] + log_transmat[j*K + k] + pstay;
            double stay0 = logsumexp_k(tmp.data(), K);
            for (int j = 0; j < K; ++j)
                tmp[j] = prev[j*2 + 1] + log_transmat[j*K + k] + pswitch;
            double switch0 = logsumexp_k(tmp.data(), K);

            // h=1: stay from h=1, or switch from h=0
            for (int j = 0; j < K; ++j)
                tmp[j] = prev[j*2 + 1] + log_transmat[j*K + k] + pstay;
            double stay1 = logsumexp_k(tmp.data(), K);
            for (int j = 0; j < K; ++j)
                tmp[j] = prev[j*2 + 0] + log_transmat[j*K + k] + pswitch;
            double switch1 = logsumexp_k(tmp.data(), K);

            curr[k*2 + 0] = lls0_seg[t*K + k] + logaddexp(stay0, switch0);
            curr[k*2 + 1] = lls1_seg[t*K + k] + logaddexp(stay1, switch1);
        }

        // Normalize
        double m = curr[0];
        for (int k = 0; k < K; ++k)
            for (int h = 0; h < 2; ++h)
                if (curr[k*2 + h] > m) m = curr[k*2 + h];
        double s = 0.0;
        for (int k = 0; k < K; ++k)
            for (int h = 0; h < 2; ++h)
                s += std::exp(curr[k*2 + h] - m);
        log_c[t] = m + std::log(s);
        for (int k = 0; k < K; ++k) {
            curr[k*2 + 0] -= log_c[t];
            curr[k*2 + 1] -= log_c[t];
        }
    }

    // ---- Backward (initialised to 0) ----
    for (int t = T-2; t >= 0; --t) {
        double pswitch = log_sw_seg[t+1];
        double pstay   = log_st_seg[t+1];
        const double* next_bwd  = &bwd[(t+1) * K * 2];
        double*       curr_bwd  = &bwd[t     * K * 2];
        const double* next_lls0 = &lls0_seg[(t+1) * K];
        const double* next_lls1 = &lls1_seg[(t+1) * K];

        for (int k = 0; k < K; ++k) {
            // bwd[t, k, 0]: stay h=0 → h=0 (pstay), or switch h=0 → h=1 (pswitch)
            for (int j = 0; j < K; ++j)
                tmp[j] = log_transmat[k*K + j] + pstay + next_lls0[j] + next_bwd[j*2 + 0];
            double stay0 = logsumexp_k(tmp.data(), K);
            for (int j = 0; j < K; ++j)
                tmp[j] = log_transmat[k*K + j] + pswitch + next_lls1[j] + next_bwd[j*2 + 1];
            double switch0 = logsumexp_k(tmp.data(), K);

            // bwd[t, k, 1]: stay h=1 → h=1 (pstay), or switch h=1 → h=0 (pswitch)
            for (int j = 0; j < K; ++j)
                tmp[j] = log_transmat[k*K + j] + pstay + next_lls1[j] + next_bwd[j*2 + 1];
            double stay1 = logsumexp_k(tmp.data(), K);
            for (int j = 0; j < K; ++j)
                tmp[j] = log_transmat[k*K + j] + pswitch + next_lls0[j] + next_bwd[j*2 + 0];
            double switch1 = logsumexp_k(tmp.data(), K);

            curr_bwd[k*2 + 0] = logaddexp(stay0, switch0) - log_c[t+1];
            curr_bwd[k*2 + 1] = logaddexp(stay1, switch1) - log_c[t+1];
        }
    }

    // ---- Posterior ----
    for (int t = 0; t < T; ++t) {
        const double* f = &fwd[t * K * 2];
        const double* b = &bwd[t * K * 2];
        double*       p = &posts_seg[t * K * 2];

        double m = f[0] + b[0];
        for (int k = 0; k < K; ++k)
            for (int h = 0; h < 2; ++h) {
                double v = f[k*2 + h] + b[k*2 + h];
                if (v > m) m = v;
            }
        double s = 0.0;
        for (int k = 0; k < K; ++k)
            for (int h = 0; h < 2; ++h)
                s += std::exp(f[k*2 + h] + b[k*2 + h] - m);
        double log_norm = m + std::log(s);
        for (int k = 0; k < K; ++k)
            for (int h = 0; h < 2; ++h)
                p[k*2 + h] = std::exp(f[k*2 + h] + b[k*2 + h] - log_norm);
    }

    double loglik = 0.0;
    for (int t = 0; t < T; ++t) loglik += log_c[t];
    return loglik;
}


double forward_backward_cpp(
    const double*  lls0,
    const double*  lls1,
    const int64_t* seg_starts,
    const int64_t* seg_lengths,
    const double*  log_startprobs,
    const double*  log_switchprobs,
    const double*  log_stayprobs,
    const double*  log_transmat,
    double*        posts,
    int N, int K, int S)
{
    double total_ll = 0.0;

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) reduction(+:total_ll)
#endif
    for (int s = 0; s < S; ++s) {
        int64_t start = seg_starts[s];
        int64_t len   = seg_lengths[s];
        double seg_ll = fwd_bwd_seg_cpp(
            lls0 + start * K,
            lls1 + start * K,
            log_startprobs,
            log_switchprobs + start,
            log_stayprobs   + start,
            log_transmat,
            posts + start * K * 2,
            (int)len, K
        );
        total_ll += seg_ll;
    }

    return total_ll;
}
