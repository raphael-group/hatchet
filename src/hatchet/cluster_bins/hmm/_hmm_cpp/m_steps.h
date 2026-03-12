#pragma once

/**
 * C++ BAF M-step: MLE for Beta-Binomial means via Boost.Math Brent's method.
 *
 * Parallelised over (k, m) pairs with OpenMP.
 *
 * alphas_mn  (M, N) C-contiguous — A-allele counts
 * betas_mn   (M, N) C-contiguous — B-allele counts
 * posts_kn2  (K, N, 2) C-contiguous — posteriors: posts_kn2[k*N*2 + n*2 + h]
 * baf_taus   (M,) — per-sample dispersion
 * p_km       (K, M) C-contiguous output — written in place
 * eps        — search bounds [eps, 1-eps]
 */
void update_baf_means_cpp(
    const double* alphas_mn,
    const double* betas_mn,
    const double* posts_kn2,
    const double* baf_taus,
    double*       p_km,
    int N, int K, int M, double eps
);

/**
 * RDR M-step: closed-form weighted Gaussian parameter updates.
 *
 * X_rdrs    (N, M) C-contiguous — RDR observations
 * posts_nk  (N, K) C-contiguous — marginal posteriors (summed over h)
 * rdr_means (K, M) C-contiguous — updated in place
 * rdr_vars  (K, M) C-contiguous — updated in place
 * min_covar — minimum variance floor
 */
void update_rdr_params_cpp(
    const double* X_rdrs,
    const double* posts_nk,
    double*       rdr_means,
    double*       rdr_vars,
    int N, int K, int M, double min_covar,
    double ig_alpha = 10.0,
    const double* ig_beta = nullptr
);

/**
 * BAF tau M-step: MLE of dispersion from cluster-0 bins via Brent in log-space.
 *
 * alphas_nm (N, M) C-contiguous — A-allele counts
 * betas_nm  (N, M) C-contiguous — B-allele counts
 * posts_nk  (N, K) C-contiguous — marginal posteriors
 * baf_taus  (M,)   — updated in place
 * p_fixed   — fixed BAF mean (= 0.5 for cluster 0)
 * min_tau, max_tau — Brent search bounds (optimisation in log-tau space)
 */
void update_baf_tau_cpp(
    const double* alphas_nm,
    const double* betas_nm,
    const double* posts_nk,
    double*       baf_taus,
    int N, int K, int M,
    double p_fixed,
    double min_tau, double max_tau
);
