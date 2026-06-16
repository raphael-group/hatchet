#pragma once

/**
 * C++ log-likelihood kernel for the 2-mixture BAF+RDR HMM.
 *
 * Computes per-bin per-cluster emission log-likelihoods under both haplotype
 * orientations (h=0 and h=1).  BAF uses a Beta-Binomial model; RDR uses a
 * Gaussian model.  Parallelised over bins with OpenMP.
 *
 * All arrays are C-contiguous (row-major) float64.
 *
 * X_rdrs    (N, M) — RDR observations
 * X_alphas  (N, M) — A-allele counts
 * X_betas   (N, M) — B-allele counts
 * X_totals  (N, M) — total allele counts
 * rdr_means (K, M) — per-cluster per-sample Gaussian RDR means
 * rdr_vars  (K, M) — per-cluster per-sample Gaussian RDR variances
 * baf_means (K, M) — per-cluster per-sample Beta-Binomial BAF means
 * baf_taus  (K, M) — per-cluster per-sample Beta-Binomial dispersion parameters
 * lls0      (N, K) — output log-likelihoods under h=0
 * lls1      (N, K) — output log-likelihoods under h=1
 */
void compute_loglik_cpp(
    const double* X_rdrs,
    const double* X_alphas,
    const double* X_betas,
    const double* X_totals,
    const double* rdr_means,
    const double* rdr_vars,
    const double* baf_means,
    const double* baf_taus,
    double*       lls0,
    double*       lls1,
    int N, int K, int M, bool share_tau = true
);
