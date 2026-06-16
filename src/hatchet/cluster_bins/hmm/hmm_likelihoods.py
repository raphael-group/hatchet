"""Numpy log-likelihood kernels for the 2-mixture BAF+RDR HMM.

Computes per-bin per-cluster emission log-likelihoods under both haplotype
orientations (h=0 and h=1).  BAF uses a Beta-Binomial model; RDR uses a
Gaussian model.  All inputs are C-contiguous float64 numpy arrays.
"""

import numpy as np
from scipy.special import betaln, gammaln


def compute_loglik(
    X_rdrs, X_alphas, X_betas, X_totals, rdr_means, rdr_vars, baf_means, baf_taus
):
    """Per-bin per-cluster log-likelihoods summed over samples.

    Args:
        X_rdrs:    (N, M) observed RDR values.
        X_alphas:  (N, M) A-allele counts.
        X_betas:   (N, M) B-allele counts.
        X_totals:  (N, M) total allele counts.
        rdr_means: (K, M) per-cluster per-sample Gaussian RDR means.
        rdr_vars:  (K, M) per-cluster per-sample Gaussian RDR variances.
        baf_means: (K, M) per-cluster per-sample Beta-Binomial BAF means.
        baf_taus:  (K, M) or (M,) Beta-Binomial dispersion parameters.

    Returns:
        lls0: (N, K) log-likelihoods under haplotype orientation h=0.
        lls1: (N, K) log-likelihoods under haplotype orientation h=1.
    """
    log_binom_const = (
        gammaln(X_totals + 1) - gammaln(X_betas + 1) - gammaln(X_alphas + 1)
    )  # (N, M)
    bt = baf_taus[None, :, :] if baf_taus.ndim == 2 else baf_taus[None, None, :]
    bb_alpha = bt * baf_means[None, :, :]  # (1, K, M)
    bb_beta = bt * (1 - baf_means[None, :, :])
    bb_delta = betaln(bb_alpha, bb_beta)
    lnB_h0 = betaln(X_alphas[:, None, :] + bb_alpha, X_betas[:, None, :] + bb_beta)
    lnB_h1 = betaln(X_betas[:, None, :] + bb_alpha, X_alphas[:, None, :] + bb_beta)
    ll_bafs_h0 = np.sum(
        log_binom_const[:, None, :] + lnB_h0 - bb_delta, axis=2
    )  # (N, K)
    ll_bafs_h1 = np.sum(log_binom_const[:, None, :] + lnB_h1 - bb_delta, axis=2)

    log_norm_const = 0.5 * np.sum(np.log(2 * np.pi * rdr_vars), axis=1)  # (K,)
    quad = 0.5 * np.einsum(
        "nkm,km->nk",
        (X_rdrs[:, None, :] - rdr_means[None, :, :]) ** 2,
        1.0 / rdr_vars,
    )  # (N, K)
    ll_rdrs = -quad - log_norm_const

    return ll_rdrs + ll_bafs_h0, ll_rdrs + ll_bafs_h1


def compute_loglik_single_cluster_batch(
    X_rdrs,
    X_alphas,
    X_betas,
    rdr_means_batch,
    rdr_vars_batch,
    baf_means_batch,
    baf_taus,
    log_binom_const,
):
    """Per-bin per-sample log-likelihoods for a batch of candidate clusters.

    Same emission model as compute_loglik but for C candidate clusters at
    once, reusing a precomputed log_binom_const to avoid redundant gammaln
    evaluations across calls.

    Args:
        X_rdrs:           (N, M) observed RDR values.
        X_alphas:         (N, M) A-allele counts.
        X_betas:          (N, M) B-allele counts.
        rdr_means_batch:  (C, M) RDR means for each candidate cluster.
        rdr_vars_batch:   (C, M) RDR variances for each candidate cluster.
        baf_means_batch:  (C, M) BAF means for each candidate cluster.
        baf_taus:         (M,)   per-sample Beta-Binomial dispersion.
        log_binom_const:  (N, M) precomputed gammaln(T+1)-gammaln(B+1)-gammaln(A+1).

    Returns:
        lls0: (N, C, M) log-likelihoods under haplotype orientation h=0.
        lls1: (N, C, M) log-likelihoods under haplotype orientation h=1.
    """
    bb_alpha = baf_taus[None, None, :] * baf_means_batch[None, :, :]  # (1, C, M)
    bb_beta = baf_taus[None, None, :] * (1 - baf_means_batch[None, :, :])  # (1, C, M)
    bb_delta = betaln(bb_alpha, bb_beta)  # (1, C, M)
    lnB_h0 = betaln(X_alphas[:, None, :] + bb_alpha, X_betas[:, None, :] + bb_beta)
    lnB_h1 = betaln(X_betas[:, None, :] + bb_alpha, X_alphas[:, None, :] + bb_beta)
    ll_baf_h0 = log_binom_const[:, None, :] + lnB_h0 - bb_delta  # (N, C, M)
    ll_baf_h1 = log_binom_const[:, None, :] + lnB_h1 - bb_delta  # (N, C, M)

    diff2 = (X_rdrs[:, None, :] - rdr_means_batch[None, :, :]) ** 2  # (N, C, M)
    ll_rdr = (
        -0.5 * diff2 / rdr_vars_batch[None, :, :]
        - 0.5 * np.log(2 * np.pi * rdr_vars_batch)[None, :, :]
    )  # (N, C, M)

    return ll_rdr + ll_baf_h0, ll_rdr + ll_baf_h1
