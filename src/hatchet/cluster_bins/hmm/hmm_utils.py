"""Model selection scoring and BAF identifiability utilities."""

import numpy as np


def score_BIC(ll: float, K: int, M: int, N: int, share_tau: bool = True):
    """Bayesian Information Criterion for the fitted HMM."""
    num_free_params = 3 * K * M  # RDR means + RDR vars + BAF means
    num_free_params += M if share_tau else K * M  # baf_taus
    num_free_params += 1  # transition diag
    num_free_params += 2 * K - 1  # start probabilities
    return -2.0 * ll + num_free_params * np.log(N)


def score_ICL(
    posts: np.ndarray, ll: float, K: int, M: int, N: int, share_tau: bool = True,
    eps=1e-15,
):
    """ICL = BIC + 2 * classification_entropy; posts: (N, K)"""
    bic = score_BIC(ll, K, M, N, share_tau=share_tau)
    entropy = -np.nansum(posts * np.log(posts + eps))
    return bic + 2.0 * entropy


def convert_mhbafs(baf_means: np.ndarray) -> np.ndarray:
    """
    "mhBAF" = minor-haplotype BAF.  For each cluster, if its mean BAF across
    samples exceeds 0.5 the row is replaced with 1 - BAF, enforcing the
    convention that the reported value is always the minor allele frequency.
    Returns a new array; does not modify baf_means in-place.
    """
    is_minor = np.mean(baf_means, axis=1) <= 0.5  # (K,) bool
    return np.where(is_minor[:, None], baf_means, 1.0 - baf_means)
