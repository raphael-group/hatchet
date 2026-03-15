"""EM M-step routines for the 2-mixture BAF+RDR HMM.
- RDR Gaussian means & variances: closed-form weighted statistics.
- BAF Beta-Binomial means: scipy bounded scalar optimization (Brent).
- BAF tau (optional): MLE from cluster-0 posterior.
- Start probabilities: posterior counts at segment starts.
- Transition parameter (optional): from xi sufficient statistics.
"""

import logging

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import betaln

from hatchet.cluster_bins.cluster_utils import mle_BB_dispersion
from hatchet.cluster_bins.hmm.hmm_utils import convert_mhbafs


def do_mstep(
    X_rdrs,  # (N, M)
    X_alphas,  # (N, M)
    X_betas,  # (N, M)
    posts,  # (N, K, 2) — posteriors from E-step
    baf_taus,  # (M,)
    baf_means_init,  # (K, M) — warm start
    X_lengths,  # (S,) segment lengths
    update_tau=False,
    min_covar=1e-3,
    tol=1e-6,
    min_tau=50,
    max_tau=100,
    baf_eps=1e-6,
    ig_alpha=10.0,
    ig_beta=0.01,
):
    """EM M-step: emission parameters + start probabilities.

    Updates RDR Gaussian parameters (means and variances) via closed-form
    weighted statistics, BAF Beta-Binomial means via scipy bounded scalar
    optimization, and start probabilities from posterior counts at segment
    starts.  BAF means are NOT folded here; the mhBAF identifiability
    constraint is applied once after EM converges (in run_hmm) so that the
    ELBO is strictly monotone throughout training.

    Returns:
        rdr_means:      (K, M) numpy array.
        rdr_vars:       (K, M) numpy array.
        baf_means:      (K, M) numpy array.
        baf_taus:       (M,)   numpy array.
        log_startprobs: (K, 2) numpy array.
    """
    N, K, _ = posts.shape
    M = X_rdrs.shape[1]

    # ---- start probabilities ----
    seg_starts = np.concatenate([[0], np.cumsum(X_lengths[:-1])])
    gamma0 = np.sum(np.maximum(posts[seg_starts], tol), axis=0)  # (K, 2)
    log_startprobs = np.log(gamma0 / np.sum(gamma0))  # (K, 2)

    # ---- RDR means and variances (closed-form) ----
    posts_marg = np.sum(posts, axis=-1)  # (N, K)
    Nk = np.sum(posts_marg, axis=0)  # (K,)

    if np.any(Nk < tol):
        logging.warning(f"Some clusters have effective Nk < {tol}")

    safe_Nk = np.maximum(Nk, tol)[:, None]  # (K, 1)
    rdr_means = np.einsum("nk,nm->km", posts_marg, X_rdrs) / safe_Nk
    weighted_var = (
        np.einsum("nk,nm->km", posts_marg, X_rdrs**2) / safe_Nk - rdr_means**2
    )
    if ig_alpha > 0:
        raw_SS = weighted_var * safe_Nk  # (K, M)
        rdr_vars = np.maximum(
            (raw_SS + 2 * ig_beta) / (Nk[:, None] + 2 * (ig_alpha + 1)),
            min_covar,
        )
    else:
        rdr_vars = np.maximum(weighted_var, min_covar)

    # ---- BAF means (scipy Brent) ----
    posts_kn2 = np.ascontiguousarray(posts.transpose(1, 0, 2))  # (K, N, 2)
    baf_means = _update_baf_means(
        baf_means_init, X_alphas.T, X_betas.T, baf_taus, posts_kn2, baf_eps
    )
    baf_means = convert_mhbafs(baf_means)

    # ---- BAF tau (optional) — estimate from k=0 (BAF=0.5) only ----
    if update_tau:
        mask_0 = np.argmax(posts_marg, axis=1) == 0
        if np.sum(mask_0) > 1:
            taus_new = baf_taus.copy()
            for m in range(M):
                taus_new[m] = mle_BB_dispersion(
                    X_alphas[mask_0, m],
                    X_betas[mask_0, m],
                    p=0.5,
                    min_tau=min_tau,
                    max_tau=max_tau,
                )
            baf_taus = taus_new

    return rdr_means, rdr_vars, baf_means, baf_taus, log_startprobs


def _update_baf_means(p0_km, alphas_mn, betas_mn, baf_taus, posts_kn2, baf_eps=1e-6):
    """MLE for BAF means via scipy bounded scalar optimization.

    For each (k, m), minimizes the posterior-weighted negative BB log-likelihood
    over p in (baf_eps, 1-baf_eps) using Brent's method.

    Args:
        p0_km:     (K, M) — warm start (unused by Brent, kept for output shape).
        alphas_mn: (M, N) — A-allele counts.
        betas_mn:  (M, N) — B-allele counts.
        baf_taus:  (M,)   — dispersion params.
        posts_kn2: (K, N, 2) — posteriors.
        baf_eps:   float  — Brent search bounds [baf_eps, 1-baf_eps].

    Returns:
        (K, M) BAF means.
    """
    K, M = p0_km.shape
    p_km = np.empty_like(p0_km)
    posts0 = posts_kn2[:, :, 0]  # (K, N)
    posts1 = posts_kn2[:, :, 1]  # (K, N)
    EPS = baf_eps

    for m in range(M):
        tau = baf_taus[m]
        alpha_m = alphas_mn[m]  # (N,)
        beta_m = betas_mn[m]  # (N,)
        for k in range(K):
            w0 = posts0[k]  # (N,)
            w1 = posts1[k]  # (N,)

            def neg_Q(p, _tau=tau, _a=alpha_m, _b=beta_m, _w0=w0, _w1=w1):
                a = _tau * p
                b = _tau * (1.0 - p)
                # h=0: BB(alpha, beta | p, tau)
                ll0 = betaln(_a + a, _b + b) - betaln(a, b)
                # h=1: BB(beta, alpha | p, tau)
                ll1 = betaln(_b + a, _a + b) - betaln(a, b)
                return -(_w0 @ ll0 + _w1 @ ll1)

            res = minimize_scalar(neg_Q, bounds=(EPS, 1.0 - EPS), method="bounded")
            p_km[k, m] = res.x

    return p_km
