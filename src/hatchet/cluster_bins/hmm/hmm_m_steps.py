"""EM M-step routines for the 2-mixture BAF+RDR HMM.
- RDR Gaussian means & variances: closed-form weighted statistics.
- BAF Beta-Binomial means: scipy bounded scalar optimization (Brent).
- BAF tau (optional): posterior-weighted MLE over all clusters.
- Start probabilities: posterior counts at segment starts.
"""

import logging

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import betaln


def do_mstep(
    X_rdrs,  # (N, M)
    X_alphas,  # (N, M)
    X_betas,  # (N, M)
    posts,  # (N, K, 2) — posteriors from E-step
    baf_taus,  # (K, M)
    baf_means_init,  # (K, M) — warm start
    X_lengths,  # (S,) segment lengths
    update_tau=False,
    share_tau=True,
    min_covar=1e-3,
    tol=1e-6,
    min_tau=50,
    max_tau=100,
    baf_eps=1e-6,
    ig_alpha=10.0,
    ig_beta=0.01,
    baf_k_start=0,
):
    """EM M-step: emission parameters + start probabilities.

    Updates RDR Gaussian parameters (means and variances) via closed-form
    weighted statistics, BAF Beta-Binomial means via scipy bounded scalar
    optimization, and start probabilities from posterior counts at segment
    starts.  BAF means are NOT folded here; the mhBAF fold is applied
    after decoding in cluster_bins.py to preserve EM monotonicity.

    Returns:
        rdr_means:      (K, M) numpy array.
        rdr_vars:       (K, M) numpy array.
        baf_means:      (K, M) numpy array.
        baf_taus:       (K, M) numpy array.
        log_startprobs: (K, 2) numpy array.
    """
    N, K, _ = posts.shape

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

    # ---- BAF tau (optional, before BAF means so p is optimal for new tau) ----
    if update_tau:
        baf_taus = _update_baf_tau(
            X_alphas,
            X_betas,
            posts,
            baf_means_init,
            baf_taus,
            min_tau=min_tau,
            max_tau=max_tau,
            share_tau=share_tau,
        )

    # ---- BAF means (scipy Brent, uses possibly updated tau) ----
    posts_kn2 = np.ascontiguousarray(posts.transpose(1, 0, 2))  # (K, N, 2)
    baf_means = _update_baf_means(
        baf_means_init,
        X_alphas.T,
        X_betas.T,
        baf_taus,
        posts_kn2,
        baf_eps,
        k_start=baf_k_start,
    )

    return rdr_means, rdr_vars, baf_means, baf_taus, log_startprobs


def _update_baf_tau(
    X_alphas,
    X_betas,
    posts,
    baf_means,
    baf_taus,
    min_tau=50,
    max_tau=500,
    share_tau=True,
):
    """MLE for BAF tau via Brent in log-tau space, maximising Q_BAF.

    With share_tau=True a single tau per sample (pooling all clusters) is fit
    and broadcast to all K rows; with share_tau=False tau is fit independently
    per (cluster, sample).

    Args:
        X_alphas:  (N, M) A-allele counts.
        X_betas:   (N, M) B-allele counts.
        posts:     (N, K, 2) full posteriors.
        baf_means: (K, M) current BAF means.
        baf_taus:  (K, M) current tau values.
        min_tau, max_tau: search bounds.
        share_tau: tie tau across clusters within a sample.

    Returns:
        (K, M) updated tau values.
    """
    N, K, _ = posts.shape
    M = X_alphas.shape[1]
    taus_new = baf_taus.copy()
    lo, hi = np.log(min_tau), np.log(max_tau)

    def neg_Q_km(log_tau, alpha, beta, w0, w1, p):
        tau = np.exp(log_tau)
        a, b = tau * p, tau * (1 - p)
        norm = betaln(a, b)
        ll0 = betaln(alpha + a, beta + b) - norm
        ll1 = betaln(beta + a, alpha + b) - norm
        return -(w0 @ ll0 + w1 @ ll1)

    for m in range(M):
        alpha_m = X_alphas[:, m]
        beta_m = X_betas[:, m]

        if share_tau:

            def neg_Q(
                log_tau, _a=alpha_m, _b=beta_m, _posts=posts, _baf=baf_means[:, m]
            ):
                tau = np.exp(log_tau)
                total = 0.0
                for k in range(K):
                    p = _baf[k]
                    a, b = tau * p, tau * (1 - p)
                    norm = betaln(a, b)
                    ll0 = betaln(_a + a, _b + b) - norm
                    ll1 = betaln(_b + a, _a + b) - norm
                    total += _posts[:, k, 0] @ ll0 + _posts[:, k, 1] @ ll1
                return -total

            res = minimize_scalar(neg_Q, bounds=(lo, hi), method="bounded")
            taus_new[:, m] = np.exp(res.x)
        else:
            for k in range(K):
                res = minimize_scalar(
                    neg_Q_km,
                    bounds=(lo, hi),
                    method="bounded",
                    args=(
                        alpha_m,
                        beta_m,
                        posts[:, k, 0],
                        posts[:, k, 1],
                        baf_means[k, m],
                    ),
                )
                taus_new[k, m] = np.exp(res.x)

    return taus_new


def _update_baf_means(
    p0_km, alphas_mn, betas_mn, baf_taus, posts_kn2, baf_eps=1e-6, k_start=0
):
    """MLE for BAF means via scipy bounded scalar optimization.

    For each (k, m) with k >= k_start, minimizes the posterior-weighted
    negative BB log-likelihood over p in (baf_eps, 1-baf_eps) using Brent.
    Clusters k < k_start retain their initial BAF means.

    Args:
        p0_km:     (K, M) — initial BAF means (preserved for k < k_start).
        alphas_mn: (M, N) — A-allele counts.
        betas_mn:  (M, N) — B-allele counts.
        baf_taus:  (K, M) — dispersion params.
        posts_kn2: (K, N, 2) — posteriors.
        baf_eps:   float  — Brent search bounds [baf_eps, 1-baf_eps].
        k_start:   int    — first cluster index to update (default 0 = all).

    Returns:
        (K, M) BAF means.
    """
    K, M = p0_km.shape
    p_km = p0_km.copy()
    posts0 = posts_kn2[:, :, 0]  # (K, N)
    posts1 = posts_kn2[:, :, 1]  # (K, N)
    EPS = baf_eps

    for m in range(M):
        alpha_m = alphas_mn[m]  # (N,)
        beta_m = betas_mn[m]  # (N,)
        for k in range(k_start, K):
            tau = baf_taus[k, m]
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
