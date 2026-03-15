from __future__ import annotations

import logging

import numpy as np
from scipy.special import gammaln
from sklearn.cluster import kmeans_plusplus

from hatchet.cluster_bins.hmm.hmm_likelihoods import compute_loglik_single_cluster_batch
from hatchet.cluster_bins.hmm.hmm_utils import convert_mhbafs
from hatchet.cluster_bins.hmm.hmm_init_utils import *


##################################################
def init_hmm_cna_plus_plus(
    X_rdrs,
    X_bafs,
    X_alphas,
    X_betas,
    X_totals,
    baf_taus0: np.ndarray,
    rdr_vars: np.ndarray,
    K: int,
    random_state=42,
    restarts=10,
    n_local_trials: int | None = None,
    log_rdr=False,
    baf_eps: float = 1e-3,
    collect_diag: bool = False,
):
    """k-means++ style initialization for the HMM emission parameters.

    Seeds K centroids sequentially: first centroid is fixed at (RDR=1, BAF=0.5);
    each subsequent centroid is drawn proportionally to how poorly the current
    centroids explain each bin.  Repeats for `restarts` independent restarts;
    all are returned so the caller can screen with short EM runs and pick the best.

    Returns:
        params_dict: dict mapping restart index → [baf_means, rdr_means, rdr_vars, potential].
    """
    logging.info(
        f"cna++ seeding, K={K}, restarts={restarts}, random_state={random_state}"
    )
    if n_local_trials is None:
        n_local_trials = max(2 + round(np.log(K)), 1)
    N, M = X_rdrs.shape
    unif_probs = np.ones(N) / N
    # Working buffer reused across _weights_from_lls calls to avoid repeated allocation
    nll_weights_buf = np.zeros(N, dtype=np.float64)
    totals_ok = np.all(X_totals > 0, axis=1)

    # Precompute data-only gammaln term once (was recomputed every call)
    log_binom_const = (
        gammaln(X_totals + 1) - gammaln(X_betas + 1) - gammaln(X_alphas + 1)
    )  # (N, M)

    def _weights_from_lls(lls0_full, lls1_full, d=2):
        """Compute sampling weights from per-sample logliks.

        Args:
            lls0_full: (N, k, M) log-likelihoods under haplotype orientation h=0.
            lls1_full: (N, k, M) log-likelihoods under haplotype orientation h=1.
            d: exponent applied to shifted NLL distances (default 2).

        Returns:
            probs_w:   (N,) sampling probability for each bin (sums to 1).
            potential: scalar total potential (sum of D_w); np.inf if degenerate.
            entropy:   scalar Shannon entropy of probs_w; 0.0 if degenerate.
        """
        # Pick the better haplotype orientation per (bin, cluster) pair (sum across samples)
        chosen_h = lls0_full.sum(axis=2) >= lls1_full.sum(axis=2)  # (N, k)
        lls_best = np.where(chosen_h[:, :, None], lls0_full, lls1_full)  # (N, k, M)

        # Chebyshev-style distance: worst sample per cluster, then best cluster per bin
        worst_nll_per_cluster = (-lls_best).max(axis=2)  # (N, k)
        min_nll_per_bin = worst_nll_per_cluster.min(axis=1)  # (N,)

        active = np.isfinite(min_nll_per_bin) & totals_ok
        nll_weights_buf[:] = 0.0
        potential = 0.0
        if np.any(active):
            D_raw = min_nll_per_bin[active]
            D_shifted = D_raw - D_raw.min()
            D_w = D_shifted**d
            potential = float(np.sum(D_w))
            nll_weights_buf[active] = D_w

        if potential <= 0 or not np.isfinite(potential):
            return unif_probs, np.inf, 0.0
        probs_w = nll_weights_buf / nll_weights_buf.sum()
        pos = probs_w > 0
        entropy = float(-np.sum(probs_w[pos] * np.log(probs_w[pos])))
        return probs_w, potential, entropy

    rng = np.random.default_rng(random_state)
    baf_means0 = np.array([[0.5] * M])
    rdr_means0 = np.array([[0.0 if log_rdr else 1.0] * M])
    rdr_vars0 = rdr_vars

    # Compute loglik for the first centroid once — shared across all restarts
    lls0_init, lls1_init = compute_loglik_single_cluster_batch(
        X_rdrs,
        X_alphas,
        X_betas,
        rdr_means0,
        rdr_vars0,
        baf_means0,
        baf_taus0,
        log_binom_const,
    )  # (N, 1, M) each
    probs0, _, _ = _weights_from_lls(lls0_init, lls1_init, d=2)

    params_dict = {}
    diag_dict = {}
    for it in range(restarts):
        logging.debug(f"init_hmm_cna_plus_plus {it}/{restarts}")
        baf_means = baf_means0.copy()
        rdr_means = rdr_means0.copy()
        # Per-cluster variance matrix grown by row-stacking as centroids are chosen
        rdr_vars_k = rdr_vars0.copy()
        probs = probs0
        # Cache per-sample logliks for clusters chosen so far: (N, k, M)
        lls0_cached = lls0_init.copy()
        lls1_cached = lls1_init.copy()
        if collect_diag:
            probs_hist = []
            centroids_hist = []
            selected_bins_hist = []
            candidates_hist = []
        it_potential = np.inf
        for k in range(1, K):
            n_draw = max(min(n_local_trials, int(np.sum(probs > 0))), 1)
            candidates = rng.choice(np.arange(N), size=n_draw, replace=False, p=probs)

            # Batch-compute logliks for all candidate clusters at once
            cand_baf = np.clip(
                X_bafs[candidates, :], baf_eps, 1 - baf_eps
            )  # (n_draw, M)
            cand_rdr = X_rdrs[candidates, :]  # (n_draw, M)
            cand_rdr_vars = np.tile(rdr_vars0[0], (n_draw, 1))  # (n_draw, M)
            cand_lls0, cand_lls1 = compute_loglik_single_cluster_batch(
                X_rdrs,
                X_alphas,
                X_betas,
                cand_rdr,
                cand_rdr_vars,
                cand_baf,
                baf_taus0,
                log_binom_const,
            )  # (N, n_draw, M) each

            # Evaluate potential for each candidate using cached logliks
            best_idx = None
            best_potential = np.inf
            best_j = 0
            for j, cand in enumerate(candidates):
                # Append this candidate's loglik column to the cached set
                full_lls0 = np.concatenate(
                    [lls0_cached, cand_lls0[:, j : j + 1, :]], axis=1
                )  # (N, k+1, M)
                full_lls1 = np.concatenate(
                    [lls1_cached, cand_lls1[:, j : j + 1, :]], axis=1
                )
                _, potential, entropy = _weights_from_lls(full_lls0, full_lls1, d=2)
                logging.debug(
                    "  cand=%d BAF=%s RDR=%s rdr_vars=%s phi=%.4f entropy=%.3f",
                    cand,
                    np.round(X_bafs[cand], 3),
                    np.round(X_rdrs[cand], 3),
                    np.round(cand_rdr_vars[j], 3),
                    potential,
                    entropy,
                )
                if potential < best_potential:
                    best_potential = potential
                    best_idx = cand
                    best_j = j

            if best_idx is None:
                # Fallback: all candidates produced degenerate potentials; take first
                best_idx = int(candidates[0])
                best_j = 0
            baf_means = np.vstack([baf_means, X_bafs[best_idx, :][None, :]])
            rdr_means = np.vstack([rdr_means, X_rdrs[best_idx, :][None, :]])
            rdr_vars_k = np.vstack([rdr_vars_k, rdr_vars0[:1]])

            # Extend cache with the winner's loglik column
            lls0_cached = np.concatenate(
                [lls0_cached, cand_lls0[:, best_j : best_j + 1, :]], axis=1
            )
            lls1_cached = np.concatenate(
                [lls1_cached, cand_lls1[:, best_j : best_j + 1, :]], axis=1
            )

            logging.debug(
                "  k=%d: BAF=%s  RDR=%s",
                k,
                np.round(np.asarray(X_bafs[best_idx, :]), 3),
                np.round(np.asarray(X_rdrs[best_idx, :]), 3),
            )

            if collect_diag:
                probs_hist.append(probs.copy())
                selected_bins_hist.append(best_idx)
                candidates_hist.append((candidates.copy(), best_idx))
                centroids_hist.append((baf_means.copy(), rdr_means.copy()))
            probs, it_potential, _ = _weights_from_lls(lls0_cached, lls1_cached, d=2)

        mhbaf_means = convert_mhbafs(baf_means)
        params_dict[it] = [mhbaf_means, rdr_means, rdr_vars_k, it_potential]
        if collect_diag:
            diag_dict[it] = {
                "probs_history": probs_hist,
                "centroids_history": centroids_hist,
                "selected_bins_history": selected_bins_hist,
                "candidates_history": candidates_hist,
                "final_baf_means": baf_means.copy(),
                "final_rdr_means": rdr_means.copy(),
            }
    return params_dict, diag_dict


##################################################
def init_hmm_kmeans_plus_plus(
    X_rdrs: np.ndarray,
    X_mhbafs: np.ndarray,
    rdr_vars: np.ndarray,
    K: int,
    random_state: int = 42,
    restarts: int = 10,
    n_local_trials: int = 3,
    baf_eps: float = 1e-3,
) -> tuple[dict, dict]:
    """k-means++ seeding in flat [RDR, mhBAF] feature space.

    Runs `restarts` independent k-means++ seedings on the concatenated
    [RDR, mhBAF] feature matrix and returns the resulting centers directly,
    without running k-means EM iterations.

    Args:
        X_rdrs: (N, M) RDR values.
        X_mhbafs: (N, M) mhBAF values, folded to [0, 0.5].
        rdr_vars: (1, M) or (K, M) initial RDR variances — tiled to (K, M) for output.
        K: number of clusters.
        random_state: base random seed; each restart uses random_state + it.
        restarts: number of independent seeding runs.
        n_local_trials: number of candidate points evaluated per seeding step.
        baf_eps: lower bound for BAF means to avoid degenerate Beta-Binomial params.

    Returns:
        params_dict: dict mapping restart index → [baf_means (K,M), rdr_means (K,M), rdr_vars (K,M), potential].
        {}: empty diagnostics dict (matches interface of init_hmm_cna_plus_plus).
    """
    logging.info(
        f"kmeans++ seeding, K={K}, restarts={restarts}, random_state={random_state}"
    )
    N, M = X_rdrs.shape
    X_feat = np.hstack([X_rdrs, X_mhbafs])  # (N, 2M)

    params_dict = {}
    for it in range(restarts):
        rng = np.random.RandomState(random_state + it)
        centers, _ = kmeans_plusplus(
            X_feat, n_clusters=K, random_state=rng, n_local_trials=n_local_trials
        )
        rdr_means = centers[:, :M]  # (K, M)
        baf_means = np.clip(centers[:, M:], baf_eps, 1.0 - baf_eps)  # (K, M)
        rdr_vars_k = (
            np.tile(rdr_vars[0], (K, 1)) if rdr_vars.shape[0] == 1 else rdr_vars.copy()
        )
        dists2 = np.sum(
            (X_feat[:, None, :] - centers[None, :, :]) ** 2, axis=2
        )  # (N, K)
        potential = float(dists2.min(axis=1).sum())
        params_dict[it] = [baf_means, rdr_means, rdr_vars_k, potential]
        logging.debug("kmeans++ restart %d/%d  potential=%.4f", it, restarts, potential)
    return params_dict, {}
