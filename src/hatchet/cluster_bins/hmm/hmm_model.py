"""Baum-Welch EM training for the phased 2-mixture Beta-Binomial BAF + Gaussian RDR HMM."""

from __future__ import annotations

import time
import logging
import numpy as np

from hatchet.cluster_bins.hmm.hmm_fwd_bwd import forward_backward
from hatchet.cluster_bins.hmm.hmm_likelihoods import compute_loglik
from hatchet.cluster_bins.hmm.hmm_decode import map_decoding, run_viterbi, decode_hmm  # noqa: F401
from hatchet.cluster_bins.hmm.hmm_utils import score_BIC, score_ICL  # noqa: F401
from hatchet.cluster_bins.hmm.hmm_m_steps import do_mstep
from hatchet.cluster_bins.cluster_utils import count_multimodal_clusters
from hatchet.cluster_bins.hmm import _USE_CPP, _cpp_run_hmm


def _check_elbo_convergence(elbo_trace, tol_ll, prefix=""):
    """Check elbo_trace for convergence and monotonicity; log warnings."""
    trace = elbo_trace[1:]  # skip initial -inf
    if len(trace) < 2:
        return
    final_delta = abs(trace[-1] - trace[-2])
    if final_delta >= tol_ll:
        logging.warning(
            f"{prefix}EM did NOT converge: final delta_ll={final_delta:.6e} >= tol_ll={tol_ll:.6e}"
        )
    drops = []
    for i in range(1, len(trace)):
        delta = trace[i] - trace[i - 1]
        if delta < -tol_ll:
            drops.append((i, delta))
    if drops:
        logging.warning(
            f"{prefix}ELBO decreased at {len(drops)} iteration(s): "
            + ", ".join(f"iter {i} (delta={d:.6e})" for i, d in drops[:5])
        )


##################################################
def _run_hmm_cpp(
    K,
    X_rdrs,
    X_alphas,
    X_betas,
    X_totals,
    X_lengths,
    log_switchprobs,
    log_stayprobs,
    log_transmat,
    rdr_means,
    rdr_vars,
    baf_means,
    baf_taus,
    n_iter,
    min_covar,
    tol_ll,
    tol,
    tau_iters,
    min_tau,
    max_tau,
    baf_eps,
    _pfx,
    ig_alpha=10.0,
    ig_beta=0.01,
) -> dict:
    """Thin Python wrapper around the C++ full EM loop."""
    t0 = time.perf_counter()

    ig_beta_arr = np.atleast_1d(np.asarray(ig_beta, dtype=np.float64))
    res = _cpp_run_hmm(
        K,
        X_rdrs,
        X_alphas,
        X_betas,
        X_totals,
        X_lengths,
        log_switchprobs,
        log_stayprobs,
        log_transmat,
        rdr_means,
        rdr_vars,
        baf_means,
        baf_taus,
        n_iter,
        min_covar,
        tol_ll,
        tol,
        tau_iters,
        min_tau,
        max_tau,
        baf_eps,
        ig_alpha,
        ig_beta_arr,
    )

    elapsed = time.perf_counter() - t0
    n_done = res["n_iters_done"]
    logging.info(
        f"{_pfx}HMM baum-welch (K={K}) finished after {n_done} iters | "
        f"loglik={res['model_ll']:.6f} | total={elapsed:.2f}s ({elapsed / n_done:.3f}s/it)"
    )
    _check_elbo_convergence(res["elbo_trace"], tol_ll, _pfx)

    full_posts = res["full_posts"]  # (N, K, 2)
    phase_posts = np.sum(full_posts, axis=1)  # (N, 2)
    cluster_posts = np.sum(full_posts, axis=2)  # (N, K)

    return {
        "RDR_means": res["RDR_means"],
        "RDR_vars": res["RDR_vars"],
        "BAF_means": res["BAF_means"],
        "BAF_taus": res["BAF_taus"],
        "elbo_trace": res["elbo_trace"],
        "obj_ll": res["obj_ll"],
        "model_ll": res["model_ll"],
        "log_startprobs": res["log_startprobs"],
        "log_transmat": log_transmat,
        "full_posts": full_posts,
        "phase_posts": phase_posts,
        "cluster_posts": cluster_posts,
        "lls0": res["lls0"],
        "lls1": res["lls1"],
    }


##################################################
def run_hmm(
    K: int,
    X_rdrs: np.ndarray,
    X_alphas: np.ndarray,
    X_betas: np.ndarray,
    X_totals: np.ndarray,
    X_lengths: np.ndarray,
    log_switchprobs: np.ndarray,
    log_stayprobs: np.ndarray,
    log_transmat: np.ndarray,
    rdr_means: np.ndarray,
    rdr_vars: np.ndarray,
    baf_means: np.ndarray,
    baf_taus: np.ndarray,
    X_rdrs_orig: np.ndarray,
    X_totals_orig: np.ndarray,
    n_iter: int = 10,
    min_covar: float = 1e-3,
    tol_ll: float = 1e-4,
    tol: float = 1e-6,
    tau_iters: int = 1,
    min_tau: float = 50,
    max_tau: float = 100,
    baf_eps: float = 1e-3,
    log_rdr: bool = True,
    restart_id: int | None = None,
    ig_alpha: float = 10.0,
    ig_beta: float | np.ndarray = 0.01,
) -> dict:
    """Run EM training for a K-state 2-mixture BAF+RDR HMM.

    All data arrays are C-contiguous float64 numpy (int64 for X_lengths).
    Alternates between:
      - E-step: numpy loglik kernel + Numba forward-backward.
      - M-step: closed-form RDR Gaussian updates + scipy Brent BAF updates.

    The first ``tau_iters`` EM iterations also update the Beta-Binomial
    dispersion parameters ``baf_taus``.

    Args:
        K:                 Number of cluster states.
        X_rdrs:         (N, M) float64 — (log-)RDR observations.
        X_alphas:       (N, M) float64 — A-allele counts.
        X_betas:        (N, M) float64 — B-allele counts.
        X_totals:       (N, M) float64 — total allele counts.
        X_lengths:         (S,)  int64    — segment lengths in bins.
        log_switchprobs: (N,) float64  — log phase-switch probabilities.
        log_stayprobs:  (N,)  float64  — log phase-stay probabilities.
        log_transmat:      (K, K) float64 — log cluster transition matrix.
        rdr_means:         (K, M) float64 — initial RDR means.
        rdr_vars:          (K, M) float64 — initial RDR variances.
        baf_means:         (K, M) float64 — initial BAF means.
        baf_taus:          (M,)   float64 — initial BB dispersion.
        X_rdrs_orig:       (N, M) float64 — linear-scale RDR for multimodal diagnostic.
        X_totals_orig:     (N, M) float64 — total counts for multimodal diagnostic.
        n_iter:            Maximum EM iterations.
        min_covar:         Minimum RDR variance floor.
        tol_ll:            Convergence threshold on log-likelihood delta.
        tol:               Minimum effective cluster size.
        tau_iters:         Number of iterations during which tau is updated.
        min_tau:           Lower bound for tau optimisation.
        max_tau:           Upper bound for tau optimisation.
        baf_eps:           Brent search bounds for BAF mean: [baf_eps, 1-baf_eps].
        log_rdr:           Whether RDR model params are in log-space.

    Returns:
        dict with keys: RDR_means, RDR_vars, BAF_means, BAF_taus,
        elbo_trace, model_ll, lls0 (N,K), lls1 (N,K), log_startprobs (K,2),
        full_posts (N,K,2), phase_posts (N,2), cluster_posts (N,K).
    """
    assert n_iter > 1

    _pfx = f"[r{restart_id}] " if restart_id is not None else ""
    N, M = X_rdrs.shape

    if _USE_CPP:
        return _run_hmm_cpp(
            K=K,
            X_rdrs=X_rdrs,
            X_alphas=X_alphas,
            X_betas=X_betas,
            X_totals=X_totals,
            X_lengths=X_lengths,
            log_switchprobs=log_switchprobs,
            log_stayprobs=log_stayprobs,
            log_transmat=log_transmat,
            rdr_means=rdr_means,
            rdr_vars=rdr_vars,
            baf_means=baf_means,
            baf_taus=baf_taus,
            n_iter=n_iter,
            min_covar=min_covar,
            tol_ll=tol_ll,
            tol=tol,
            tau_iters=tau_iters,
            min_tau=min_tau,
            max_tau=max_tau,
            baf_eps=baf_eps,
            _pfx=_pfx,
            ig_alpha=ig_alpha,
            ig_beta=ig_beta,
        )

    log_startprobs = np.log(np.full((K, 2), 1.0 / (2 * K), dtype=np.float64))

    elbo_trace = [-np.inf]
    t_loglik_sum = 0.0
    t_fwdbwd_sum = 0.0
    t_mstep_sum = 0.0
    for it in range(n_iter):
        t0 = time.perf_counter()

        lls0, lls1 = compute_loglik(
            X_rdrs,
            X_alphas,
            X_betas,
            X_totals,
            rdr_means,
            rdr_vars,
            baf_means,
            baf_taus,
        )
        t1_loglik = time.perf_counter()

        posts, loglik = forward_backward(
            lls0,
            lls1,
            X_lengths,
            log_startprobs,
            log_switchprobs,
            log_stayprobs,
            log_transmat,
        )
        t2_fwdbwd = time.perf_counter()

        rdr_means, rdr_vars, baf_means, baf_taus, log_startprobs = do_mstep(
            X_rdrs,
            X_alphas,
            X_betas,
            posts,
            baf_taus,
            baf_means,
            X_lengths,
            update_tau=(it < tau_iters),
            min_covar=min_covar,
            tol=tol,
            min_tau=min_tau,
            max_tau=max_tau,
            baf_eps=baf_eps,
            ig_alpha=ig_alpha,
            ig_beta=ig_beta,
        )
        t3_mstep = time.perf_counter()

        t_loglik_sum += t1_loglik - t0
        t_fwdbwd_sum += t2_fwdbwd - t1_loglik
        t_mstep_sum += t3_mstep - t2_fwdbwd

        if ig_alpha > 0:
            ig_log_prior = np.sum(
                -(ig_alpha + 1) * np.log(rdr_vars) - ig_beta / rdr_vars
            )
            loglik_penalized = loglik + ig_log_prior
        else:
            loglik_penalized = loglik

        delta_ll = loglik_penalized - elbo_trace[-1]
        elbo_trace.append(loglik_penalized)

        # Per-iter multimodal diagnostic (cheap MAP decode from posteriors)
        cluster_posts_it = np.sum(posts, axis=2)  # (N, K)
        phase_posts_it = np.sum(posts, axis=1)  # (N, 2)
        labels_it = np.argmax(cluster_posts_it, axis=1)
        phases_it = np.argmax(phase_posts_it, axis=1)
        betas_phased = (
            X_alphas * (1 - phases_it[:, None]) + X_betas * phases_it[:, None]
        )
        bafs_it = betas_phased / X_totals_orig
        n_multi, multi_ids = count_multimodal_clusters(
            labels_it, X_rdrs_orig, bafs_it, log_rdr
        )
        n_used = len(np.unique(labels_it))
        multi_info = f" | {n_multi}/{n_used} multimodal {multi_ids}"

        logging.info(
            f"{_pfx}Iter {it:03d} | Q={loglik: .6f} | delta={delta_ll: .6f}{multi_info}"
        )

        if logging.getLogger().isEnabledFor(logging.DEBUG):
            tau_str = " ".join(f"{baf_taus[m]:8.1f}" for m in range(M))
            tau_hdr = "  ".join(f"{'tau' + str(m + 1):>8s}" for m in range(M))
            logging.debug(f"  {'tau':>3s}        {tau_hdr}")
            logging.debug(f"       {' ' * 5}  {tau_str}")
            sample_hdr = "  ".join(
                f"{'baf' + str(m + 1):>8s} {'rdr' + str(m + 1):>8s} {'rdr_var' + str(m + 1):>8s}"
                for m in range(M)
            )
            nks = cluster_posts_it.sum(axis=0)
            total_nk = nks.sum()
            logging.debug(f"  {'k':>3s}  {'pi_k':>6s}  {sample_hdr}")
            for k in range(K):
                nk = float(nks[k])
                if nk < 0.5:
                    continue
                per_sample = "  ".join(
                    f"{baf_means[k, m]:8.4f} {rdr_means[k, m]:8.4f} {rdr_vars[k, m]:8.4f}"
                    for m in range(M)
                )
                logging.debug(f"  {k:3d}  {nk / total_nk:6.3f}  {per_sample}")

        if abs(delta_ll) < tol_ll:
            logging.info(f"{_pfx}Converged at iteration {it}")
            break

    n_done = it + 1
    logging.info(
        f"{_pfx}EM profile ({n_done} iters): "
        f"loglik={t_loglik_sum:.2f}s ({t_loglik_sum / n_done:.3f}s/it) | "
        f"fwdbwd={t_fwdbwd_sum:.2f}s ({t_fwdbwd_sum / n_done:.3f}s/it) | "
        f"mstep={t_mstep_sum:.2f}s ({t_mstep_sum / n_done:.3f}s/it)"
    )
    _check_elbo_convergence(elbo_trace, tol_ll, _pfx)

    obj_ll = elbo_trace[-1]
    model_ll = loglik
    phase_posts = np.sum(posts, axis=1)  # (N, 2)
    cluster_posts = np.sum(posts, axis=2)  # (N, K)

    lls0_final, lls1_final = compute_loglik(
        X_rdrs,
        X_alphas,
        X_betas,
        X_totals,
        rdr_means,
        rdr_vars,
        baf_means,
        baf_taus,
    )

    return {
        "RDR_means": rdr_means,
        "RDR_vars": rdr_vars,
        "BAF_means": baf_means,
        "BAF_taus": baf_taus,
        "elbo_trace": elbo_trace,
        "obj_ll": obj_ll,
        "model_ll": model_ll,
        "log_startprobs": log_startprobs,
        "log_transmat": log_transmat,
        "full_posts": posts,
        "phase_posts": phase_posts,
        "cluster_posts": cluster_posts,
        "lls0": lls0_final,
        "lls1": lls1_final,
    }
