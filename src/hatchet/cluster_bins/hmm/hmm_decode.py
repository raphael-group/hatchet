"""Decoding routines for the 2-mixture HMM.

Provides MAP decoding (from marginal posteriors) and Viterbi decoding
(numba-accelerated, jointly over K cluster states x 2 phase states).
"""

import numpy as np
import numba


def map_decoding(phase_posts: np.ndarray, cluster_posts: np.ndarray):
    """Maximum-a-posteriori (MAP) decoding from marginal posteriors.

    Args:
        phase_posts:   (N, 2) marginal posteriors over haplotype orientation.
        cluster_posts: (N, K) marginal posteriors over cluster states.

    Returns:
        cluster_labels: (N,) int32 — MAP cluster assignment per bin.
        phase_labels:   (N,) int8  — MAP haplotype orientation (0 or 1) per bin.
    """
    phase_labels = np.argmax(phase_posts, axis=1)
    cluster_labels = np.argmax(cluster_posts, axis=1)
    return cluster_labels, phase_labels


@numba.njit(cache=True)
def _viterbi_segment(
    lls0_seg,
    lls1_seg,
    log_startprobs,
    log_switchprobs_seg,
    log_stayprobs_seg,
    log_transmat,
):
    """Numba-accelerated Viterbi for one segment.

    Args:
        lls0_seg, lls1_seg: (T, K) log-emission likelihoods for h=0 and h=1.
        log_startprobs:     (K, 2) log start probabilities.
        log_switchprobs_seg: (T,) log phase-switch probabilities.
        log_stayprobs_seg:  (T,) log phase-stay probabilities.
        log_transmat:       (K, K) log cluster transition matrix.

    Returns:
        path_k: (T,) int64 — Viterbi cluster path.
        path_h: (T,) int64 — Viterbi haplotype-orientation path.
    """
    nobs = lls0_seg.shape[0]
    K = lls0_seg.shape[1]

    delta = np.empty((nobs, K, 2))
    psi_k = np.empty((nobs, K, 2), dtype=np.int64)
    psi_h = np.empty((nobs, K, 2), dtype=np.int64)

    for k in range(K):
        delta[0, k, 0] = lls0_seg[0, k] + log_startprobs[k, 0]
        delta[0, k, 1] = lls1_seg[0, k] + log_startprobs[k, 1]

    for t in range(1, nobs):
        pswitch = log_switchprobs_seg[t]
        pstay = log_stayprobs_seg[t]

        for k in range(K):
            # h=0: arrive from (prev_k, h=0) via stay or (prev_k, h=1) via switch
            best_val = -np.inf
            best_k = 0
            best_h = 0
            for prev_k in range(K):
                v = delta[t - 1, prev_k, 0] + log_transmat[prev_k, k] + pstay
                if v > best_val:
                    best_val = v
                    best_k = prev_k
                    best_h = 0
                v = delta[t - 1, prev_k, 1] + log_transmat[prev_k, k] + pswitch
                if v > best_val:
                    best_val = v
                    best_k = prev_k
                    best_h = 1
            delta[t, k, 0] = lls0_seg[t, k] + best_val
            psi_k[t, k, 0] = best_k
            psi_h[t, k, 0] = best_h

            # h=1: arrive from (prev_k, h=1) via stay or (prev_k, h=0) via switch
            best_val = -np.inf
            best_k = 0
            best_h = 0
            for prev_k in range(K):
                v = delta[t - 1, prev_k, 1] + log_transmat[prev_k, k] + pstay
                if v > best_val:
                    best_val = v
                    best_k = prev_k
                    best_h = 1
                v = delta[t - 1, prev_k, 0] + log_transmat[prev_k, k] + pswitch
                if v > best_val:
                    best_val = v
                    best_k = prev_k
                    best_h = 0
            delta[t, k, 1] = lls1_seg[t, k] + best_val
            psi_k[t, k, 1] = best_k
            psi_h[t, k, 1] = best_h

    # Backtrace
    best_val = -np.inf
    best_t_k = 0
    best_t_h = 0
    for k in range(K):
        for h in range(2):
            if delta[nobs - 1, k, h] > best_val:
                best_val = delta[nobs - 1, k, h]
                best_t_k = k
                best_t_h = h

    path_k = np.empty(nobs, dtype=np.int64)
    path_h = np.empty(nobs, dtype=np.int64)
    path_k[nobs - 1] = best_t_k
    path_h[nobs - 1] = best_t_h
    for t in range(nobs - 2, -1, -1):
        pk = psi_k[t + 1, path_k[t + 1], path_h[t + 1]]
        ph = psi_h[t + 1, path_k[t + 1], path_h[t + 1]]
        path_k[t] = pk
        path_h[t] = ph

    return path_k, path_h


def run_viterbi(
    lls0: np.ndarray,
    lls1: np.ndarray,
    X_lengths: np.ndarray,
    log_startprobs: np.ndarray,
    log_switchprobs: np.ndarray,
    log_stayprobs: np.ndarray,
    log_transmat: np.ndarray,
    K: int,
    N: int,
):
    """Viterbi decoding over all segments (numba-accelerated).

    Runs the classic Viterbi algorithm independently for each chromosomal
    segment, jointly over K cluster states and 2 haplotype-orientation states.
    The joint state space is (K x 2); transitions factor as:
        P(k', h' | k, h) = transmat[k, k'] x (stay if h'==h else switch).

    Args:
        lls0, lls1:      (N, K) log-emission likelihoods for h=0 and h=1.
        X_lengths:       (S,) number of bins per segment.
        log_startprobs:  (K, 2) log start probabilities (tied across segments).
        log_switchprobs: (N,) log probability of phase switch at each bin.
        log_stayprobs:   (N,) log probability of staying in same phase.
        log_transmat:    (K, K) log cluster transition matrix.
        K:               Number of cluster states.
        N:               Total number of bins.

    Returns:
        cluster_labels: (N,) int32 — Viterbi cluster path.
        phase_labels:   (N,) int8  — Viterbi haplotype-orientation path.
    """
    # Ensure contiguous float64 arrays for numba
    lls0 = np.ascontiguousarray(lls0, dtype=np.float64)
    lls1 = np.ascontiguousarray(lls1, dtype=np.float64)
    log_startprobs = np.ascontiguousarray(log_startprobs, dtype=np.float64)
    log_switchprobs = np.ascontiguousarray(log_switchprobs, dtype=np.float64)
    log_stayprobs = np.ascontiguousarray(log_stayprobs, dtype=np.float64)
    log_transmat = np.ascontiguousarray(log_transmat, dtype=np.float64)

    cluster_labels = np.empty(N, dtype=np.int32)
    phase_labels = np.empty(N, dtype=np.int8)
    start = 0
    for s, nobs in enumerate(X_lengths):
        end = start + nobs
        path_k, path_h = _viterbi_segment(
            lls0[start:end],
            lls1[start:end],
            log_startprobs,
            log_switchprobs[start:end],
            log_stayprobs[start:end],
            log_transmat,
        )
        cluster_labels[start:end] = path_k
        phase_labels[start:end] = path_h
        start = end

    return cluster_labels, phase_labels


def decode_hmm(
    sol, decode_method, X_lengths, log_switchprobs, log_stayprobs, log_transmat
):
    """Decode cluster and phase labels from a fitted HMM solution.

    Args:
        sol:             Dict returned by run_hmm.
        decode_method:   "map" or "viterbi".
        X_lengths:       (S,) segment lengths in bins.
        log_switchprobs: (N,) log phase-switch probabilities.
        log_stayprobs:   (N,) log phase-stay probabilities.
        log_transmat:    (K, K) log cluster transition matrix.

    Returns:
        cluster_labels: (N,) int32
        phase_labels:   (N,) int8
    """
    assert decode_method in ["map", "viterbi"]
    if decode_method == "map":
        return map_decoding(sol["phase_posts"], sol["cluster_posts"])
    K = sol["cluster_posts"].shape[1]
    N = sol["cluster_posts"].shape[0]
    return run_viterbi(
        sol["lls0"],
        sol["lls1"],
        X_lengths,
        sol["log_startprobs"],
        log_switchprobs,
        log_stayprobs,
        log_transmat,
        K,
        N,
    )
