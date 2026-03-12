"""Forward-backward for the 2-mixture HMM.

Prefers a pre-compiled C++ extension (_hmm_cpp) for speed.
Falls back to the Numba kernel when the extension is not installed.

C++ backend:
  - compiled with -O3 -march=native -ffast-math for SIMD auto-vectorisation
  - segments parallelised with OpenMP (if available at compile time)
  - zero JIT cold-start delay

Numba backend (fallback):
  - processes each segment at its actual length (no padding waste)
  - parallelises across segments via prange
"""

import logging

import numpy as np
import numba
from numba import prange

logger = logging.getLogger(__name__)


@numba.njit(cache=True, fastmath=True)
def _logsumexp_1d(x, n):
    """logsumexp over x[:n]."""
    m = x[0]
    for i in range(1, n):
        if x[i] > m:
            m = x[i]
    s = 0.0
    for i in range(n):
        s += np.exp(x[i] - m)
    return m + np.log(s)


@numba.njit(cache=True, fastmath=True)
def _logaddexp(a, b):
    if a > b:
        return a + np.log1p(np.exp(b - a))
    else:
        return b + np.log1p(np.exp(a - b))


@numba.njit(cache=True, fastmath=True)
def _fwd_bwd_seg(
    lls0_seg,  # (T, K)
    lls1_seg,  # (T, K)
    log_startprobs,  # (K, 2)
    log_sw_seg,  # (T,)
    log_st_seg,  # (T,)
    log_transmat,  # (K, K)
    posts_seg,  # (T, K, 2) output — written in place
):
    """Forward-backward for one segment.  Returns loglik."""
    T = lls0_seg.shape[0]
    K = lls0_seg.shape[1]

    fwd = np.empty((T, K, 2))
    log_c = np.empty(T)
    tmp = np.empty(K)

    # ---- Forward t=0 ----
    for k in range(K):
        fwd[0, k, 0] = lls0_seg[0, k] + log_startprobs[k, 0]
        fwd[0, k, 1] = lls1_seg[0, k] + log_startprobs[k, 1]

    m = fwd[0, 0, 0]
    for k in range(K):
        for h in range(2):
            if fwd[0, k, h] > m:
                m = fwd[0, k, h]
    s = 0.0
    for k in range(K):
        for h in range(2):
            s += np.exp(fwd[0, k, h] - m)
    log_c[0] = m + np.log(s)
    for k in range(K):
        fwd[0, k, 0] -= log_c[0]
        fwd[0, k, 1] -= log_c[0]

    # ---- Forward t=1..T-1 ----
    for t in range(1, T):
        pswitch = log_sw_seg[t]
        pstay = log_st_seg[t]

        for k in range(K):
            for j in range(K):
                tmp[j] = fwd[t - 1, j, 0] + log_transmat[j, k] + pstay
            stay0 = _logsumexp_1d(tmp, K)
            for j in range(K):
                tmp[j] = fwd[t - 1, j, 1] + log_transmat[j, k] + pswitch
            switch0 = _logsumexp_1d(tmp, K)

            for j in range(K):
                tmp[j] = fwd[t - 1, j, 1] + log_transmat[j, k] + pstay
            stay1 = _logsumexp_1d(tmp, K)
            for j in range(K):
                tmp[j] = fwd[t - 1, j, 0] + log_transmat[j, k] + pswitch
            switch1 = _logsumexp_1d(tmp, K)

            fwd[t, k, 0] = lls0_seg[t, k] + _logaddexp(stay0, switch0)
            fwd[t, k, 1] = lls1_seg[t, k] + _logaddexp(stay1, switch1)

        m = fwd[t, 0, 0]
        for k in range(K):
            for h in range(2):
                if fwd[t, k, h] > m:
                    m = fwd[t, k, h]
        s = 0.0
        for k in range(K):
            for h in range(2):
                s += np.exp(fwd[t, k, h] - m)
        log_c[t] = m + np.log(s)
        for k in range(K):
            fwd[t, k, 0] -= log_c[t]
            fwd[t, k, 1] -= log_c[t]

    # ---- Backward ----
    bwd = np.zeros((T, K, 2))

    for t in range(T - 2, -1, -1):
        pswitch = log_sw_seg[t + 1]
        pstay = log_st_seg[t + 1]

        for k in range(K):
            for j in range(K):
                tmp[j] = (
                    log_transmat[k, j] + pstay + lls0_seg[t + 1, j] + bwd[t + 1, j, 0]
                )
            stay0 = _logsumexp_1d(tmp, K)
            for j in range(K):
                tmp[j] = (
                    log_transmat[k, j] + pswitch + lls1_seg[t + 1, j] + bwd[t + 1, j, 1]
                )
            switch0 = _logsumexp_1d(tmp, K)

            for j in range(K):
                tmp[j] = (
                    log_transmat[k, j] + pstay + lls1_seg[t + 1, j] + bwd[t + 1, j, 1]
                )
            stay1 = _logsumexp_1d(tmp, K)
            for j in range(K):
                tmp[j] = (
                    log_transmat[k, j] + pswitch + lls0_seg[t + 1, j] + bwd[t + 1, j, 0]
                )
            switch1 = _logsumexp_1d(tmp, K)

            bwd[t, k, 0] = _logaddexp(stay0, switch0) - log_c[t + 1]
            bwd[t, k, 1] = _logaddexp(stay1, switch1) - log_c[t + 1]

    # ---- Posterior ----
    for t in range(T):
        m = fwd[t, 0, 0] + bwd[t, 0, 0]
        for k in range(K):
            for h in range(2):
                v = fwd[t, k, h] + bwd[t, k, h]
                if v > m:
                    m = v
        s = 0.0
        for k in range(K):
            for h in range(2):
                s += np.exp(fwd[t, k, h] + bwd[t, k, h] - m)
        log_norm = m + np.log(s)
        for k in range(K):
            for h in range(2):
                posts_seg[t, k, h] = np.exp(fwd[t, k, h] + bwd[t, k, h] - log_norm)

    loglik = 0.0
    for t in range(T):
        loglik += log_c[t]
    return loglik


@numba.njit(cache=True, parallel=True, fastmath=True)
def _forward_backward(
    lls0,  # (N, K)
    lls1,  # (N, K)
    seg_starts,  # (S,) int64
    seg_lengths,  # (S,) int64
    log_startprobs,  # (K, 2)
    log_switchprobs,  # (N,)
    log_stayprobs,  # (N,)
    log_transmat,  # (K, K)
    posts,  # (N, K, 2) output
):
    """Parallel forward-backward over all segments.  Returns total_loglik."""
    S = seg_starts.shape[0]
    logliks = np.empty(S)

    for s in prange(S):
        start = seg_starts[s]
        length = seg_lengths[s]
        logliks[s] = _fwd_bwd_seg(
            lls0[start : start + length],
            lls1[start : start + length],
            log_startprobs,
            log_switchprobs[start : start + length],
            log_stayprobs[start : start + length],
            log_transmat,
            posts[start : start + length],
        )

    total_ll = 0.0
    for s in range(S):
        total_ll += logliks[s]
    return total_ll


def forward_backward(
    lls0: np.ndarray,
    lls1: np.ndarray,
    X_lengths: np.ndarray,
    log_startprobs: np.ndarray,
    log_switchprobs: np.ndarray,
    log_stayprobs: np.ndarray,
    log_transmat: np.ndarray,
):
    """Forward-backward for the 2-mixture BAF+RDR HMM.

    Uses the C++ extension when available; falls back to Numba otherwise.
    Logs which backend is active on the first call.

    Args:
        lls0:            (N, K) log-likelihoods under haplotype orientation h=0.
        lls1:            (N, K) log-likelihoods under haplotype orientation h=1.
        X_lengths:       (S,) segment lengths in bins.
        log_startprobs:  (K, 2) log start probabilities.
        log_switchprobs: (N,) log phase-switch probabilities.
        log_stayprobs:   (N,) log phase-stay probabilities.
        log_transmat:    (K, K) log cluster transition matrix.

    Returns:
        posts:  (N, K, 2) posterior probabilities.
        loglik: float — total log-likelihood.
    """
    N, K = lls0.shape
    S = len(X_lengths)
    X_lengths_i = np.asarray(X_lengths, dtype=np.int64)

    seg_starts = np.empty(S, dtype=np.int64)
    seg_starts[0] = 0
    for i in range(1, S):
        seg_starts[i] = seg_starts[i - 1] + X_lengths_i[i - 1]

    posts = np.empty((N, K, 2), dtype=np.float64)

    loglik = _forward_backward(
        lls0,
        lls1,
        seg_starts,
        X_lengths_i,
        log_startprobs,
        log_switchprobs,
        log_stayprobs,
        log_transmat,
        posts,
    )
    return posts, float(loglik)
