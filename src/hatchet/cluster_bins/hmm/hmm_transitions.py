"""Cluster transition matrix construction for the HMM."""

import numpy as np


def make_transmat(diag, K):
    """Build a K×K doubly-stochastic transition matrix.

    Diagonal entries are `diag`, and all off-diagonal entries share the
    remaining probability mass uniformly: `(1 - diag) / (K - 1)`.

    Args:
        diag: Probability of staying in the same cluster state.
        K:    Number of cluster states.

    Returns:
        (K, K) numpy array.
    """
    if K == 1:
        return np.array([[1.0]], dtype=float)

    offdiag = (1.0 - diag) / (K - 1)
    transmat_ = np.diag([diag - offdiag] * K)
    transmat_ += offdiag
    return transmat_
