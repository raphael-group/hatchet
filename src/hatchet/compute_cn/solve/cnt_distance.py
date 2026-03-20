"""Copy Number Transformation (CNT) distance computation.

Implements the O(n)-per-chromosome algorithm from Zeira et al. 2017
(PMID 28837352). The CNT distance is the minimum number of contiguous
amplification/deletion events needed to transform one copy-number
profile into another.
"""

import numpy as np


def _cnt_distance_1d(d):
    """Compute CNT distance for a single-allele, single-chromosome difference vector.

    Parameters
    ----------
    d : np.ndarray, shape (n,)
        Difference vector q - p for one allele on one chromosome,
        where p and q are integer CN profiles in genomic order.

    Returns
    -------
    int
        Number of contiguous amplification + deletion events.
    """
    d_pos = np.maximum(d, 0)
    d_neg = np.maximum(-d, 0)

    amp = d_pos[0] + np.sum(np.maximum(np.diff(d_pos), 0))
    dele = d_neg[0] + np.sum(np.maximum(np.diff(d_neg), 0))

    return int(amp + dele)


def compute_cnt_distances(cA, cB, chr_boundaries):
    """Compute pairwise CNT distances between all clone pairs.

    Parameters
    ----------
    cA : np.ndarray, shape (m, n)
        Allele-A integer copy numbers; m segments in genomic order, n clones.
    cB : np.ndarray, shape (m, n)
        Allele-B integer copy numbers; same shape as cA.
    chr_boundaries : np.ndarray, shape (m,), dtype bool
        True at the first segment of each chromosome.

    Returns
    -------
    np.ndarray, shape (n, n)
        Symmetric pairwise CNT distance matrix.
    """
    m, n = cA.shape
    chr_starts = np.where(chr_boundaries)[0]
    chr_ends = np.append(chr_starts[1:], m)

    dist = np.zeros((n, n), dtype=int)

    for i in range(n):
        for j in range(i + 1, n):
            total = 0
            for cs, ce in zip(chr_starts, chr_ends):
                dA = cA[cs:ce, j] - cA[cs:ce, i]
                dB = cB[cs:ce, j] - cB[cs:ce, i]
                total += _cnt_distance_1d(dA) + _cnt_distance_1d(dB)
            dist[i, j] = total
            dist[j, i] = total

    return dist
