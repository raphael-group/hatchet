"""Copy Number Transformation (CNT) distance computation.

Unweighted CND: O(n)-per-chromosome algorithm from Zeira et al. 2017
(PMID 28837352).
"""

import numpy as np


def _cnt_distance_1d(d):
    """Compute CNT distance for a single-allele, single-chromosome difference vector.

    Args:
        d: Difference vector q - p for one allele on one chromosome,
            where p and q are integer CN profiles in genomic order. Shape (n,).

    Returns:
        Number of contiguous amplification + deletion events.
    """
    d_pos = np.maximum(d, 0)
    d_neg = np.maximum(-d, 0)

    amp = d_pos[0] + np.sum(np.maximum(np.diff(d_pos), 0))
    dele = d_neg[0] + np.sum(np.maximum(np.diff(d_neg), 0))

    return int(amp + dele)


def compute_cnt_distances(cA, cB, chr_boundaries):
    """Compute pairwise CNT distances between all clone pairs.

    A pair is infeasible (distance = inf) when any chromosome has a segment
    where one clone has CN=0 and the other has CN>0, since the CNT model
    requires a non-zero source CN to amplify from.

    Args:
        cA: Allele-A integer copy numbers, shape (m, n) — m segments, n clones.
        cB: Allele-B integer copy numbers, same shape as cA.
        chr_boundaries: Boolean array of shape (m,); True at the first segment
            of each chromosome.

    Returns:
        Asymmetric pairwise CNT distance matrix of shape (n, n).
        ``dist[i, j]`` is the distance from clone i to clone j.
        Returns ``np.inf`` for infeasible (source allele=0, target>0).
        Off-diagonal entries are inf when the pair is infeasible.
    """
    m, n = cA.shape
    chr_starts = np.where(chr_boundaries)[0]
    chr_ends = np.append(chr_starts[1:], m)

    dist = np.zeros((n, n), dtype=float)

    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            total = 0.0
            feasible = True
            for cs, ce in zip(chr_starts, chr_ends):
                sA, tA = cA[cs:ce, i], cA[cs:ce, j]
                sB, tB = cB[cs:ce, i], cB[cs:ce, j]
                # Infeasible: source allele is 0 but target > 0
                if np.any((sA == 0) & (tA > 0)) or np.any((sB == 0) & (tB > 0)):
                    feasible = False
                    break
                total += _cnt_distance_1d(tA - sA) + _cnt_distance_1d(tB - sB)
            dist[i, j] = total if feasible else np.inf

    return dist
