"""Copy Number Transformation (CNT) distance computation.

Unweighted CND: O(n)-per-chromosome algorithm from Zeira et al. 2017
(PMID 28837352).

Weighted CND: LP-based algorithm from Zeira & Raphael, Bioinformatics 2020
(doi:10.1093/bioinformatics/btaa470). Events have weights based on their
length, position, and type (amplification vs deletion). The semi-ordered
CNT (3 phases: del, amp, del) is formulated as an LP whose constraint
matrix is totally unimodular, guaranteeing integer optimal solutions.
"""

import numpy as np
from scipy.optimize import linprog


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
        Symmetric pairwise CNT distance matrix of shape (n, n).
        Off-diagonal entries are inf when the pair is infeasible.
    """
    m, n = cA.shape
    chr_starts = np.where(chr_boundaries)[0]
    chr_ends = np.append(chr_starts[1:], m)

    dist = np.zeros((n, n), dtype=float)

    for i in range(n):
        for j in range(i + 1, n):
            total = 0.0
            feasible = True
            for cs, ce in zip(chr_starts, chr_ends):
                aiA, ajA = cA[cs:ce, i], cA[cs:ce, j]
                aiB, ajB = cB[cs:ce, i], cB[cs:ce, j]
                # Infeasible if either direction has source=0, target>0
                if (
                    np.any((aiA == 0) & (ajA > 0))
                    or np.any((ajA == 0) & (aiA > 0))
                    or np.any((aiB == 0) & (ajB > 0))
                    or np.any((ajB == 0) & (aiB > 0))
                ):
                    feasible = False
                    break
                total += _cnt_distance_1d(ajA - aiA) + _cnt_distance_1d(ajB - aiB)
            dist[i, j] = total if feasible else np.inf
            dist[j, i] = dist[i, j]

    return dist


# ---------------------------------------------------------------------------
# Weighted CNT distance  (Zeira & Raphael 2020)
# ---------------------------------------------------------------------------


def make_length_weight(seg_lengths=None, beta=1.0, amp_weight=1.0, del_weight=1.0):
    """Create an event weight function based on event length.

    Args:
        seg_lengths: Genomic length of each segment. If None, event length is
            measured in number of segments (k - l + 1).
        beta: Exponent applied to event length.
        amp_weight: Multiplicative factor for amplification events.
        del_weight: Multiplicative factor for deletion events.

    Returns:
        Callable ``w(l, k, tau) -> float`` where *l*, *k* are 0-based segment
        indices and *tau* is +1 (amp) or -1 (del).
    """
    if seg_lengths is not None:
        cum = np.concatenate([[0], np.cumsum(seg_lengths)])

        def w(l, k, tau):
            length = cum[k + 1] - cum[l]
            return (amp_weight if tau == 1 else del_weight) * length**beta

    else:

        def w(l, k, tau):
            length = k - l + 1
            return (amp_weight if tau == 1 else del_weight) * length**beta

    return w


def _weighted_cnt_distance_1d(s, t, w_func):
    """Weighted CNT distance for a single allele on a single chromosome.

    Solves LP 1 from Zeira & Raphael (2020).  The semi-ordered CNT has
    three phases: (1) deletions, (2) amplifications, (3) deletions.

    LP variables: x^j_{lk} — number of events spanning positions l..k in phase j.

    Objective: min  Σ_j Σ_{l≤k} w(l, k, τ_j) · x^j_{lk}

    Constraints (per position i, defining D1/A2/D3 as total phase coverage):
        balance:   s_i − D1_i + A2_i − D3_i = t_i
        (1.1)  t_i = 0  ⇒  D1_i ≥ s_i        (fully delete)
        (1.2)  t_i > 0  ⇒  D1_i ≤ s_i − 1    (don't fully delete)
        (1.4)  phase coverage ≤ B = max(max(S), max(T))
        (1.5)  x ≥ 0

    Args:
        s: Source CN profile, integer array of shape (n,).
        t: Target CN profile, integer array of shape (n,).
        w_func: Weight function ``w(l, k, tau) -> float``.

    Returns:
        Minimum weighted distance, or ``np.inf`` if infeasible (s has a zero
        where t is positive — de novo amplification is undefined in this model).
    """
    n = len(s)
    if n == 0 or np.array_equal(s, t):
        return 0.0
    if np.any((s == 0) & (t > 0)):
        return np.inf

    B = int(max(s.max(), t.max()))

    # Enumerate (l, k) span pairs with l <= k; there are n*(n+1)/2 such pairs
    pairs = [(l, k) for l in range(n) for k in range(l, n)]
    n_pairs = len(pairs)
    n_vars = 3 * n_pairs  # one variable block per phase (del, amp, del)

    # For each position i, collect indices of pairs whose span covers i
    covers = [[] for _ in range(n)]
    for pair_idx, (l, k) in enumerate(pairs):
        for i in range(l, k + 1):
            covers[i].append(pair_idx)

    # Objective: c[var] = weight of that event
    c = np.zeros(n_vars)
    for pair_idx, (l, k) in enumerate(pairs):
        w_del = w_func(l, k, -1)
        w_amp = w_func(l, k, +1)
        c[pair_idx] = w_del  # phase 0 (del)
        c[n_pairs + pair_idx] = w_amp  # phase 1 (amp)
        c[2 * n_pairs + pair_idx] = w_del  # phase 2 (del)

    eq_rows, eq_rhs = [], []
    ub_rows, ub_rhs = [], []

    for i in range(n):
        si, ti = int(s[i]), int(t[i])
        if si == 0 and ti == 0:
            continue

        # Balance constraint: -D1 + A2 - D3 = t_i - s_i
        row = np.zeros(n_vars)
        for pair_idx in covers[i]:
            row[pair_idx] = -1  # phase 0
            row[n_pairs + pair_idx] = 1  # phase 1
            row[2 * n_pairs + pair_idx] = -1  # phase 2
        eq_rows.append(row)
        eq_rhs.append(ti - si)

        if ti == 0:
            # Constraint (1.1): D1_i >= s_i  ⟺  -D1_i <= -s_i
            row = np.zeros(n_vars)
            for pair_idx in covers[i]:
                row[pair_idx] = -1
            ub_rows.append(row)
            ub_rhs.append(-si)
        elif si > 0:
            # Constraint (1.2): D1_i <= s_i - 1
            row = np.zeros(n_vars)
            for pair_idx in covers[i]:
                row[pair_idx] = 1
            ub_rows.append(row)
            ub_rhs.append(si - 1)

    # Constraint (1.4): phase coverage at each position <= B
    for i in range(n):
        for j in range(3):
            row = np.zeros(n_vars)
            for pair_idx in covers[i]:
                row[j * n_pairs + pair_idx] = 1
            ub_rows.append(row)
            ub_rhs.append(B)

    A_eq = np.array(eq_rows) if eq_rows else None
    b_eq = np.array(eq_rhs) if eq_rhs else None
    A_ub = np.array(ub_rows) if ub_rows else None
    b_ub = np.array(ub_rhs) if ub_rhs else None

    res = linprog(
        c,
        A_ub=A_ub,
        b_ub=b_ub,
        A_eq=A_eq,
        b_eq=b_eq,
        bounds=[(0, None)] * n_vars,
        method="highs",
    )
    return res.fun if res.success else np.inf


def compute_weighted_cnt_distances(
    cA,
    cB,
    chr_boundaries,
    seg_lengths=None,
    beta=1.0,
    amp_weight=1.0,
    del_weight=1.0,
):
    """Compute pairwise weighted CNT distances between all clone pairs.

    Because weighted CNT distance is asymmetric, the returned matrix is the
    symmetrised average: ``D(i,j) = (d(i→j) + d(j→i)) / 2``.

    Args:
        cA: Allele-A integer copy numbers, shape (m, n) — m segments, n clones.
        cB: Allele-B integer copy numbers, same shape as cA.
        chr_boundaries: Boolean array of shape (m,); True at the first segment
            of each chromosome.
        seg_lengths: Genomic length of each segment, shape (m,). When None the
            weight function counts segments instead of base-pairs.
        beta: Exponent on event length in the weight function.
        amp_weight: Multiplicative weight for amplification events.
        del_weight: Multiplicative weight for deletion events.

    Returns:
        Symmetric pairwise weighted CNT distance matrix of shape (n, n).
    """
    m, n = cA.shape
    chr_starts = np.where(chr_boundaries)[0]
    chr_ends = np.append(chr_starts[1:], m)

    dist = np.zeros((n, n))

    for i in range(n):
        for j in range(i + 1, n):
            d_ij, d_ji = 0.0, 0.0
            for cs, ce in zip(chr_starts, chr_ends):
                chr_lens = seg_lengths[cs:ce] if seg_lengths is not None else None
                w = make_length_weight(chr_lens, beta, amp_weight, del_weight)

                d_ij += _weighted_cnt_distance_1d(cA[cs:ce, i], cA[cs:ce, j], w)
                d_ij += _weighted_cnt_distance_1d(cB[cs:ce, i], cB[cs:ce, j], w)
                d_ji += _weighted_cnt_distance_1d(cA[cs:ce, j], cA[cs:ce, i], w)
                d_ji += _weighted_cnt_distance_1d(cB[cs:ce, j], cB[cs:ce, i], w)

            dist[i, j] = (d_ij + d_ji) / 2
            dist[j, i] = dist[i, j]

    return dist
