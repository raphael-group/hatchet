import logging
import numpy as np


def split_by_chromosome(inputs) -> list[np.ndarray]:
    """Split segment indices by chromosome using chr_boundaries."""
    S = inputs.m
    if inputs.chr_boundaries is None:
        return [np.arange(S)]
    chrom_groups = []
    start = 0
    for i in range(1, S):
        if inputs.chr_boundaries[i]:
            chrom_groups.append(np.arange(start, i))
            start = i
    chrom_groups.append(np.arange(start, S))
    return chrom_groups


def dedup_pool_instances(pool_instances, u_atol=1e-3):
    """Remove duplicate solutions that differ only by clone ordering.

    Two solutions are duplicates if their CN profiles (cA, cB) are identical
    after sorting clones into a canonical order and their clone proportions (u)
    agree within ``u_atol``.

    Args:
        pool_instances: dict of {sol_id: {"imf_obj": ..., "cA": ..., "cB": ..., "u": ...}}.
        u_atol: Absolute tolerance for comparing clone proportions.

    Returns:
        Deduplicated dict (order preserved, first occurrence kept).
    """
    if len(pool_instances) <= 1:
        return pool_instances

    def _canon(cA, cB, u):
        arr = np.array(cA + cB)
        u_arr = np.array(u)
        order = np.lexsort(arr[::-1])
        return arr[:, order], u_arr[order]

    keep = {}
    canon_cache = []
    for sol_id, sol in pool_instances.items():
        cn_c, u_c = _canon(sol["cA"], sol["cB"], sol["u"])
        is_dup = False
        for cn_k, u_k in canon_cache:
            if np.array_equal(cn_c, cn_k) and np.allclose(u_c, u_k, atol=u_atol):
                is_dup = True
                break
        if not is_dup:
            keep[sol_id] = sol
            canon_cache.append((cn_c, u_c))
    n_removed = len(pool_instances) - len(keep)
    if n_removed > 0:
        logging.info(
            f"removed {n_removed} duplicate solutions (same CN up to clone reordering)"
        )
    return keep


def compute_individual_objs(pname, weights, fA, fB, cA, cB, u):
    """Compute the IMF and regularisation objective values for one solution.

    Returns [imf_obj, reg_obj].
    """
    w_ = weights.to_numpy().reshape((len(weights), 1))
    fA_ = fA.to_numpy()
    fB_ = fB.to_numpy()
    cA_ = np.array(cA)
    cB_ = np.array(cB)
    u_ = np.array(u)

    imf_obj = compute_obj_IMF(w_, fA_, fB_, cA_, cB_, u_)
    reg_objs = {
        "MAXCN": compute_obj_MAXCN,
        "DBOX_L1": compute_obj_DBOX_L1,
        "DBOX_L0": compute_obj_DBOX_L0,
        "DROOT_SUM": compute_obj_DROOT_SUM,
        "DADJ_SUM": compute_obj_DADJ_SUM,
    }
    sub_obj = reg_objs[pname](w_, fA_, fB_, cA_, cB_, u_) if pname in reg_objs else 0.0
    return [imf_obj, sub_obj]


def compute_obj_IMF(weights, fA, fB, cA, cB, u):
    """Compute the weighted IMF (integer matrix factorisation) objective."""
    leftA_w = weights * np.abs(fA - cA @ u)
    leftB_w = weights * np.abs(fB - cB @ u)
    obj = np.sum(leftA_w) + np.sum(leftB_w)
    return obj


def compute_obj_DROOT_SUM(weights, _fA, _fB, cA, cB, _u):
    """Compute weighted sum of CN distance from the diploid root (1,1) for tumor clones.

    All regularisation objective functions share the same call signature
    ``(weights, fA, fB, cA, cB, u)`` so they can be dispatched uniformly via a
    dict in ``compute_individual_objs``.  Parameters not needed by this objective
    are prefixed with ``_``.
    """
    distA = weights * np.abs(cA[:, 1:] - cA[:, :1])
    distB = weights * np.abs(cB[:, 1:] - cB[:, :1])
    obj = np.sum(distA) + np.sum(distB)
    return obj


def compute_obj_DADJ_SUM(weights, _fA, _fB, cA, cB, _u):
    """Compute weighted pairwise Hamming distance between clone CN states.

    All regularisation objective functions share the same call signature
    ``(weights, fA, fB, cA, cB, u)`` so they can be dispatched uniformly via a
    dict in ``compute_individual_objs``.  Parameters not needed by this objective
    are prefixed with ``_``.
    """
    obj = 0
    num_clusters, num_clones = cA.shape
    for cluster_idx in range(num_clusters):
        obj_cluster = 0.0
        for clone_i in range(num_clones - 1):
            for clone_j in range(clone_i + 1, num_clones):
                obj_cluster += abs(cA[cluster_idx, clone_i] - cA[cluster_idx, clone_j])
                obj_cluster += abs(cB[cluster_idx, clone_i] - cB[cluster_idx, clone_j])
        obj += weights[cluster_idx, 0] * obj_cluster
    return obj


def compute_obj_MAXCN(weights, _fA, _fB, cA, cB, _u):
    """Compute weighted sum of maximum CN per cluster across tumor clones.

    All regularisation objective functions share the same call signature
    ``(weights, fA, fB, cA, cB, u)`` so they can be dispatched uniformly via a
    dict in ``compute_individual_objs``.  Parameters not needed by this objective
    are prefixed with ``_``.
    """
    maxA_w = np.dot(np.max(cA[:, 1:], axis=1), weights)[0]
    maxB_w = np.dot(np.max(cB[:, 1:], axis=1), weights)[0]
    return maxA_w + maxB_w


def compute_obj_DBOX_L1(weights, _fA, _fB, cA, cB, _u):
    """Compute weighted sum of allelic span (max-min) per cluster across tumor clones."""
    spanA = np.max(cA[:, 1:], axis=1) - np.min(cA[:, 1:], axis=1)
    spanB = np.max(cB[:, 1:], axis=1) - np.min(cB[:, 1:], axis=1)
    return float(np.dot(spanA + spanB, weights)[0])


def compute_obj_DBOX_L0(weights, _fA, _fB, cA, cB, _u):
    """L0((a_max-a_min)+(b_max-b_min)) per cluster: 1 if any allelic span > 0."""
    spanA = np.max(cA[:, 1:], axis=1) - np.min(cA[:, 1:], axis=1)
    spanB = np.max(cB[:, 1:], axis=1) - np.min(cB[:, 1:], axis=1)
    is_sub = ((spanA + spanB) > 0).astype(float)
    return float(np.dot(is_sub, weights)[0])
