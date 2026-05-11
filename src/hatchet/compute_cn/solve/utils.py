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
