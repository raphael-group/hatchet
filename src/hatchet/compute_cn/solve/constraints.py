"""Constraint builders for the copy-number ILP model.

Each function takes a Pyomo model (with variables already created) and a
SolverParams, and adds constraints to model.constraints (a ConstraintList).
"""

from __future__ import annotations

import math

from pyomo import environ as pe

from hatchet.compute_cn.solve.variables import SolverParams


def add_l1_constraints(model, p: SolverParams):
    """L1 linearisation: yA >= |f_a_obs - fA| for all (m, k)."""
    f_a_vals = p.f_a.values
    f_b_vals = p.f_b.values
    for _m in range(p.m):
        for _k in range(p.k):
            fa_obs = float(f_a_vals[_m, _k])
            fb_obs = float(f_b_vals[_m, _k])
            model.constraints.add(fa_obs - model.fA[_m, _k] <= model.yA[_m, _k])
            model.constraints.add(model.fA[_m, _k] - fa_obs <= model.yA[_m, _k])
            model.constraints.add(fb_obs - model.fB[_m, _k] <= model.yB[_m, _k])
            model.constraints.add(model.fB[_m, _k] - fb_obs <= model.yB[_m, _k])


def add_mixture_constraints(
    model, p: SolverParams, fixed_u=None, fixed_cA=None, fixed_cB=None
):
    """Define fA/fB = Σ cA/cB × u, mode-dependent.

    FULL:   bilinear via McCormick (bit-encoded cA × u)
    CARCH:  linear in cA/cB (u fixed via fixed_u)
    UARCH:  linear in u (cA/cB fixed via fixed_cA/fixed_cB)
    """
    m, n, k = p.m, p.n, p.k
    _M = p.M

    if p.mode == "FULL":
        for _m in range(m):
            for _k in range(k):
                for fX, cX, bitcX, vX in [
                    (model.fA, model.cA, model.bitcA, model.vA),
                    (model.fB, model.cB, model.bitcB, model.vB),
                ]:
                    if _m in p.fixed_rows:
                        s = sum(cX[_m, _n] * model.u[_n, _k] for _n in range(n))
                    else:
                        s = cX[_m, 0] * model.u[0, _k]
                        for _n in p.tumor_clones:
                            for _b in range(_M):
                                s += vX[_b, _m, _n, _k] * math.pow(2, _b)
                                model.constraints.add(
                                    vX[_b, _m, _n, _k] <= bitcX[_b, _m, _n]
                                )
                                model.constraints.add(
                                    vX[_b, _m, _n, _k] <= model.u[_n, _k]
                                )
                                model.constraints.add(
                                    vX[_b, _m, _n, _k]
                                    >= bitcX[_b, _m, _n] + model.u[_n, _k] - 1
                                )
                    model.constraints.add(fX[_m, _k] == s)

        # Non-zero CN: if clone has proportion, must have CN somewhere
        for _n in p.tumor_clones:
            for _k in range(k):
                _sum = sum(model.cA[_m, _n] + model.cB[_m, _n] for _m in range(m))
                model.constraints.add(_sum >= model.u[_n, _k])

    elif p.mode == "CARCH":
        for _m in range(m):
            for _k in range(k):
                sA = sum(
                    model.cA[_m, _n] * fixed_u[_n][_k]
                    for _n in range(n)
                    if fixed_u[_n][_k] >= p.minprop - p.tol
                )
                sB = sum(
                    model.cB[_m, _n] * fixed_u[_n][_k]
                    for _n in range(n)
                    if fixed_u[_n][_k] >= p.minprop - p.tol
                )
                model.constraints.add(model.fA[_m, _k] == sA)
                model.constraints.add(model.fB[_m, _k] == sB)

    elif p.mode == "UARCH":
        for _m in range(m):
            for _k in range(k):
                sA = sum(int(fixed_cA[_m][_n]) * model.u[_n, _k] for _n in range(n))
                sB = sum(int(fixed_cB[_m][_n]) * model.u[_n, _k] for _n in range(n))
                model.constraints.add(model.fA[_m, _k] == sA)
                model.constraints.add(model.fB[_m, _k] == sB)


def add_bit_encoding(model, p: SolverParams):
    """Bit decomposition: cA[m,n] = Σ_b bitcA[b,m,n] × 2^b."""
    if not ((p.mode == "FULL") or (p.max_ncns_seg > 0 and p.mode == "CARCH")):
        return
    _M = p.M
    for _m in p.free_rows:
        for _n in p.tumor_clones:
            for cX, bitcX in [(model.cA, model.bitcA), (model.cB, model.bitcB)]:
                s = sum(bitcX[_b, _m, _n] * math.pow(2, _b) for _b in range(_M))
                model.constraints.add(cX[_m, _n] == s)


def add_proportion_constraints(model, p: SolverParams):
    """Clone proportions sum to 1; minprop enforcement."""
    if p.mode not in ("FULL", "UARCH"):
        return
    for _k in range(p.k):
        model.constraints.add(sum(model.u[_n, _k] for _n in range(p.n)) == 1)
    if p.minprop > 0:
        for _n in p.tumor_clones:
            for _k in range(p.k):
                model.constraints.add(model.x[_n, _k] >= model.u[_n, _k])
                model.constraints.add(model.u[_n, _k] >= p.minprop * model.x[_n, _k])


def add_domain_constraints(model, p: SolverParams):
    """Upper bound, zero CN, amp/del, balanced, MRCA LOH, symmetry breaking."""
    if p.mode not in ("FULL", "CARCH"):
        return
    # CN upper bound
    for _m in p.free_rows:
        cid = p.cluster_ids[_m]
        bound = max(sum(p.copy_numbers.get(cid, (0, 0))), p.cn_max)
        for _n in p.tumor_clones:
            model.constraints.add(model.cA[_m, _n] + model.cB[_m, _n] <= bound)

    # Zero CN prevention
    w_total = sum(p.w)
    for _m in p.free_rows:
        cid = p.cluster_ids[_m]
        if p.w[cid] / w_total >= p.zero_cn_thres:
            for _n in p.tumor_clones:
                model.constraints.add(model.cA[_m, _n] + model.cB[_m, _n] >= 1)

    # Amp/del
    if p.ampdel:
        for _m in p.free_rows:
            for _n in p.tumor_clones:
                model.constraints.add(
                    model.cA[_m, _n]
                    <= p.cn_max * model.adA[_m] + p.base - p.base * model.adA[_m]
                )
                model.constraints.add(model.cA[_m, _n] >= p.base * model.adA[_m])
                model.constraints.add(
                    model.cB[_m, _n]
                    <= p.cn_max * model.adB[_m] + p.base - p.base * model.adB[_m]
                )
                model.constraints.add(model.cB[_m, _n] >= p.base * model.adB[_m])

    # Balanced: cA == cB
    if p.balanced_clusters:
        bal_set = set(p.balanced_clusters)
        for _m in p.free_rows:
            if p.cluster_ids[_m] in bal_set:
                for _n in p.tumor_clones:
                    model.constraints.add(model.cA[_m, _n] == model.cB[_m, _n])

    # MRCA LOH
    if p.mrca and p.n >= 3:
        for _m in p.free_rows:
            for _n in range(2, p.n):
                model.constraints.add(model.cA[_m, _n] <= model.cA[_m, 1] * p.cn_max)
                model.constraints.add(model.cB[_m, _n] <= model.cB[_m, 1] * p.cn_max)

    # Symmetry breaking
    for i in range(1, p.n - 1):
        sum1 = sum(
            (model.cA[_m, i] + model.cB[_m, i]) * p.symmCoeff(_m) for _m in range(p.m)
        )
        sum2 = sum(
            (model.cA[_m, i + 1] + model.cB[_m, i + 1]) * p.symmCoeff(_m)
            for _m in range(p.m)
        )
        model.constraints.add(sum1 <= sum2)


def add_ncns_seg_constraints(model, p: SolverParams):
    """Max-NCNS-Seg grouping constraints."""
    if not ((p.mode in ("FULL", "CARCH")) and p.max_ncns_seg > 0):
        return
    _M = p.M
    for _m in p.free_rows:
        for _n in p.tumor_clones:
            model.constraints.add(
                sum(model.z[_m, _n, _d] for _d in range(p.max_ncns_seg)) == 1
            )

    for _m in p.free_rows:
        for _b in range(_M):
            for _d in range(p.max_ncns_seg):
                for _i in range(1, p.n - 1):
                    for _j in range(1, p.n):
                        for bitcX in [model.bitcA, model.bitcB]:
                            model.constraints.add(
                                bitcX[_b, _m, _i] - bitcX[_b, _m, _j]
                                <= 2 - model.z[_m, _i, _d] - model.z[_m, _j, _d]
                            )
                            model.constraints.add(
                                bitcX[_b, _m, _j] - bitcX[_b, _m, _i]
                                <= 2 - model.z[_m, _i, _d] - model.z[_m, _j, _d]
                            )

    for _m in p.free_rows:
        for _d in range(p.max_ncns_seg - 1):
            _sum_l = _sum_l1 = 0
            for _n in p.tumor_clones:
                _sum_l += model.z[_m, _n, _d] * p.symmCoeff(_n)
                _sum_l1 += model.z[_m, _n, _d + 1] * p.symmCoeff(_n)
                model.constraints.add(_sum_l <= _sum_l1)
