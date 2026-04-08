"""Objective builders (Pyomo) and post-solve objective computation (numpy)."""

from __future__ import annotations

import numpy as np
from pyomo import environ as pe

from hatchet.compute_cn.solve.variables import SolverParams


# ---- Pyomo objective builders (used inside build_model) ----


def build_imf_objective(model, p: SolverParams):
    """Weighted L1 deviation: Σ_{m,k} (yA + yB) × w[m]."""
    obj = 0
    for _m in range(p.m):
        cid = p.cluster_ids[_m]
        for _k in range(p.k):
            obj += (model.yA[_m, _k] + model.yB[_m, _k]) * p.w[cid]
    return obj


def build_final_objective(model, obj_imf, obj_reg, penalty_param):
    """Set model.obj combining IMF and regularization."""
    pname = penalty_param[0]
    model.obj_imf = pe.Expression(expr=obj_imf)
    model.obj_reg = pe.Expression(expr=obj_reg)
    if pname == "RAW":
        model.obj = pe.Objective(expr=obj_imf, sense=pe.minimize)
    else:
        model.obj = pe.Objective(
            expr=(1 - model.pparam) * obj_imf + model.pparam * obj_reg,
            sense=pe.minimize,
        )


# ---- Post-solve numpy objective recomputation ----


def compute_obj_IMF(cA, cB, u, f_a, f_b, w):
    """Weighted L1 deviation (numpy)."""
    cA, cB, u = np.array(cA), np.array(cB), np.array(u)
    pred_a = cA @ u
    pred_b = cB @ u
    residual = np.abs(f_a.values - pred_a) + np.abs(f_b.values - pred_b)
    return float(sum(w[cid] * residual[_m].sum() for _m, cid in enumerate(f_a.index)))


def compute_obj_DROOT_SUM(cA, cB, w, base, cluster_ids):
    cA, cB = np.array(cA), np.array(cB)
    obj = 0.0
    for _m, cid in enumerate(cluster_ids):
        for _n in range(1, cA.shape[1]):
            obj += w[cid] * (abs(cA[_m, _n] - base) + abs(cB[_m, _n] - base))
    return obj


def compute_obj_DADJ_SUM(cA, cB, w, base, cluster_ids):
    cA, cB = np.array(cA), np.array(cB)
    n = cA.shape[1]
    obj = 0.0
    for _m, cid in enumerate(cluster_ids):
        for _n1 in range(n - 1):
            for _n2 in range(_n1 + 1, n):
                obj += w[cid] * (
                    abs(cA[_m, _n1] - cA[_m, _n2]) + abs(cB[_m, _n1] - cB[_m, _n2])
                )
    return obj


def compute_obj_MAXCN(cA, cB, w, base, cluster_ids):
    cA, cB = np.array(cA), np.array(cB)
    obj = 0.0
    for _m, cid in enumerate(cluster_ids):
        obj += w[cid] * (max(cA[_m, 1:].max(), base) + max(cB[_m, 1:].max(), base))
    return obj


def compute_obj_DSPAN(cA, cB, w, cluster_ids):
    cA, cB = np.array(cA), np.array(cB)
    obj = 0.0
    for _m, cid in enumerate(cluster_ids):
        tumor = slice(1, None)
        obj += w[cid] * (
            cA[_m, tumor].max()
            - cA[_m, tumor].min()
            + cB[_m, tumor].max()
            - cB[_m, tumor].min()
        )
    return obj


def compute_obj_tree(cA, cB, w, cluster_ids, tree_edges):
    """Total tree edge length (L1 distances along edges)."""
    cA, cB = np.array(cA), np.array(cB)
    obj = 0.0
    for child, parent in tree_edges.items():
        for _m, cid in enumerate(cluster_ids):
            obj += w[cid] * (
                abs(cA[_m, child] - cA[_m, parent])
                + abs(cB[_m, child] - cB[_m, parent])
            )
    return obj
