"""Regularization term builders for the copy-number ILP model.

Each function adds auxiliary variables/constraints to *model* and returns a
Pyomo expression for the regularization objective component.
"""

from __future__ import annotations

import numpy as np
from pyomo import environ as pe

from hatchet.compute_cn.solve.datatypes import SolverParams, SolverInputs


def build_reg_maxcn(model, mode: str, params: SolverParams, inputs: SolverInputs):
    """MAXCN: penalise max CN per cluster."""
    obj = 0
    hcn_idx = [(m_, ab) for m_ in inputs.free_rows for ab in ("a", "b")]
    model.hcn = pe.Var(hcn_idx, bounds=(0, np.inf), domain=pe.Reals)
    for _m in inputs.free_rows:
        for _n in params.tumor_clones:
            model.constraints.add(model.cA[_m, _n] <= model.hcn[_m, "a"])
            model.constraints.add(model.cB[_m, _n] <= model.hcn[_m, "b"])
        cid = inputs.cluster_ids[_m]
        obj += inputs.w[cid] * (model.hcn[_m, "a"] + model.hcn[_m, "b"])
    for _m in inputs.fixed_rows:
        cid = inputs.cluster_ids[_m]
        ca, cb = inputs.copy_numbers[cid]
        obj += inputs.w[cid] * (max(ca, params.base) + max(cb, params.base))
    return obj


def build_reg_droot_sum(model, mode: str, params: SolverParams, inputs: SolverInputs):
    """DROOT_SUM: L1 distance from normal clone (n=0)."""
    n = params.n
    obj = 0
    md_idx = [
        (_m, _n, ab)
        for _m in inputs.free_rows
        for _n in params.tumor_clones
        for ab in ("a", "b")
    ]
    model.md_root = pe.Var(md_idx, bounds=(0, np.inf), domain=pe.Reals)
    for _m in inputs.free_rows:
        cid = inputs.cluster_ids[_m]
        for _n in params.tumor_clones:
            model.constraints.add(
                model.cA[_m, _n] - model.cA[_m, 0] <= model.md_root[_m, _n, "a"]
            )
            model.constraints.add(
                model.cA[_m, 0] - model.cA[_m, _n] <= model.md_root[_m, _n, "a"]
            )
            model.constraints.add(
                model.cB[_m, _n] - model.cB[_m, 0] <= model.md_root[_m, _n, "b"]
            )
            model.constraints.add(
                model.cB[_m, 0] - model.cB[_m, _n] <= model.md_root[_m, _n, "b"]
            )
            obj += inputs.w[cid] * (
                model.md_root[_m, _n, "a"] + model.md_root[_m, _n, "b"]
            )
    for _m in inputs.fixed_rows:
        cid = inputs.cluster_ids[_m]
        ca, cb = inputs.copy_numbers[cid]
        obj += inputs.w[cid] * (n - 1) * (abs(ca - params.base) + abs(cb - params.base))
    return obj


def build_reg_dadj_sum(model, mode: str, params: SolverParams, inputs: SolverInputs):
    """DADJ_SUM: pairwise L1 distance between all clone pairs."""
    n = params.n
    obj = 0
    md_idx = [
        (_m, _n1, _n2, ab)
        for _m in inputs.free_rows
        for _n1 in range(n - 1)
        for _n2 in range(_n1 + 1, n)
        for ab in ("a", "b")
    ]
    model.md_adj = pe.Var(md_idx, bounds=(0, np.inf), domain=pe.Reals)
    for _m in inputs.free_rows:
        cid = inputs.cluster_ids[_m]
        for _n1 in range(n - 1):
            for _n2 in range(_n1 + 1, n):
                model.constraints.add(
                    model.cA[_m, _n1] - model.cA[_m, _n2]
                    <= model.md_adj[_m, _n1, _n2, "a"]
                )
                model.constraints.add(
                    model.cA[_m, _n2] - model.cA[_m, _n1]
                    <= model.md_adj[_m, _n1, _n2, "a"]
                )
                model.constraints.add(
                    model.cB[_m, _n1] - model.cB[_m, _n2]
                    <= model.md_adj[_m, _n1, _n2, "b"]
                )
                model.constraints.add(
                    model.cB[_m, _n2] - model.cB[_m, _n1]
                    <= model.md_adj[_m, _n1, _n2, "b"]
                )
                obj += inputs.w[cid] * (
                    model.md_adj[_m, _n1, _n2, "a"] + model.md_adj[_m, _n1, _n2, "b"]
                )
    for _m in inputs.fixed_rows:
        cid = inputs.cluster_ids[_m]
        ca, cb = inputs.copy_numbers[cid]
        obj += inputs.w[cid] * (n - 1) * (abs(ca - params.base) + abs(cb - params.base))
    return obj


def build_reg_dbox_l1(model, mode: str, params: SolverParams, inputs: SolverInputs):
    """DBOX_L1: penalise CN range (max - min) per allele per cluster."""
    obj = 0
    sp_idx = [
        (_m, tag) for _m in inputs.free_rows for tag in ("maxA", "minA", "maxB", "minB")
    ]
    model.span = pe.Var(sp_idx, bounds=(0, np.inf), domain=pe.Reals)
    for _m in inputs.free_rows:
        for _n in params.tumor_clones:
            model.constraints.add(model.cA[_m, _n] <= model.span[_m, "maxA"])
            model.constraints.add(model.cA[_m, _n] >= model.span[_m, "minA"])
            model.constraints.add(model.cB[_m, _n] <= model.span[_m, "maxB"])
            model.constraints.add(model.cB[_m, _n] >= model.span[_m, "minB"])
        cid = inputs.cluster_ids[_m]
        obj += inputs.w[cid] * (
            model.span[_m, "maxA"]
            - model.span[_m, "minA"]
            + model.span[_m, "maxB"]
            - model.span[_m, "minB"]
        )
    return obj


def build_reg_dbox_l0(model, mode: str, params: SolverParams, inputs: SolverInputs):
    """DBOX_L0: L0((a_max - a_min) + (b_max - b_min)) per cluster.

    Same span as DBOX_L1 but binary: 1 if any allelic span > 0, 0 otherwise.
    Uses the same max/min auxiliary vars, then big-M links span to a binary indicator.
    """
    obj = 0
    big_M = 2 * params.cn_max
    sp_idx = [
        (_m, tag) for _m in inputs.free_rows for tag in ("maxA", "minA", "maxB", "minB")
    ]
    model.span_l0 = pe.Var(sp_idx, bounds=(0, np.inf), domain=pe.Reals)
    model.is_subclonal = pe.Var(inputs.free_rows, bounds=(0, 1), domain=pe.Binary)

    for _m in inputs.free_rows:
        for _n in params.tumor_clones:
            model.constraints.add(model.cA[_m, _n] <= model.span_l0[_m, "maxA"])
            model.constraints.add(model.cA[_m, _n] >= model.span_l0[_m, "minA"])
            model.constraints.add(model.cB[_m, _n] <= model.span_l0[_m, "maxB"])
            model.constraints.add(model.cB[_m, _n] >= model.span_l0[_m, "minB"])

        # span = (maxA - minA) + (maxB - minB)
        # is_subclonal = 1 iff span > 0
        # Linearised: span <= big_M * is_subclonal  AND  span >= is_subclonal (if span>0 then >=1)
        span_expr = (
            model.span_l0[_m, "maxA"]
            - model.span_l0[_m, "minA"]
            + model.span_l0[_m, "maxB"]
            - model.span_l0[_m, "minB"]
        )
        model.constraints.add(span_expr <= big_M * model.is_subclonal[_m])
        model.constraints.add(span_expr >= model.is_subclonal[_m])

        cid = inputs.cluster_ids[_m]
        obj += inputs.w[cid] * model.is_subclonal[_m]

    return obj


def build_regularization(model, mode: str, params: SolverParams, inputs: SolverInputs):
    """Dispatch to the appropriate regularization builder."""
    pname = params.reg_name
    model.pparam = pe.Param(mutable=True, initialize=0.0)

    obj_reg = 0
    if mode not in ("FULL", "CARCH") or pname == "RAW":
        return obj_reg

    builders = {
        "MAXCN": build_reg_maxcn,
        "DROOT_SUM": build_reg_droot_sum,
        "DADJ_SUM": build_reg_dadj_sum,
        "DBOX_L1": build_reg_dbox_l1,
        "DBOX_L0": build_reg_dbox_l0,
    }
    if pname in builders:
        obj_reg = builders[pname](model, mode, params, inputs)

    return obj_reg
