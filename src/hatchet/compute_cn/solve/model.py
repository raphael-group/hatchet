"""Model builder: model assembly and initialization helpers."""

from __future__ import annotations

import numpy as np
from pyomo import environ as pe

from hatchet.compute_cn.solve.datatypes import SolverParams, SolverInputs
from hatchet.compute_cn.solve.constraints import (
    add_l1_constraints,
    add_ci_hinge_constraints,
    add_mixture_constraints,
    add_bit_encoding,
    add_proportion_constraints,
    add_domain_constraints,
    add_ncns_seg_constraints,
    update_fixed_u,  # noqa: F401
    update_fixed_cn,  # noqa: F401
)
from hatchet.compute_cn.solve.regularization import build_regularization
from hatchet.compute_cn.solve.objectives import (
    build_imf_objective,
    build_ci_violation_objective,
    build_final_objective,
)


def first_hot_start(params: SolverParams, inputs: SolverInputs):
    """Generate initial CN guess by rounding observed fractional CN."""
    if params.max_ncns_seg > 0:
        targetA = np.empty((inputs.m, params.max_ncns_seg))
        targetB = np.empty_like(targetA)
        for _m in range(inputs.m):
            targetA[_m, :] = np.random.choice(
                inputs.f_a.iloc[_m].values, params.max_ncns_seg, replace=False
            )
            targetB[_m, :] = np.random.choice(
                inputs.f_b.iloc[_m].values, params.max_ncns_seg, replace=False
            )
        targetA = targetA.round()
        targetB = targetB.round()
    else:
        targetA = np.round(inputs.f_a.values)
        targetB = np.round(inputs.f_b.values)

    hcA = np.zeros((inputs.m, params.n))
    hcB = np.zeros((inputs.m, params.n))

    for _m, cid in enumerate(inputs.cluster_ids):
        adA = adB = 0
        hcA[_m][0] = 1
        hcB[_m][0] = 1
        for _n in range(1, params.n):
            if cid in inputs.copy_numbers:
                hcA[_m][_n] = inputs.copy_numbers[cid][0]
                hcB[_m][_n] = inputs.copy_numbers[cid][1]
                continue
            mod = min(params.n, inputs.k)
            a = min(targetA[_m][_n % mod], params.cn_max)
            b = min(targetB[_m][_n % mod], params.cn_max)
            if params.ampdel:
                base = params.base
                if adA == 0 and a > base:
                    adA = 1
                if adA == 0 and a < base:
                    adA = -1
                if adB == 0 and b > base:
                    adB = 1
                if adB == 0 and b < base:
                    adB = -1
                a = max(a, base) if adA >= 0 else min(a, base)
                b = max(b, base) if adB >= 0 else min(b, base)
                if a + b > params.cn_max:
                    a = params.cn_max - base
                    b = params.cn_max - a
            hcA[_m][_n] = a
            hcB[_m][_n] = b if a + b <= params.cn_max else max(params.cn_max - a, 0)
    return hcA, hcB


def hot_start(model, params: SolverParams, inputs: SolverInputs, _cA=None, _cB=None):
    """Set warm-start values on model.cA/cB from a CN guess."""
    if _cA is None:
        _cA, _cB = first_hot_start(params, inputs)
    m, n = len(_cA), len(_cA[0])

    rank = np.zeros(n)
    rank[0] = -1
    for _m in range(m):
        for _n in range(1, n):
            rank[_n] += (_cA[_m][_n] + _cB[_m][_n]) * params.symmCoeff(_m)
    rank_indices = np.argsort(rank)

    for _m in range(m):
        if _m in inputs.fixed_rows:
            continue
        for _n in range(1, n):
            model.cA[_m, rank_indices[_n]].value = _cA[_m][_n]
            model.cB[_m, rank_indices[_n]].value = _cB[_m][_n]


def _random_tumor_props_dirichlet(n_tumor, alpha, minprop):
    t = np.random.dirichlet(alpha * np.ones(n_tumor))
    t[t < minprop] = 0
    if t.sum() > 0:
        t = t / t.sum()
    else:
        t = np.zeros(n_tumor)
        t[np.random.randint(n_tumor)] = 1.0
    return t


def _random_tumor_props_bubble(n_tumor, minprop, size_bubbles=100):
    t = np.zeros(n_tumor)
    n0 = np.random.randint(1, n_tumor + 1)
    n1 = np.random.randint(1, n_tumor + 1)
    n_parts = min(max(n0, n1), min(size_bubbles, n_tumor))
    positions = np.random.choice(n_tumor, n_parts, replace=False)
    breaks = (
        np.sort(
            np.random.choice(np.arange(1, size_bubbles), n_parts - 1, replace=False)
        )
        / size_bubbles
    )
    t[positions] = np.diff(breaks, prepend=0, append=1)
    t[(minprop - 1e-6 <= t) & (t < minprop)] = minprop
    t[t < minprop] = 0
    if t.sum() > 0:
        t = t / t.sum()
    else:
        t = np.zeros(n_tumor)
        t[np.random.randint(n_tumor)] = 1.0
    return t


def build_random_u(
    params: SolverParams, inputs: SolverInputs, method="dirichlet", alpha=0.3
):
    """Generate random U initialization matrix (n * k)."""
    U = np.empty((params.n, inputs.k))
    n_tumor = params.n - 1
    for _k in range(inputs.k):
        sid = inputs.sample_ids[_k]
        if inputs.purities is not None and sid in inputs.purities:
            purity = inputs.purities[sid]
            U[0, _k] = 1 - purity
            if n_tumor == 1:
                U[1, _k] = purity
            else:
                if method == "bubble":
                    t = _random_tumor_props_bubble(n_tumor, params.minprop)
                else:
                    t = _random_tumor_props_dirichlet(n_tumor, alpha, params.minprop)
                U[1:, _k] = purity * t
        else:
            if method == "bubble":
                t = _random_tumor_props_bubble(params.n, params.minprop)
            else:
                t = _random_tumor_props_dirichlet(params.n, alpha, params.minprop)
            U[:, _k] = t
    return U


def _build_variables(model, mode: str, params: SolverParams, inputs: SolverInputs):
    """Create all indexed Pyomo Vars on *model*. Fix constants via .fix()."""
    m, n, k = inputs.m, params.n, inputs.k
    _M = params.M

    model.fA = pe.Var(range(m), range(k), bounds=(0, params.cn_max), domain=pe.Reals)
    model.fB = pe.Var(range(m), range(k), bounds=(0, params.cn_max), domain=pe.Reals)
    if params.obj_type != "ci":
        model.yA = pe.Var(range(m), range(k), bounds=(0, np.inf), domain=pe.Reals)
        model.yB = pe.Var(range(m), range(k), bounds=(0, np.inf), domain=pe.Reals)
    if params.obj_type != "imf":
        model.hA = pe.Var(range(m), range(k), bounds=(0, np.inf), domain=pe.Reals)
        model.hB = pe.Var(range(m), range(k), bounds=(0, np.inf), domain=pe.Reals)

    if mode in ("FULL", "CARCH"):
        model.cA = pe.Var(
            range(m), range(n), bounds=(0, params.cn_max), domain=pe.Integers
        )
        model.cB = pe.Var(
            range(m), range(n), bounds=(0, params.cn_max), domain=pe.Integers
        )
        model.cA[:, 0].fix(1)
        model.cB[:, 0].fix(1)
        for _m in inputs.fixed_rows:
            ca, cb = inputs.copy_numbers[inputs.cluster_ids[_m]]
            for _n in params.tumor_clones:
                model.cA[_m, _n].fix(ca)
                model.cB[_m, _n].fix(cb)
        if params.ampdel:
            model.adA = pe.Var(inputs.free_rows, bounds=(0, 1), domain=pe.Binary)
            model.adB = pe.Var(inputs.free_rows, bounds=(0, 1), domain=pe.Binary)

    if (mode == "FULL") or (params.max_ncns_seg > 0 and mode == "CARCH"):
        bit_idx = [
            (_b, _m, _n)
            for _b in range(_M)
            for _m in inputs.free_rows
            for _n in params.tumor_clones
        ]
        model.bitcA = pe.Var(bit_idx, bounds=(0, 1), domain=pe.Binary)
        model.bitcB = pe.Var(bit_idx, bounds=(0, 1), domain=pe.Binary)

    if mode in ("FULL", "UARCH"):
        model.u = pe.Var(range(n), range(k), bounds=(0, 1), domain=pe.Reals)
        if inputs.purities is not None:
            for _k, sid in enumerate(inputs.sample_ids):
                if sid in inputs.purities:
                    model.u[0, _k].fix(1 - inputs.purities[sid])

    if mode == "FULL":
        v_idx = [
            (_b, _m, _n, _k)
            for _b in range(_M)
            for _m in inputs.free_rows
            for _n in params.tumor_clones
            for _k in range(k)
        ]
        model.vA = pe.Var(v_idx, bounds=(0, 1), domain=pe.Reals)
        model.vB = pe.Var(v_idx, bounds=(0, 1), domain=pe.Reals)

    if mode in ("FULL", "UARCH") and params.minprop > 0:
        model.x = pe.Var(range(1, n), range(k), domain=pe.Binary)

    if mode in ("FULL", "CARCH") and params.max_ncns_seg > 0:
        z_idx = [
            (_m, _n, _d)
            for _m in inputs.free_rows
            for _n in params.tumor_clones
            for _d in range(params.max_ncns_seg)
        ]
        model.z = pe.Var(z_idx, bounds=(0, 1), domain=pe.Binary)


def build_model(
    mode: str,
    params: SolverParams,
    inputs: SolverInputs,
    fixed_u=None,
    fixed_cA=None,
    fixed_cB=None,
):
    """Build a complete Pyomo ConcreteModel."""
    model = pe.ConcreteModel()
    _build_variables(model, mode, params, inputs)
    model.constraints = pe.ConstraintList()

    if params.obj_type != "ci":
        add_l1_constraints(model, mode, params, inputs)
    if params.obj_type != "imf":
        add_ci_hinge_constraints(model, mode, params, inputs)
    add_mixture_constraints(model, mode, params, inputs, fixed_u, fixed_cA, fixed_cB)
    add_bit_encoding(model, mode, params, inputs)
    add_proportion_constraints(model, mode, params, inputs)
    add_domain_constraints(model, mode, params, inputs)
    add_ncns_seg_constraints(model, mode, params, inputs)

    if params.obj_type == "ci":
        obj_imf = build_ci_violation_objective(model, mode, params, inputs)
    else:
        obj_imf = build_imf_objective(model, mode, params, inputs)
    obj_reg = build_regularization(model, mode, params, inputs)
    build_final_objective(model, obj_imf, obj_reg, params)

    if mode == "FULL":
        hot_start(model, params, inputs)

    return model
