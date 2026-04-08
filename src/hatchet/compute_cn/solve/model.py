"""Model builder: assembles Pyomo model from variables, constraints, and objectives."""

from __future__ import annotations

import numpy as np
from pyomo import environ as pe

# Stack-based random seeding for reproducibility
_random_states = []


class Random:
    """Context manager that pushes/pops numpy random state for reproducibility."""

    def __init__(self, seed=None):
        self.seed = seed

    def __enter__(self):
        if self.seed is not None:
            _random_states.append(np.random.get_state())
            np.random.set_state(np.random.RandomState(self.seed).get_state())

    def __exit__(self, *args):
        if self.seed is not None:
            np.random.set_state(_random_states.pop())


from hatchet.compute_cn.solve.variables import SolverParams, build_variables
from hatchet.compute_cn.solve.constraints import (
    add_l1_constraints,
    add_mixture_constraints,
    add_bit_encoding,
    add_proportion_constraints,
    add_domain_constraints,
    add_ncns_seg_constraints,
)
from hatchet.compute_cn.solve.regularization import build_regularization
from hatchet.compute_cn.solve.objectives import (
    build_imf_objective,
    build_final_objective,
)


def first_hot_start(p: SolverParams):
    """Generate initial CN guess by rounding observed fractional CN."""
    if p.max_ncns_seg > 0:
        targetA = np.empty((p.m, p.max_ncns_seg))
        targetB = np.empty_like(targetA)
        for _m in range(p.m):
            targetA[_m, :] = np.random.choice(
                p.f_a.iloc[_m].values, p.max_ncns_seg, replace=False
            )
            targetB[_m, :] = np.random.choice(
                p.f_b.iloc[_m].values, p.max_ncns_seg, replace=False
            )
        targetA = targetA.round()
        targetB = targetB.round()
    else:
        targetA = np.round(p.f_a.values)
        targetB = np.round(p.f_b.values)

    hcA = np.zeros((p.m, p.n))
    hcB = np.zeros((p.m, p.n))

    for _m, cid in enumerate(p.cluster_ids):
        adA = adB = 0
        hcA[_m][0] = 1
        hcB[_m][0] = 1
        for _n in range(1, p.n):
            if cid in p.copy_numbers:
                hcA[_m][_n] = p.copy_numbers[cid][0]
                hcB[_m][_n] = p.copy_numbers[cid][1]
                continue
            mod = min(p.n, p.k)
            a = min(targetA[_m][_n % mod], p.cn_max)
            b = min(targetB[_m][_n % mod], p.cn_max)
            if p.ampdel:
                base = p.base
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
                if a + b > p.cn_max:
                    a = p.cn_max - base
                    b = p.cn_max - a
            hcA[_m][_n] = a
            hcB[_m][_n] = b if a + b <= p.cn_max else max(p.cn_max - a, 0)
    return hcA, hcB


def hot_start(model, p: SolverParams, _cA=None, _cB=None):
    """Set warm-start values on model.cA/cB from a CN guess."""
    if _cA is None:
        _cA, _cB = first_hot_start(p)
    m, n = len(_cA), len(_cA[0])

    rank = np.zeros(n)
    rank[0] = -1
    for _m in range(m):
        for _n in range(1, n):
            rank[_n] += (_cA[_m][_n] + _cB[_m][_n]) * p.symmCoeff(_m)
    rank_indices = np.argsort(rank)

    for _m in range(m):
        if _m in p.fixed_rows:
            continue
        for _n in range(1, n):
            model.cA[_m, rank_indices[_n]].value = _cA[_m][_n]
            model.cB[_m, rank_indices[_n]].value = _cB[_m][_n]


def build_random_u(p: SolverParams, method="dirichlet", alpha=0.3):
    """Generate random U initialization matrix (n × k)."""
    U = np.empty((p.n, p.k))
    n_tumor = p.n - 1
    for _k in range(p.k):
        sid = p.sample_ids[_k]
        if p.purities is not None and sid in p.purities:
            purity = p.purities[sid]
            U[0, _k] = 1 - purity
            if n_tumor == 1:
                U[1, _k] = purity
            else:
                t = np.random.dirichlet(alpha * np.ones(n_tumor))
                t[t < p.minprop] = 0
                if t.sum() > 0:
                    t = t / t.sum()
                else:
                    t = np.zeros(n_tumor)
                    t[np.random.randint(n_tumor)] = 1.0
                U[1:, _k] = purity * t
        else:
            t = np.random.dirichlet(alpha * np.ones(p.n))
            t[t < p.minprop] = 0
            if t.sum() > 0:
                t = t / t.sum()
            else:
                t = np.zeros(p.n)
                t[0] = 1.0
            U[:, _k] = t
    return U


def build_model(
    params: SolverParams, penalty_param, fixed_u=None, fixed_cA=None, fixed_cB=None
):
    """Build a complete Pyomo ConcreteModel for the CN deconvolution ILP.

    Returns:
        (model, var_z): Pyomo model and DRMST topology vars (None if not DRMST).
    """
    model = pe.ConcreteModel()
    build_variables(model, params)
    model.constraints = pe.ConstraintList()

    add_l1_constraints(model, params)
    add_mixture_constraints(model, params, fixed_u, fixed_cA, fixed_cB)
    add_bit_encoding(model, params)
    add_proportion_constraints(model, params)
    add_domain_constraints(model, params)
    add_ncns_seg_constraints(model, params)

    obj_imf = build_imf_objective(model, params)
    obj_reg, var_z = build_regularization(model, params, penalty_param)
    build_final_objective(model, obj_imf, obj_reg, penalty_param)

    if params.mode == "FULL":
        hot_start(model, params)

    return model, var_z
