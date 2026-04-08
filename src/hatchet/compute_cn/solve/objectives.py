"""Pyomo objective builders for the copy-number ILP model."""

from __future__ import annotations

from pyomo import environ as pe

from hatchet.compute_cn.solve.variables import SolverParams


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
