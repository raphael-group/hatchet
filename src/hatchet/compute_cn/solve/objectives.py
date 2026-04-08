"""Pyomo objective builders for the copy-number ILP model."""

from __future__ import annotations

from pyomo import environ as pe

from hatchet.compute_cn.solve.datatypes import SolverParams, SolverInputs


def build_imf_objective(model, mode: str, params: SolverParams, inputs: SolverInputs):
    """Weighted L1 deviation: Σ_{m,k} (yA + yB) * w[m]."""
    obj = 0
    for _m in range(inputs.m):
        cid = inputs.cluster_ids[_m]
        for _k in range(inputs.k):
            obj += (model.yA[_m, _k] + model.yB[_m, _k]) * inputs.w[cid]
    return obj


def build_final_objective(model, obj_imf, obj_reg, params: SolverParams):
    """Set model.obj combining IMF and regularization."""
    pname = params.reg_name
    model.obj_imf = pe.Expression(expr=obj_imf)
    model.obj_reg = pe.Expression(expr=obj_reg)
    if pname == "RAW":
        model.obj = pe.Objective(expr=obj_imf, sense=pe.minimize)
    else:
        model.obj = pe.Objective(
            expr=(1 - model.pparam) * obj_imf + model.pparam * obj_reg,
            sense=pe.minimize,
        )
