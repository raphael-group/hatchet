"""Pyomo variable creation and solver parameter/input containers."""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from pyomo import environ as pe


@dataclass
class SolverInputs:
    """Data inputs — observed FCN, weights, cluster structure."""

    f_a: pd.DataFrame = field(repr=False)
    f_b: pd.DataFrame = field(repr=False)
    w: pd.Series
    cluster_ids: list
    sample_ids: list
    copy_numbers: dict
    free_rows: list
    fixed_rows: set
    purities: dict | None
    balanced_clusters: list | None

    @property
    def m(self):
        return len(self.cluster_ids)

    @property
    def k(self):
        return len(self.sample_ids)


@dataclass
class SolverParams:
    """Solver configuration — clone count, CN bounds, regularization."""

    n: int
    cn_max: int
    base: int  # 1 (diploid) or 2 (tetraploid)
    ampdel: bool
    minprop: float
    max_ncns_seg: int
    mrca: bool
    max_degree: int
    tol: float
    zero_cn_thres: float
    reg_name: str = "RAW"
    reg_lambda: float = 0.0

    @property
    def M(self):
        return math.floor(math.log2(self.cn_max)) + 1

    @property
    def tumor_clones(self):
        return range(1, self.n)

    def symmCoeff(self, i):
        return i + 1


def build_variables(model, mode: str, params: SolverParams, inputs: SolverInputs):
    """Create all indexed Pyomo Vars on *model*. Fix constants via .fix()."""
    m, n, k = inputs.m, params.n, inputs.k
    _M = params.M

    model.yA = pe.Var(range(m), range(k), bounds=(0, np.inf), domain=pe.Reals)
    model.yB = pe.Var(range(m), range(k), bounds=(0, np.inf), domain=pe.Reals)
    model.fA = pe.Var(range(m), range(k), bounds=(0, params.cn_max), domain=pe.Reals)
    model.fB = pe.Var(range(m), range(k), bounds=(0, params.cn_max), domain=pe.Reals)

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
