"""Pyomo variable creation and solver parameter container."""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from pyomo import environ as pe


@dataclass
class SolverParams:
    """Immutable snapshot of solver configuration, passed to all build functions."""

    m: int
    n: int
    k: int
    cn_max: int
    mode: str  # "FULL", "CARCH", "UARCH"
    base: int  # 1 (diploid) or 2 (tetraploid)
    free_rows: list
    fixed_rows: set
    copy_numbers: dict
    cluster_ids: list
    sample_ids: list
    w: pd.Series
    ampdel: bool
    minprop: float
    max_ncns_seg: int
    purities: dict | None
    mrca: bool
    max_degree: int
    balanced_clusters: list | None
    tol: float
    zero_cn_thres: float
    f_a: pd.DataFrame = field(repr=False)
    f_b: pd.DataFrame = field(repr=False)

    @property
    def M(self):
        return math.floor(math.log2(self.cn_max)) + 1

    @property
    def tumor_clones(self):
        return range(1, self.n)

    def symmCoeff(self, i):
        return i + 1


def build_variables(model, p: SolverParams):
    """Create all indexed Pyomo Vars on *model*. Fix constants via .fix()."""
    _M = p.M

    # L1 slack and predicted fractional CN
    model.yA = pe.Var(range(p.m), range(p.k), bounds=(0, np.inf), domain=pe.Reals)
    model.yB = pe.Var(range(p.m), range(p.k), bounds=(0, np.inf), domain=pe.Reals)
    model.fA = pe.Var(range(p.m), range(p.k), bounds=(0, p.cn_max), domain=pe.Reals)
    model.fB = pe.Var(range(p.m), range(p.k), bounds=(0, p.cn_max), domain=pe.Reals)

    # Integer allele-specific CN
    if p.mode in ("FULL", "CARCH"):
        model.cA = pe.Var(
            range(p.m), range(p.n), bounds=(0, p.cn_max), domain=pe.Integers
        )
        model.cB = pe.Var(
            range(p.m), range(p.n), bounds=(0, p.cn_max), domain=pe.Integers
        )
        model.cA[:, 0].fix(1)
        model.cB[:, 0].fix(1)
        for _m in p.fixed_rows:
            ca, cb = p.copy_numbers[p.cluster_ids[_m]]
            for _n in p.tumor_clones:
                model.cA[_m, _n].fix(ca)
                model.cB[_m, _n].fix(cb)
        if p.ampdel:
            model.adA = pe.Var(p.free_rows, bounds=(0, 1), domain=pe.Binary)
            model.adB = pe.Var(p.free_rows, bounds=(0, 1), domain=pe.Binary)

    # Bit encoding
    if (p.mode == "FULL") or (p.max_ncns_seg > 0 and p.mode == "CARCH"):
        bit_idx = [
            (_b, _m, _n)
            for _b in range(_M)
            for _m in p.free_rows
            for _n in p.tumor_clones
        ]
        model.bitcA = pe.Var(bit_idx, bounds=(0, 1), domain=pe.Binary)
        model.bitcB = pe.Var(bit_idx, bounds=(0, 1), domain=pe.Binary)

    # Clone proportions
    if p.mode in ("FULL", "UARCH"):
        model.u = pe.Var(range(p.n), range(p.k), bounds=(0, 1), domain=pe.Reals)
        if p.purities is not None:
            for _k, sid in enumerate(p.sample_ids):
                if sid in p.purities:
                    model.u[0, _k].fix(1 - p.purities[sid])

    # McCormick product variables
    if p.mode == "FULL":
        v_idx = [
            (_b, _m, _n, _k)
            for _b in range(_M)
            for _m in p.free_rows
            for _n in p.tumor_clones
            for _k in range(p.k)
        ]
        model.vA = pe.Var(v_idx, bounds=(0, 1), domain=pe.Reals)
        model.vB = pe.Var(v_idx, bounds=(0, 1), domain=pe.Reals)

    # Minprop indicators
    if p.mode in ("FULL", "UARCH") and p.minprop > 0:
        model.x = pe.Var(range(1, p.n), range(p.k), domain=pe.Binary)

    # Max-NCNS-Seg group assignments
    if p.mode in ("FULL", "CARCH") and p.max_ncns_seg > 0:
        z_idx = [
            (_m, _n, _d)
            for _m in p.free_rows
            for _n in p.tumor_clones
            for _d in range(p.max_ncns_seg)
        ]
        model.z = pe.Var(z_idx, bounds=(0, 1), domain=pe.Binary)
