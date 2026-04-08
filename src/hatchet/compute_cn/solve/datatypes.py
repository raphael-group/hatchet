"""Solver dataclasses: SolverParams and SolverInputs."""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import pandas as pd


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

    @property
    def M(self):
        return math.floor(math.log2(self.cn_max)) + 1

    @property
    def tumor_clones(self):
        return range(1, self.n)

    def symmCoeff(self, i):
        return i + 1
