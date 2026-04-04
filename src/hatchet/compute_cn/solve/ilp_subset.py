import textwrap
import math
import logging
import numpy as np
import pandas as pd
from pyomo import environ as pe

from hatchet.compute_cn.solve.base_solver import BaseSolver
from hatchet.compute_cn.solve.utils import Random


class ILPSubset(BaseSolver):
    """ILP solver for the copy-number deconvolution problem.

    Models each genomic cluster's copy numbers (cA, cB) across n clones and
    their mixture proportions (u) across k samples. The objective minimizes a
    weighted L1 deviation between observed fractional copy numbers (f_a, f_b)
    and the mixture-model predictions. Supports optional regularization terms
    (MAXCN, DSPAN, DROOT_SUM, DADJ_SUM, DRMST) and a coordinate-descent alternation
    mode (CARCH / UARCH) in which either c or u is held fixed.
    """

    def __init__(
        self,
        n: int,
        cn_max: int,
        max_ncns_seg: int,
        minprop: float,
        ampdel: bool,
        copy_numbers: dict,
        fcn_data: dict,
        w: pd.Series,
        purities: dict,
        penalty_param: list,
        base: int = 1,
        zero_cn_thres=0.005,
        tol=0.001,
        balanced_clusters=None,
        mrca=False,
        max_degree=3,
    ):
        # Deep-copy fa/fb so ILPSubset can mutate freely
        fcn_data = dict(fcn_data)
        fcn_data["fa"] = fcn_data["fa"].copy(deep=True)
        fcn_data["fb"] = fcn_data["fb"].copy(deep=True)

        super().__init__(
            n=n,
            cn_max=cn_max,
            fcn_data=fcn_data,
            w=w,
            copy_numbers=copy_numbers,
            ampdel=ampdel,
            base=base,
            minprop=minprop,
            zero_cn_thres=zero_cn_thres,
            tol=tol,
            balanced_clusters=balanced_clusters,
            mrca=mrca,
        )

        # ILPSubset-specific fields
        self.max_ncns_seg = max_ncns_seg
        self.penalty_param = penalty_param
        self.purities = purities
        self.max_degree = max_degree
        self._var_z = None  # DRMST topology vars (set during create_model)

        self.mode = "FULL"

        # Values we want to optimize for
        self.cA = [[np.nan for _ in range(self.n)] for _ in range(self.m)]
        self.cB = [[np.nan for _ in range(self.n)] for _ in range(self.m)]
        self.u = [[np.nan for _ in range(self.k)] for _ in range(self.n)]

        # Fixed values of cA/cB/u
        self._fixed_cA = [[np.nan for _ in range(self.n)] for _ in range(self.m)]
        self._fixed_cB = [[np.nan for _ in range(self.n)] for _ in range(self.m)]
        self._fixed_u = [[np.nan for _ in range(self.k)] for _ in range(self.n)]

    def __copy__(self):
        new = self.__class__.__new__(self.__class__)
        # Share read-only data from BaseSolver
        new.fcn_data = self.fcn_data
        new.f_a = self.f_a
        new.f_b = self.f_b
        new.m, new.k = self.m, self.k
        new.cluster_ids = self.cluster_ids
        new.sample_ids = self.sample_ids
        new.n = self.n
        new.cn_max = self.cn_max
        new._base = self._base
        new.w = self.w
        new.copy_numbers = self.copy_numbers
        new.ampdel = self.ampdel
        new.minprop = self.minprop
        new.zero_cn_thres = self.zero_cn_thres
        new.tol = self.tol
        new.balanced_clusters = self.balanced_clusters
        new.mrca = self.mrca
        new.warmstart = False
        new.model = None
        # ILPSubset-specific
        new.max_ncns_seg = self.max_ncns_seg
        new.penalty_param = self.penalty_param
        new.purities = self.purities
        new.max_degree = self.max_degree
        new._var_z = None
        new.mode = "FULL"
        new.cA = [[np.nan for _ in range(new.n)] for _ in range(new.m)]
        new.cB = [[np.nan for _ in range(new.n)] for _ in range(new.m)]
        new.u = [[np.nan for _ in range(new.k)] for _ in range(new.n)]
        new._fixed_cA = [[np.nan for _ in range(new.n)] for _ in range(new.m)]
        new._fixed_cB = [[np.nan for _ in range(new.n)] for _ in range(new.m)]
        new._fixed_u = [[np.nan for _ in range(new.k)] for _ in range(new.n)]
        return new

    def __str__(self):
        # Pyomo pprint gives us too much information - too unwieldy for large models
        # This method is implemented to supply the bare-minimum but useful model information.
        if self.model is None:
            return ""
        else:
            return textwrap.dedent(
                f"""
                # ------------------------------------------
                #   Problem Information
                # ------------------------------------------
                #     Number of constraints: {self.model.nconstraints()}
                #     Number of variables: {self.model.nvariables()}
                # ------------------------------------------
            """
            )

    @property
    def M(self):
        return math.floor(math.log2(self.cn_max)) + 1

    @property
    def optimized_cA(self):
        if self.mode == "UARCH":
            return self._fixed_cA
        else:
            return self.cA

    @property
    def optimized_cB(self):
        if self.mode == "UARCH":
            return self._fixed_cB
        else:
            return self.cB

    @property
    def optimized_u(self):
        if self.mode == "CARCH":
            return self._fixed_u
        else:
            return self.u

    def create_model(self, pprint=False):
        m, n, k = self.m, self.n, self.k
        f_a, f_b = self.f_a, self.f_b
        cn_max = self.cn_max
        ampdel = self.ampdel
        copy_numbers = self.copy_numbers
        mode_t = self.mode
        max_ncns_seg = self.max_ncns_seg
        _M = self.M  # compute binary length
        _base = self.base
        zero_cn_thres = self.zero_cn_thres

        model = pe.ConcreteModel()

        # Identify fixed-CN cluster rows: store constants, skip variable creation
        fixed_rows = set()
        if mode_t in ("FULL", "CARCH"):
            for _m, cid in enumerate(self.cluster_ids):
                if cid in copy_numbers:
                    fixed_rows.add(_m)
                    ca, cb = copy_numbers[cid]
                    self.cA[_m][0] = self._base
                    self.cB[_m][0] = self._base
                    for _n in range(1, n):
                        self.cA[_m][_n] = ca
                        self.cB[_m][_n] = cb
        self._fixed_rows = fixed_rows
        free_rows = [_m for _m in range(m) if _m not in fixed_rows]

        # 0-cn_max var per m*k
        fA = {}
        fB = {}

        # 0-inf var per m*k
        yA = {}
        yB = {}

        # 0-1 var per m, ampdel only
        adA = {}
        adB = {}

        # upper bound for solver
        cAB_bounds = {}
        for _m in range(m):
            cluster_id = f_a.index[_m]
            cAB_bounds[_m] = max(sum(copy_numbers.get(cluster_id, (0, 0))), cn_max)

        for _m, _k in np.ndindex((m, k)):
            yA[(_m, _k)] = pe.Var(bounds=(0, np.inf), domain=pe.Reals)
            model.add_component(f"yA_{_m + 1}_{_k + 1}", yA[(_m, _k)])
            yB[(_m, _k)] = pe.Var(bounds=(0, np.inf), domain=pe.Reals)
            model.add_component(f"yB_{_m + 1}_{_k + 1}", yB[(_m, _k)])
            fA[(_m, _k)] = pe.Var(bounds=(0, cAB_bounds[_m]), domain=pe.Reals)
            model.add_component(f"fA_{_m + 1}_{_k + 1}", fA[(_m, _k)])
            fB[(_m, _k)] = pe.Var(bounds=(0, cAB_bounds[_m]), domain=pe.Reals)
            model.add_component(f"fB_{_m + 1}_{_k + 1}", fB[(_m, _k)])

        if mode_t in ("FULL", "CARCH"):
            for _m, _n in np.ndindex((m, n)):
                if _m in fixed_rows:
                    continue
                self.cA[_m][_n] = pe.Var(bounds=(0, cAB_bounds[_m]), domain=pe.Integers)
                model.add_component(f"cA_{_m + 1}_{_n + 1}", self.cA[_m][_n])
                self.cB[_m][_n] = pe.Var(bounds=(0, cAB_bounds[_m]), domain=pe.Integers)
                model.add_component(f"cB_{_m + 1}_{_n + 1}", self.cB[_m][_n])

            if ampdel:
                for _m in range(m):
                    cluster_id = f_a.index[_m]
                    if cluster_id not in copy_numbers:
                        adA[_m] = pe.Var(bounds=(0, 1), domain=pe.Binary)
                        model.add_component(f"adA_{_m + 1}", adA[_m])
                        adB[_m] = pe.Var(bounds=(0, 1), domain=pe.Binary)
                        model.add_component(f"adB_{_m + 1}", adB[_m])

        # 0-1 var per M * m * n; M=floor(log(cn_max)) + 1
        bitcA = {}
        bitcB = {}
        if (mode_t == "FULL") or (max_ncns_seg > 0 and mode_t == "CARCH"):
            for _b, _m, _n in np.ndindex((_M, m, n)):
                if _m in fixed_rows:
                    continue
                bitcA[(_b, _m, _n)] = pe.Var(bounds=(0, 1), domain=pe.Binary)
                model.add_component(
                    f"bitcA_{_b + 1}_{_m + 1}_{_n + 1}",
                    bitcA[(_b, _m, _n)],
                )
                bitcB[(_b, _m, _n)] = pe.Var(bounds=(0, 1), domain=pe.Binary)
                model.add_component(
                    f"bitcB_{_b + 1}_{_m + 1}_{_n + 1}",
                    bitcB[(_b, _m, _n)],
                )

        if mode_t in ("FULL", "UARCH"):
            for _n, _k in np.ndindex((n, k)):
                self.u[_n][_k] = pe.Var(bounds=(0, 1), domain=pe.Reals)
                model.add_component(f"u_{_n + 1}_{_k + 1}", self.u[_n][_k])

        vA = {}
        vB = {}
        if mode_t == "FULL":
            for _b, _m, _n, _k in np.ndindex((_M, m, n, k)):
                if _m in fixed_rows:
                    continue
                vA[(_b, _m, _n, _k)] = pe.Var(bounds=(0, 1), domain=pe.Reals)
                model.add_component(
                    f"vA_{_b + 1}_{_m + 1}_{_n + 1}_{_k + 1}",
                    vA[(_b, _m, _n, _k)],
                )
                vB[(_b, _m, _n, _k)] = pe.Var(bounds=(0, 1), domain=pe.Reals)
                model.add_component(
                    f"vB_{_b + 1}_{_m + 1}_{_n + 1}_{_k + 1}",
                    vB[(_b, _m, _n, _k)],
                )

        x = {}
        if (mode_t in ("FULL", "UARCH")) and (self.minprop > 0):
            for _n, _k in np.ndindex((n, k)):
                x[(_n, _k)] = pe.Var(domain=pe.Binary)
                model.add_component(f"x_{_n + 1}_{_k + 1}", x[(_n, _k)])

        # buildOptionalVariables
        z = {}
        if (mode_t in ("FULL", "CARCH")) and max_ncns_seg > 0:
            for _m, _n, _d in np.ndindex((m, n, max_ncns_seg)):
                if _n != 0 and _m not in fixed_rows:
                    z[(_m, _n, _d)] = pe.Var(bounds=(0, 1), domain=pe.Binary)
                    model.add_component(
                        f"z_{_m + 1}_{_n + 1}_{_d + 1}", z[(_m, _n, _d)]
                    )

        # CONSTRAINTS
        model.constraints = pe.ConstraintList()

        self._add_l1_constraints(
            model,
            lambda _m, _k: fA[(_m, _k)],
            lambda _m, _k: fB[(_m, _k)],
            lambda _m, _k: yA[(_m, _k)],
            lambda _m, _k: yB[(_m, _k)],
        )

        if mode_t == "FULL":
            for _m, _k in np.ndindex((m, k)):
                if _m in fixed_rows:
                    # Fixed cluster: cA/cB are constants, use linear fA = Σ cA*u
                    sum_a = sum(self.cA[_m][_n] * self.u[_n][_k] for _n in range(n))
                    model.constraints.add(fA[(_m, _k)] == sum_a)
                else:
                    sum_a = 0
                    for _n, _b in np.ndindex((n, _M)):
                        sum_a += vA[(_b, _m, _n, _k)] * math.pow(2, _b)
                        model.constraints.add(
                            vA[(_b, _m, _n, _k)] <= bitcA[(_b, _m, _n)]
                        )
                        model.constraints.add(vA[(_b, _m, _n, _k)] <= self.u[_n][_k])
                        model.constraints.add(
                            vA[(_b, _m, _n, _k)]
                            >= bitcA[(_b, _m, _n)] + self.u[_n][_k] - 1
                        )
                    model.constraints.add(fA[(_m, _k)] == sum_a)

            for _m, _k in np.ndindex((m, k)):
                if _m in fixed_rows:
                    sum_b = sum(self.cB[_m][_n] * self.u[_n][_k] for _n in range(n))
                    model.constraints.add(fB[(_m, _k)] == sum_b)
                else:
                    sum_b = 0
                    for _n, _b in np.ndindex((n, _M)):
                        sum_b += vB[(_b, _m, _n, _k)] * math.pow(2, _b)
                        model.constraints.add(
                            vB[(_b, _m, _n, _k)] <= bitcB[(_b, _m, _n)]
                        )
                        model.constraints.add(vB[(_b, _m, _n, _k)] <= self.u[_n][_k])
                        model.constraints.add(
                            vB[(_b, _m, _n, _k)]
                            >= bitcB[(_b, _m, _n)] + self.u[_n][_k] - 1
                        )
                    model.constraints.add(fB[(_m, _k)] == sum_b)

            for _n, _k in np.ndindex((n, k)):
                _sum = 0
                for _m in range(m):
                    if _m in fixed_rows:
                        _sum += self.cA[_m][_n] + self.cB[_m][_n]
                    else:
                        for _b in range(_M):
                            _sum += bitcA[(_b, _m, _n)] + bitcB[(_b, _m, _n)]
                model.constraints.add(_sum >= self.u[_n][_k])

        if (mode_t == "FULL") or (max_ncns_seg > 0 and mode_t == "CARCH"):
            for _m, _n in np.ndindex((m, n)):
                if _m in fixed_rows:
                    continue
                sum_a = 0
                sum_b = 0
                for _b in range(_M):
                    sum_a += bitcA[(_b, _m, _n)] * math.pow(2, _b)
                    sum_b += bitcB[(_b, _m, _n)] * math.pow(2, _b)

                model.constraints.add(self.cA[_m][_n] == sum_a)
                model.constraints.add(self.cB[_m][_n] == sum_b)

        if mode_t == "CARCH":
            # TODO: These loops can be collapsed once validation against C++ is complete
            for _m, _k in np.ndindex((m, k)):
                _sumA = 0
                _sumB = 0
                for _n in range(n):
                    if self._fixed_u[_n][_k] >= self.minprop - self.tol:
                        _sumA += self.cA[_m][_n] * self._fixed_u[_n][_k]
                        _sumB += self.cB[_m][_n] * self._fixed_u[_n][_k]
                model.constraints.add(fA[(_m, _k)] == _sumA)
                model.constraints.add(fB[(_m, _k)] == _sumB)

        if mode_t in ("FULL", "CARCH"):
            get_cA = lambda _m, _n: self.cA[_m][_n]
            get_cB = lambda _m, _n: self.cB[_m][_n]
            self._add_cAB_upper_bound(model, get_cA, get_cB, rows=free_rows)
            self._add_normal_clone_constraints(model, get_cA, get_cB, rows=free_rows)
            self._add_zero_cn_constraints(model, get_cA, get_cB, rows=free_rows)

            if ampdel:
                for _m in free_rows:
                    cluster_id = f_a.index[_m]
                    if cluster_id not in copy_numbers:
                        for _n in range(1, n):
                            model.constraints.add(
                                self.cA[_m][_n]
                                <= cn_max * adA[_m] + _base - _base * adA[_m]
                            )
                            model.constraints.add(self.cA[_m][_n] >= _base * adA[_m])
                            model.constraints.add(
                                self.cB[_m][_n]
                                <= cn_max * adB[_m] + _base - _base * adB[_m]
                            )
                            model.constraints.add(self.cB[_m][_n] >= _base * adB[_m])

        if mode_t == "UARCH":
            # TODO: These loops can be collapsed once validation against C++ is complete
            for _m, _k in np.ndindex((m, k)):
                _sumA = 0
                _sumB = 0
                for _n in range(n):
                    _sumA += int(self._fixed_cA[_m][_n]) * self.u[_n][_k]
                    _sumB += int(self._fixed_cB[_m][_n]) * self.u[_n][_k]
                model.constraints.add(fA[(_m, _k)] == _sumA)
                model.constraints.add(fB[(_m, _k)] == _sumB)

        if mode_t in ("FULL", "UARCH"):
            for _k in range(k):
                _sum = sum(self.u[_n][_k] for _n in range(n))
                model.constraints.add(_sum == 1)

        if (mode_t in ("FULL", "UARCH")) and self.minprop > 0:
            for _k in range(k):
                for _n in range(1, n):
                    model.constraints.add(x[(_n, _k)] >= self.u[_n][_k])
                    model.constraints.add(self.u[_n][_k] >= self.minprop * x[(_n, _k)])

        # buildOptionalConstraints
        if (mode_t in ("FULL", "CARCH")) and max_ncns_seg > 0:
            for _m in free_rows:
                for _n in range(1, n):
                    _sum = 0
                    for _d in range(max_ncns_seg):
                        _sum += z[(_m, _n, _d)]
                    model.constraints.add(_sum == 1)

            for _m in free_rows:
                for _b, _d in np.ndindex((_M, max_ncns_seg)):
                    for _i in range(1, n - 1):
                        for _j in range(1, n):
                            model.constraints.add(
                                bitcA[(_b, _m, _i)] - bitcA[(_b, _m, _j)]
                                <= 2 - z[(_m, _i, _d)] - z[(_m, _j, _d)]
                            )
                            model.constraints.add(
                                bitcA[(_b, _m, _j)] - bitcA[(_b, _m, _i)]
                                <= 2 - z[(_m, _i, _d)] - z[(_m, _j, _d)]
                            )
                            model.constraints.add(
                                bitcB[(_b, _m, _i)] - bitcB[(_b, _m, _j)]
                                <= 2 - z[(_m, _i, _d)] - z[(_m, _j, _d)]
                            )
                            model.constraints.add(
                                bitcB[(_b, _m, _j)] - bitcB[(_b, _m, _i)]
                                <= 2 - z[(_m, _i, _d)] - z[(_m, _j, _d)]
                            )

            for _m in free_rows:
                for _d in range(max_ncns_seg - 1):
                    _sum_l = _sum_l1 = 0
                    for _n in range(1, self.n):
                        _sum_l += z[(_m, _n, _d)] * self.symmCoeff(_n)
                        _sum_l1 += z[(_m, _n, _d + 1)] * self.symmCoeff(_n)
                        model.constraints.add(_sum_l <= _sum_l1)  # FIXME wired term

        if mode_t in ("FULL", "CARCH"):
            get_cA = lambda _m, _n: self.cA[_m][_n]
            get_cB = lambda _m, _n: self.cB[_m][_n]
            # Symmetry breaking uses all rows (fixed rows contribute constants)
            self._add_symmetry_breaking(model, get_cA, get_cB)
            # Fixed CN constraints are no-ops for fixed rows (already constants)
            self._add_fixed_cn_constraints(model, get_cA, get_cB, rows=free_rows)
            self._add_balanced_constraints(model, get_cA, get_cB, rows=free_rows)
            self._add_mrca_loh_constraints(model, get_cA, get_cB, rows=free_rows)

        # Purities are hot-start seeds in CD (build_random_u), not constraints.

        # add objective & regularization terms
        # Format: IMF_loss + pparam * reg_term
        obj_imf = 0
        for _m, _k in np.ndindex((m, k)):
            obj_imf += (yA[(_m, _k)] + yB[(_m, _k)]) * self.w[self.cluster_ids[_m]]

        obj_reg = 0
        [pname, init_val] = self.penalty_param
        pparam = pe.Param(mutable=True, initialize=init_val)
        model.pparam = pparam

        if mode_t in ("FULL", "CARCH") and pname != "RAW":
            if pname == "MAXCN":
                hcn_vars = {}
                for _m in free_rows:
                    hcn_vars[(_m, "a")] = pe.Var(bounds=(0, np.inf), domain=pe.Reals)
                    model.add_component(f"HCA_{_m}", hcn_vars[(_m, "a")])
                    hcn_vars[(_m, "b")] = pe.Var(bounds=(0, np.inf), domain=pe.Reals)
                    model.add_component(f"HCB_{_m}", hcn_vars[(_m, "b")])
                    for _n in range(1, n):
                        model.constraints.add(self.cA[_m][_n] <= hcn_vars[(_m, "a")])
                        model.constraints.add(self.cB[_m][_n] <= hcn_vars[(_m, "b")])
                for _m in free_rows:
                    cluster_id = self.cluster_ids[_m]
                    obj_reg += self.w[cluster_id] * hcn_vars[(_m, "a")]
                    obj_reg += self.w[cluster_id] * hcn_vars[(_m, "b")]
                # Fixed rows: constant contribution = max CN across clones
                for _m in fixed_rows:
                    cluster_id = self.cluster_ids[_m]
                    ca, cb = copy_numbers[cluster_id]
                    obj_reg += self.w[cluster_id] * max(ca, self._base)
                    obj_reg += self.w[cluster_id] * max(cb, self._base)

            elif pname == "DROOT_SUM":
                manhat_vars = {}
                for _m in free_rows:
                    for _n in range(1, n):
                        manhat_vars[(_m, _n, "a")] = pe.Var(
                            bounds=(0, np.inf), domain=pe.Reals
                        )
                        model.add_component(
                            f"MDA_{_m}_{_n}", manhat_vars[(_m, _n, "a")]
                        )
                        manhat_vars[(_m, _n, "b")] = pe.Var(
                            bounds=(0, np.inf), domain=pe.Reals
                        )
                        model.add_component(
                            f"MDB_{_m}_{_n}", manhat_vars[(_m, _n, "b")]
                        )
                        model.constraints.add(
                            self.cA[_m][_n] - self.cA[_m][0]
                            <= manhat_vars[(_m, _n, "a")]
                        )
                        model.constraints.add(
                            self.cA[_m][0] - self.cA[_m][_n]
                            <= manhat_vars[(_m, _n, "a")]
                        )
                        model.constraints.add(
                            self.cB[_m][_n] - self.cB[_m][0]
                            <= manhat_vars[(_m, _n, "b")]
                        )
                        model.constraints.add(
                            self.cB[_m][0] - self.cB[_m][_n]
                            <= manhat_vars[(_m, _n, "b")]
                        )
                for _m in free_rows:
                    cluster_id = self.cluster_ids[_m]
                    for _n in range(1, n):
                        obj_reg += self.w[cluster_id] * manhat_vars[(_m, _n, "a")]
                        obj_reg += self.w[cluster_id] * manhat_vars[(_m, _n, "b")]
                # Fixed rows: constant |cA - base| + |cB - base| per clone
                for _m in fixed_rows:
                    cluster_id = self.cluster_ids[_m]
                    ca, cb = copy_numbers[cluster_id]
                    obj_reg += (
                        self.w[cluster_id]
                        * (n - 1)
                        * (abs(ca - self._base) + abs(cb - self._base))
                    )

            elif pname == "DADJ_SUM":
                manhat_vars = {}
                for _m in free_rows:
                    for _n1 in range(n - 1):
                        for _n2 in range(_n1 + 1, n):
                            manhat_vars[(_m, _n1, _n2, "a")] = pe.Var(
                                bounds=(0, np.inf), domain=pe.Reals
                            )
                            model.add_component(
                                f"MDA_{_m}_{_n1}_{_n2}",
                                manhat_vars[(_m, _n1, _n2, "a")],
                            )
                            manhat_vars[(_m, _n1, _n2, "b")] = pe.Var(
                                bounds=(0, np.inf), domain=pe.Reals
                            )
                            model.add_component(
                                f"MDB_{_m}_{_n1}_{_n2}",
                                manhat_vars[(_m, _n1, _n2, "b")],
                            )
                            model.constraints.add(
                                self.cA[_m][_n1] - self.cA[_m][_n2]
                                <= manhat_vars[(_m, _n1, _n2, "a")]
                            )
                            model.constraints.add(
                                self.cA[_m][_n2] - self.cA[_m][_n1]
                                <= manhat_vars[(_m, _n1, _n2, "a")]
                            )
                            model.constraints.add(
                                self.cB[_m][_n1] - self.cB[_m][_n2]
                                <= manhat_vars[(_m, _n1, _n2, "b")]
                            )
                            model.constraints.add(
                                self.cB[_m][_n2] - self.cB[_m][_n1]
                                <= manhat_vars[(_m, _n1, _n2, "b")]
                            )
                for _m in free_rows:
                    cluster_id = self.cluster_ids[_m]
                    for _n1 in range(n - 1):
                        for _n2 in range(_n1 + 1, n):
                            obj_reg += (
                                self.w[cluster_id] * manhat_vars[(_m, _n1, _n2, "a")]
                            )
                            obj_reg += (
                                self.w[cluster_id] * manhat_vars[(_m, _n1, _n2, "b")]
                            )
                # Fixed rows: all clones identical, so all adj distances = 0
                # except clone 0 (normal) vs clone _n: |ca - base| + |cb - base|
                for _m in fixed_rows:
                    cluster_id = self.cluster_ids[_m]
                    ca, cb = copy_numbers[cluster_id]
                    # (n-1) pairs involve clone 0: (0,1), (0,2), ..., (0,n-1)
                    obj_reg += (
                        self.w[cluster_id]
                        * (n - 1)
                        * (abs(ca - self._base) + abs(cb - self._base))
                    )

            elif pname == "DSPAN":
                span_vars = {}
                for _m in free_rows:
                    span_vars[(_m, "maxA")] = pe.Var(
                        bounds=(0, np.inf), domain=pe.Reals
                    )
                    model.add_component(f"SPA_MAX_{_m}", span_vars[(_m, "maxA")])
                    span_vars[(_m, "minA")] = pe.Var(
                        bounds=(0, np.inf), domain=pe.Reals
                    )
                    model.add_component(f"SPA_MIN_{_m}", span_vars[(_m, "minA")])
                    span_vars[(_m, "maxB")] = pe.Var(
                        bounds=(0, np.inf), domain=pe.Reals
                    )
                    model.add_component(f"SPB_MAX_{_m}", span_vars[(_m, "maxB")])
                    span_vars[(_m, "minB")] = pe.Var(
                        bounds=(0, np.inf), domain=pe.Reals
                    )
                    model.add_component(f"SPB_MIN_{_m}", span_vars[(_m, "minB")])
                    for _n in range(1, n):
                        model.constraints.add(
                            self.cA[_m][_n] <= span_vars[(_m, "maxA")]
                        )
                        model.constraints.add(
                            self.cA[_m][_n] >= span_vars[(_m, "minA")]
                        )
                        model.constraints.add(
                            self.cB[_m][_n] <= span_vars[(_m, "maxB")]
                        )
                        model.constraints.add(
                            self.cB[_m][_n] >= span_vars[(_m, "minB")]
                        )
                for _m in free_rows:
                    cluster_id = self.cluster_ids[_m]
                    obj_reg += self.w[cluster_id] * (
                        span_vars[(_m, "maxA")]
                        - span_vars[(_m, "minA")]
                        + span_vars[(_m, "maxB")]
                        - span_vars[(_m, "minB")]
                    )
                # Fixed rows: all clones have same CN → span = 0, no contribution

            elif pname == "DRMST":
                big_M_tree = 2 * cn_max
                is_mrca = self.mrca

                param_deg = pe.Param(initialize=float(self.max_degree))
                model.add_component("p_max_degree", param_deg)

                # Topology variables z[i,j]: binary edges in r-arborescence
                var_z = {}
                for i in range(1, n):
                    for j in range(i):
                        if is_mrca and i >= 2 and j == 0:
                            continue
                        var_z[(i, j)] = pe.Var(bounds=(0, 1), domain=pe.Binary)
                        model.add_component(f"tz_{i}_{j}", var_z[(i, j)])
                # One parent per non-root clone
                for i in range(1, n):
                    model.constraints.add(
                        sum(var_z[(i, j)] for j in range(i) if (i, j) in var_z) == 1
                    )

                # Degree constraints
                ch0 = [var_z[(i, 0)] for i in range(1, n) if (i, 0) in var_z]
                if ch0:
                    model.add_component(
                        "con_degree_0",
                        pe.Constraint(expr=sum(ch0) <= param_deg),
                    )
                for c in range(1, n):
                    ch = [var_z[(i, c)] for i in range(c + 1, n) if (i, c) in var_z]
                    if ch:
                        model.add_component(
                            f"con_degree_{c}",
                            pe.Constraint(expr=1 + sum(ch) <= param_deg),
                        )

                # Distance variables M[m,i] + big-M linearization
                var_M = {}
                tree_clones = range(2, n) if is_mrca else range(1, n)
                for _m in range(m):
                    for i in range(1, n):
                        for j in range(i):
                            if (i, j) not in var_z or i not in tree_clones:
                                continue
                            if (_m, i) not in var_M:
                                var_M[(_m, i)] = pe.Var(
                                    bounds=(0, None), domain=pe.Reals
                                )
                                model.add_component(f"tM_{_m}_{i}", var_M[(_m, i)])
                            dA = pe.Var(bounds=(0, None), domain=pe.Reals)
                            model.add_component(f"tdA_{_m}_{i}_{j}", dA)
                            dB = pe.Var(bounds=(0, None), domain=pe.Reals)
                            model.add_component(f"tdB_{_m}_{i}_{j}", dB)
                            cA_i, cA_j = get_cA(_m, i), get_cA(_m, j)
                            cB_i, cB_j = get_cB(_m, i), get_cB(_m, j)
                            zv = var_z[(i, j)]
                            model.constraints.add(
                                dA >= cA_i - cA_j - big_M_tree * (1 - zv)
                            )
                            model.constraints.add(
                                dA >= cA_j - cA_i - big_M_tree * (1 - zv)
                            )
                            model.constraints.add(
                                dB >= cB_i - cB_j - big_M_tree * (1 - zv)
                            )
                            model.constraints.add(
                                dB >= cB_j - cB_i - big_M_tree * (1 - zv)
                            )
                            model.constraints.add(
                                var_M[(_m, i)] >= dA + dB - big_M_tree * (1 - zv)
                            )

                # Dynamic LOH constraints — split into two loops:
                # 1) create binary LOH indicator vars (lA=1 iff cA=0)
                for _m in range(m):
                    if _m in fixed_rows:
                        continue
                    for _n in range(n):
                        lA = pe.Var(bounds=(0, 1), domain=pe.Binary)
                        model.add_component(f"tlA_{_m}_{_n}", lA)
                        lB = pe.Var(bounds=(0, 1), domain=pe.Binary)
                        model.add_component(f"tlB_{_m}_{_n}", lB)
                        model.constraints.add(get_cA(_m, _n) >= 1 - big_M_tree * lA)
                        model.constraints.add(get_cA(_m, _n) <= big_M_tree * (1 - lA))
                        model.constraints.add(get_cB(_m, _n) >= 1 - big_M_tree * lB)
                        model.constraints.add(get_cB(_m, _n) <= big_M_tree * (1 - lB))

                # 2) if parent lost allele, child must also have 0 on that allele
                for _m in range(m):
                    if _m in fixed_rows:
                        continue
                    for i in range(1, n):
                        for j in range(i):
                            if (i, j) not in var_z:
                                continue
                            lA_j = model.find_component(f"tlA_{_m}_{j}")
                            lB_j = model.find_component(f"tlB_{_m}_{j}")
                            model.constraints.add(
                                get_cA(_m, i)
                                <= cn_max * (1 - lA_j)
                                + big_M_tree * (1 - var_z[(i, j)])
                            )
                            model.constraints.add(
                                get_cB(_m, i)
                                <= cn_max * (1 - lB_j)
                                + big_M_tree * (1 - var_z[(i, j)])
                            )

                # Tree regularization objective
                for _m in range(m):
                    for i in tree_clones:
                        if (_m, i) in var_M:
                            obj_reg += self.w[self.cluster_ids[_m]] * var_M[(_m, i)]

                self._var_z = var_z

        if mode_t == "FULL":
            self.hot_start()

        # Store objective components for post-solve extraction
        model.obj_imf = pe.Expression(expr=obj_imf)
        model.obj_reg = pe.Expression(expr=obj_reg)

        # Objective: (1-pparam)*IMF + pparam*reg, pparam ∈ [0,1]
        if pname == "RAW":
            model.obj = pe.Objective(expr=obj_imf, sense=pe.minimize)
        else:
            model.obj = pe.Objective(
                expr=(1 - pparam) * obj_imf + pparam * obj_reg,
                sense=pe.minimize,
            )
        self.model = model

        if pprint:
            logging.info(str(self))

    def first_hot_start(self):
        """
        ran by cd init step and also ILP init step
        """
        if self.max_ncns_seg > 0:
            targetA = np.empty((self.m, self.max_ncns_seg))
            targetB = np.empty_like(targetA)
            for _m in range(self.m):
                targetA[_m, :] = np.random.choice(
                    self.f_a.iloc[_m].values, self.max_ncns_seg, replace=False
                )
                targetB[_m, :] = np.random.choice(
                    self.f_b.iloc[_m].values, self.max_ncns_seg, replace=False
                )
            targetA = targetA.round()
            targetB = targetB.round()
        else:
            targetA = np.round(self.f_a.values)
            targetB = np.round(self.f_b.values)

        hcA = np.zeros((self.m, self.n))
        hcB = np.zeros((self.m, self.n))

        for _m, cluster_id in enumerate(self.cluster_ids):
            adA = 0
            adB = 0
            hcA[_m][0] = 1
            hcB[_m][0] = 1
            for _n in range(1, self.n):
                if cluster_id in self.copy_numbers:
                    hcA[_m][_n] = self.copy_numbers[cluster_id][0]
                    hcB[_m][_n] = self.copy_numbers[cluster_id][1]
                    continue
                mod = min(self.n, self.k)
                a = min(targetA[_m][_n % mod], self.cn_max)
                b = min(targetB[_m][_n % mod], self.cn_max)

                if self.ampdel:
                    base = self.base
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

                    if a + b > self.cn_max:
                        a = self.cn_max - base
                        b = self.cn_max - a

                hcA[_m][_n] = a
                if a + b <= self.cn_max:
                    hcB[_m][_n] = b
                else:
                    hcB[_m][_n] = max(self.cn_max - a, 0)

                assert hcA[_m][_n] + hcB[_m][_n] <= self.cn_max
                assert hcA[_m][_n] >= 0
                assert hcB[_m][_n] >= 0

        return hcA, hcB

    def hot_start(self, _cA=None, _cB=None):
        def get_rank(hcA, hcB):
            m, n = len(hcA), len(hcA[0])
            rank = np.zeros(n)
            rank[0] = -1
            for _m in range(m):
                for _n in range(1, n):
                    rank[_n] += hcA[_m][_n] * self.symmCoeff(_m) + hcB[_m][
                        _n
                    ] * self.symmCoeff(_m)

            return rank

        if _cA is None and _cB is None:
            _cA, _cB = self.first_hot_start()
        m, n = len(_cA), len(_cA[0])
        assert self.m == m and self.n == n and n >= 1

        rank = get_rank(_cA, _cB)
        rank_indices = np.argsort(rank)

        fixed_rows = getattr(self, "_fixed_rows", set())
        for _m, _n in np.ndindex((m, n)):
            if _m in fixed_rows:
                continue
            self.cA[_m][rank_indices[_n]].value = _cA[_m][_n]
            self.cB[_m][rank_indices[_n]].value = _cB[_m][_n]
        self.warmstart = True

    def get_tree_edges(self):
        """Extract tree edges from DRMST topology variables. Returns dict {child: parent}."""
        if self._var_z is None:
            return None
        tree_edges = {}
        for i in range(1, self.n):
            for j in range(i):
                if (i, j) not in self._var_z:
                    continue
                z_val = self._var_z[(i, j)].value
                if z_val is not None and round(z_val) == 1:
                    tree_edges[i] = j
                    break
            if i not in tree_edges:
                logging.warning(
                    f"DRMST: clone {i} has no parent assigned, defaulting to root"
                )
                tree_edges[i] = 0
        return tree_edges

    def fix_u(self, u):
        self._fixed_u[:] = u
        self.mode = "CARCH"
        # TODO: sanity checks in fixU

    def fix_c(self, cA, cB):
        self._fixed_cA = cA
        self._fixed_cB = cB
        self.mode = "UARCH"
        # TODO: sanity checks in fixC

    def build_random_u(self, random_seed=None, method="dirichlet", dir_alpha=0.3):
        """Generate a random U initialization matrix (n_clones x n_samples).

        Args:
            random_seed: Optional seed for reproducibility.
            method: ``"dirichlet"`` or ``"bubble"``.
            dir_alpha: Dirichlet concentration parameter (lower = sparser).
        """
        with Random(random_seed):
            if method == "dirichlet":
                return self._build_random_u_dirichlet(dir_alpha)
            return self._build_random_u_bubble()

    def _build_random_u_dirichlet(self, alpha=0.3):
        U = np.empty((self.n, self.k))
        n_tumor = self.n - 1
        for _k in range(self.k):
            sid = self.sample_ids[_k]
            if self.purities is not None and sid in self.purities:
                purity = self.purities[sid]
                U[0, _k] = 1 - purity
                if n_tumor == 1:
                    U[1, _k] = purity
                else:
                    t = np.random.dirichlet(alpha * np.ones(n_tumor))
                    t[t < self.minprop] = 0
                    if t.sum() > 0:
                        t = t / t.sum()
                    else:
                        t = np.zeros(n_tumor)
                        t[np.random.randint(n_tumor)] = 1.0
                    U[1:, _k] = purity * t
            else:
                t = np.random.dirichlet(alpha * np.ones(self.n))
                t[t < self.minprop] = 0
                if t.sum() > 0:
                    t = t / t.sum()
                else:
                    t = np.zeros(self.n)
                    t[0] = 1.0
                U[:, _k] = t
        return U

    def _build_random_u_bubble(self):
        def _calculate_size_bubbles(minprop):
            if minprop <= 0.1:
                return 10
            elif minprop <= 0.15:
                return 6
            elif minprop <= 0.2:
                return 5
            else:
                return 3

        def _build_partition_vector(n, n_parts, size_bubbles, minprop=0.03):
            result = np.zeros(n)
            positions = np.random.choice(np.arange(n), n_parts, replace=False)
            bubbles = (
                np.sort(
                    np.random.choice(
                        np.arange(1, size_bubbles), n_parts - 1, replace=False
                    )
                )
                / size_bubbles
            )
            _result = np.diff(bubbles, prepend=0, append=1)
            result[positions] = _result
            result[(minprop - self.tol <= result) & (result < minprop)] = minprop
            return result

        size_bubbles = _calculate_size_bubbles(self.minprop)
        U = np.empty((self.n, self.k))
        for _k in range(self.k):
            sid = self.sample_ids[_k]
            if self.purities is not None and sid in self.purities:
                purity = self.purities[sid]
                U[0, _k] = 1 - purity
                n_tumor = self.n - 1
                if n_tumor == 1:
                    U[1, _k] = purity
                else:
                    _n0 = np.random.randint(1, n_tumor + 1)
                    _n1 = np.random.randint(1, n_tumor + 1)
                    n_parts = min(max(_n0, _n1), size_bubbles)
                    v = _build_partition_vector(
                        n_tumor, n_parts, size_bubbles, minprop=self.minprop
                    )
                    U[1:, _k] = purity * v
            else:
                _n0 = np.random.randint(1, self.n + 1)
                _n1 = np.random.randint(1, self.n + 1)
                n_parts = min(max(_n0, _n1), size_bubbles)
                v = _build_partition_vector(
                    self.n, n_parts, size_bubbles, minprop=self.minprop
                )
                U[:, _k] = v
        return U

    def run(
        self,
        solver_type="gurobi",
        timelimit=None,
        write_path=None,
        solver=None,
        pool_size=1,
        pool_gap=None,
    ):
        if solver is None:
            solver = self._create_solver(solver_type)
            if pool_size > 1 and solver_type in ("gurobipy", "gurobi"):
                solver.options["PoolSolutions"] = pool_size
                solver.options["PoolSearchMode"] = 0
                if pool_gap is not None:
                    solver.options["PoolGap"] = pool_gap

        # Store for pool extraction
        self._pyomo_solver = solver
        self._solver_type = solver_type

        kwargs = self._build_solve_kwargs(solver, self.warmstart, timelimit)
        results = solver.solve(self.model, **kwargs)
        if not self._check_solver_status(results):
            return None

        if write_path is not None:
            self.model.write(write_path)

        return (
            self.model.obj(),
            [[int(getattr(x, "value", x)) for x in row] for row in self.optimized_cA],
            [[int(getattr(x, "value", x)) for x in row] for row in self.optimized_cB],
            [[getattr(x, "value", x) for x in row] for row in self.optimized_u],
        )

    def get_pool_solutions(self, pool_size=10):
        """Extract additional solutions from Gurobi's solution pool.

        Must be called after run() with a Gurobi solver and pool_size > 1.
        Returns a list of (obj, cA, cB, u) tuples for pool solutions beyond
        the optimal (index 0), which is already returned by run().

        Uses the Pyomo solver's internal variable map to reliably resolve
        Pyomo variables to their Gurobi counterparts (avoids name-mismatch
        issues with getVarByName).
        """
        if not hasattr(self, "_pyomo_solver") or self._solver_type not in (
            "gurobi",
            "gurobipy",
        ):
            return []

        try:
            grb_model = self._pyomo_solver._solver_model
            var_map = self._pyomo_solver._pyomo_var_to_solver_var_map
        except AttributeError:
            return []

        n_pool = grb_model.SolCount
        if n_pool <= 1:
            return []

        def _grb_var(pyomo_var):
            """Look up the Gurobi variable for a Pyomo Var via the solver map."""
            return var_map.get(id(pyomo_var))

        solutions = []
        for sol_idx in range(1, min(n_pool, pool_size)):
            grb_model.setParam("SolutionNumber", sol_idx)
            pool_obj = grb_model.PoolObjVal

            cA = [[0] * self.n for _ in range(self.m)]
            cB = [[0] * self.n for _ in range(self.m)]
            u = [[0.0] * self.k for _ in range(self.n)]

            for _m in range(self.m):
                for _n in range(self.n):
                    pyomo_cA = self.optimized_cA[_m][_n]
                    grb_cA = _grb_var(pyomo_cA) if hasattr(pyomo_cA, "value") else None
                    if grb_cA is not None:
                        cA[_m][_n] = int(round(grb_cA.Xn))
                    else:
                        cA[_m][_n] = int(getattr(pyomo_cA, "value", pyomo_cA))

                    pyomo_cB = self.optimized_cB[_m][_n]
                    grb_cB = _grb_var(pyomo_cB) if hasattr(pyomo_cB, "value") else None
                    if grb_cB is not None:
                        cB[_m][_n] = int(round(grb_cB.Xn))
                    else:
                        cB[_m][_n] = int(getattr(pyomo_cB, "value", pyomo_cB))

            for _n in range(self.n):
                for _k in range(self.k):
                    pyomo_u = self.optimized_u[_n][_k]
                    grb_u = _grb_var(pyomo_u) if hasattr(pyomo_u, "value") else None
                    if grb_u is not None:
                        u[_n][_k] = grb_u.Xn
                    else:
                        u[_n][_k] = getattr(pyomo_u, "value", pyomo_u)

            solutions.append((pool_obj, cA, cB, u))

        return solutions
