"""Per-row lexicographic C-step solver for coordinate descent.

Extends ILPSubset with methods that build and solve single-row sub-problems
(one per cluster) using lexicographic optimization:
  Level-1: minimise PI violations (count of (allele, sample) pairs outside PI)
  Level-2: minimise weighted L1 loss + regularisation
"""

import math

import numpy as np
from pyomo import environ as pe

from hatchet.compute_cn.solve.ilp_subset import ILPSubset


class RowILP(ILPSubset):
    """Per-row lexicographic C-step solver for coordinate descent.

    Extends ILPSubset with methods that build and solve single-row
    sub-problems (one per cluster) using lexicographic optimization:
    Level-1 minimizes PI violations, Level-2 minimizes L1 loss.
    """

    def create_row_lexi_model(self, row_idx):
        """Build a per-row Pyomo model with lexicographic objectives.

        Must call ``fix_u()`` before this method. The resulting model is
        stored in ``self.model`` with two objective expressions:
        ``model.obj_z`` (PI violation count) and ``model.obj_d`` (L1 + reg).

        Args:
            row_idx: 0-based index into the cluster rows (range ``self.m``).
        """
        assert self.mode == "CARCH", "create_row_lexi_model requires fix_u() first"
        _m = row_idx
        n, k = self.n, self.k
        cn_max = self.cn_max
        _base = self.base
        _M = self.M
        zero_cn_thres = self.zero_cn_thres
        max_ncns_seg = self.max_ncns_seg

        cluster_id = self.f_a.index[_m]
        cAB_bound = max(sum(self.copy_numbers.get(cluster_id, (0, 0))), cn_max)

        model = pe.ConcreteModel()
        model.constraints = pe.ConstraintList()

        # -- decision variables: cA[n], cB[n] for this row only -----------
        row_cA = [pe.Var(bounds=(0, cAB_bound), domain=pe.Integers) for _ in range(n)]
        row_cB = [pe.Var(bounds=(0, cAB_bound), domain=pe.Integers) for _ in range(n)]
        for _n in range(n):
            model.add_component(f"cA_{_n}", row_cA[_n])
            model.add_component(f"cB_{_n}", row_cB[_n])

        # -- auxiliary variables: fA[k], fB[k], yA[k], yB[k] -------------
        fA, fB, yA, yB = {}, {}, {}, {}
        for _k in range(k):
            yA[_k] = pe.Var(bounds=(0, np.inf), domain=pe.Reals)
            model.add_component(f"yA_{_k}", yA[_k])
            yB[_k] = pe.Var(bounds=(0, np.inf), domain=pe.Reals)
            model.add_component(f"yB_{_k}", yB[_k])
            fA[_k] = pe.Var(bounds=(0, cAB_bound), domain=pe.Reals)
            model.add_component(f"fA_{_k}", fA[_k])
            fB[_k] = pe.Var(bounds=(0, cAB_bound), domain=pe.Reals)
            model.add_component(f"fB_{_k}", fB[_k])

        # -- bit-decomposition variables (max_ncns_seg) -------------------
        bitcA, bitcB, z = {}, {}, {}
        if max_ncns_seg > 0:
            for _b, _n in np.ndindex((_M, n)):
                bitcA[(_b, _n)] = pe.Var(bounds=(0, 1), domain=pe.Binary)
                model.add_component(f"bitcA_{_b}_{_n}", bitcA[(_b, _n)])
                bitcB[(_b, _n)] = pe.Var(bounds=(0, 1), domain=pe.Binary)
                model.add_component(f"bitcB_{_b}_{_n}", bitcB[(_b, _n)])
            for _n in range(1, n):
                for _d in range(max_ncns_seg):
                    z[(_n, _d)] = pe.Var(bounds=(0, 1), domain=pe.Binary)
                    model.add_component(f"z_{_n}_{_d}", z[(_n, _d)])

        # -- PI violation indicators (binary) -----------------------------
        vA, vB = {}, {}
        for _k in range(k):
            vA[_k] = pe.Var(bounds=(0, 1), domain=pe.Binary)
            model.add_component(f"vA_{_k}", vA[_k])
            vB[_k] = pe.Var(bounds=(0, 1), domain=pe.Binary)
            model.add_component(f"vB_{_k}", vB[_k])

        # ---- CONSTRAINTS (via BaseSolver helpers, single-row) ---------------

        # Variable accessors for this single row
        get_cA = lambda m, _n: row_cA[_n]  # noqa: E731
        get_cB = lambda m, _n: row_cB[_n]  # noqa: E731
        row = [_m]

        # L1 linearisation: yA >= |f_a_obs - fA|
        self._add_l1_constraints(
            model,
            lambda m, _k: fA[_k],
            lambda m, _k: fB[_k],
            lambda m, _k: yA[_k],
            lambda m, _k: yB[_k],
            rows=row,
        )

        # Mixture: fA[k] == sum_n cA[n] * u_fixed[n][k]
        for _k in range(k):
            _sumA = 0
            _sumB = 0
            for _n in range(n):
                if self._fixed_u[_n][_k] >= self.minprop - self.tol:
                    _sumA += row_cA[_n] * self._fixed_u[_n][_k]
                    _sumB += row_cB[_n] * self._fixed_u[_n][_k]
            model.constraints.add(fA[_k] == _sumA)
            model.constraints.add(fB[_k] == _sumB)

        self._add_cAB_upper_bound(model, get_cA, get_cB, rows=row)
        self._add_normal_clone_constraints(model, get_cA, get_cB, rows=row)
        self._add_zero_cn_constraints(model, get_cA, get_cB, rows=row)
        self._add_ampdel_constraints(model, get_cA, get_cB, rows=row)
        self._add_fixed_cn_constraints(model, get_cA, get_cB, rows=row)

        # Bit decomposition + max_ncns_seg constraints
        if max_ncns_seg > 0:
            for _n in range(n):
                sum_a = 0
                sum_b = 0
                for _b in range(_M):
                    sum_a += bitcA[(_b, _n)] * math.pow(2, _b)
                    sum_b += bitcB[(_b, _n)] * math.pow(2, _b)
                model.constraints.add(row_cA[_n] == sum_a)
                model.constraints.add(row_cB[_n] == sum_b)

            for _n in range(1, n):
                _sum = 0
                for _d in range(max_ncns_seg):
                    _sum += z[(_n, _d)]
                model.constraints.add(_sum == 1)

            for _b in range(_M):
                for _d in range(max_ncns_seg):
                    for _i in range(1, n - 1):
                        for _j in range(1, n):
                            model.constraints.add(
                                bitcA[(_b, _i)] - bitcA[(_b, _j)]
                                <= 2 - z[(_i, _d)] - z[(_j, _d)]
                            )
                            model.constraints.add(
                                bitcA[(_b, _j)] - bitcA[(_b, _i)]
                                <= 2 - z[(_i, _d)] - z[(_j, _d)]
                            )
                            model.constraints.add(
                                bitcB[(_b, _i)] - bitcB[(_b, _j)]
                                <= 2 - z[(_i, _d)] - z[(_j, _d)]
                            )
                            model.constraints.add(
                                bitcB[(_b, _j)] - bitcB[(_b, _i)]
                                <= 2 - z[(_i, _d)] - z[(_j, _d)]
                            )

            for _d in range(max_ncns_seg - 1):
                _sum_l = _sum_l1 = 0
                for _n in range(1, n):
                    _sum_l += z[(_n, _d)] * self.symmCoeff(_n)
                    _sum_l1 += z[(_n, _d + 1)] * self.symmCoeff(_n)
                    model.constraints.add(_sum_l <= _sum_l1)

        # ---- PI violation indicator constraints (big-M) -----------------
        # vA[k]=0 implies fA[k] in [fa_lo, fa_hi]; vA[k]=1 if violated
        fcn = self.fcn_data
        sample_ids = self.f_a.columns
        bigM = float(cAB_bound + 1)
        for _k in range(k):
            sid = sample_ids[_k]
            fa_lo_val = float(fcn["fa_lo"].loc[cluster_id, sid])
            fa_hi_val = float(fcn["fa_hi"].loc[cluster_id, sid])
            fb_lo_val = float(fcn["fb_lo"].loc[cluster_id, sid])
            fb_hi_val = float(fcn["fb_hi"].loc[cluster_id, sid])
            # fA[k] >= fa_lo - bigM * vA[k]  (if vA=0, fA must be >= fa_lo)
            model.constraints.add(fA[_k] >= fa_lo_val - bigM * vA[_k])
            # fA[k] <= fa_hi + bigM * vA[k]  (if vA=0, fA must be <= fa_hi)
            model.constraints.add(fA[_k] <= fa_hi_val + bigM * vA[_k])
            # same for B
            model.constraints.add(fB[_k] >= fb_lo_val - bigM * vB[_k])
            model.constraints.add(fB[_k] <= fb_hi_val + bigM * vB[_k])

        # ---- OBJECTIVES ------------------------------------------------
        # Level-1: PI violation count (to minimize first)
        obj_z = sum(vA[_k] + vB[_k] for _k in range(k))
        model.obj_z = pe.Expression(expr=obj_z)

        # Level-2: weighted L1 loss + regularisation (to minimize second)
        w_m = self.w[cluster_id]
        obj_d = sum((yA[_k] + yB[_k]) * w_m for _k in range(k))

        pname = self.penalty_param[0]
        init_val = self.penalty_param[1]
        if pname != "RAW":
            pparam = pe.Param(mutable=True, initialize=init_val)
            model.pparam = pparam

            if pname == "MAXCN":
                hcA_var = pe.Var(bounds=(0, np.inf), domain=pe.Reals)
                model.add_component("HCA", hcA_var)
                hcB_var = pe.Var(bounds=(0, np.inf), domain=pe.Reals)
                model.add_component("HCB", hcB_var)
                for _n in range(1, n):
                    model.constraints.add(row_cA[_n] <= hcA_var)
                    model.constraints.add(row_cB[_n] <= hcB_var)
                obj_d += pparam * w_m * hcA_var + pparam * w_m * hcB_var

            elif pname == "DROOT_SUM":
                for _n in range(1, n):
                    mdA = pe.Var(bounds=(0, np.inf), domain=pe.Reals)
                    model.add_component(f"MDA_{_n}", mdA)
                    mdB = pe.Var(bounds=(0, np.inf), domain=pe.Reals)
                    model.add_component(f"MDB_{_n}", mdB)
                    model.constraints.add(row_cA[_n] - row_cA[0] <= mdA)
                    model.constraints.add(row_cA[0] - row_cA[_n] <= mdA)
                    model.constraints.add(row_cB[_n] - row_cB[0] <= mdB)
                    model.constraints.add(row_cB[0] - row_cB[_n] <= mdB)
                    obj_d += pparam * w_m * mdA + pparam * w_m * mdB

            elif pname == "DMRCA_SUM" and n >= 3:
                for _n in range(2, n):
                    mdA = pe.Var(bounds=(0, np.inf), domain=pe.Reals)
                    model.add_component(f"MMDA_{_n}", mdA)
                    mdB = pe.Var(bounds=(0, np.inf), domain=pe.Reals)
                    model.add_component(f"MMDB_{_n}", mdB)
                    model.constraints.add(row_cA[_n] - row_cA[1] <= mdA)
                    model.constraints.add(row_cA[1] - row_cA[_n] <= mdA)
                    model.constraints.add(row_cB[_n] - row_cB[1] <= mdB)
                    model.constraints.add(row_cB[1] - row_cB[_n] <= mdB)
                    obj_d += pparam * w_m * mdA + pparam * w_m * mdB
                for _n in range(2, n):
                    model.constraints.add(row_cA[_n] <= row_cA[1] * cn_max)
                    model.constraints.add(row_cB[_n] <= row_cB[1] * cn_max)

            elif pname == "DADJ_SUM":
                for _n1 in range(n - 1):
                    for _n2 in range(_n1 + 1, n):
                        mdA = pe.Var(bounds=(0, np.inf), domain=pe.Reals)
                        model.add_component(f"MDA_{_n1}_{_n2}", mdA)
                        mdB = pe.Var(bounds=(0, np.inf), domain=pe.Reals)
                        model.add_component(f"MDB_{_n1}_{_n2}", mdB)
                        model.constraints.add(row_cA[_n1] - row_cA[_n2] <= mdA)
                        model.constraints.add(row_cA[_n2] - row_cA[_n1] <= mdA)
                        model.constraints.add(row_cB[_n1] - row_cB[_n2] <= mdB)
                        model.constraints.add(row_cB[_n2] - row_cB[_n1] <= mdB)
                        obj_d += pparam * w_m * mdA + pparam * w_m * mdB

        model.obj_d = pe.Expression(expr=obj_d)
        self.model = model
        self._row_cA = row_cA
        self._row_cB = row_cB
        self._row_idx = _m

    def hot_start_row(self, cA_row, cB_row):
        """Warm-start the per-row model from a previous iteration's CN values.

        Args:
            cA_row: List of length n with integer cA values for this cluster row.
            cB_row: List of length n with integer cB values for this cluster row.
        """
        for _n in range(self.n):
            self._row_cA[_n].value = cA_row[_n]
            self._row_cB[_n].value = cB_row[_n]
        self.warmstart = True

    def run_row_lexi(
        self, solver_type="gurobi", timelimit=None, solver=None, pparam=0.0
    ):
        """Solve the per-row lexicographic sub-problem.

        Level-1: minimise PI violation count (obj_z).
        Level-2: minimise L1 loss + regularisation (obj_d), without worsening Level-1.

        For Gurobi: uses native multi-objective API (setObjectiveN).
        For CBC/other: two sequential solves.

        Returns:
            ``(obj_d_value, cA_row, cB_row)`` on success, or ``None``.
        """
        if solver is None:
            if solver_type in ("gurobipy", "gurobi"):
                solver = pe.SolverFactory("gurobi", solver_io="python")
                solver.options["OutputFlag"] = 0
                solver.options["LogToConsole"] = 0
                solver.options["LogFile"] = ""
            else:
                solver = pe.SolverFactory(solver_type)

        # Set regularization parameter
        if hasattr(self.model, "pparam"):
            self.model.pparam = pparam

        if solver_type in ("gurobipy", "gurobi"):
            return self._run_row_lexi_gurobi(solver, timelimit)
        else:
            return self._run_row_lexi_sequential(solver, solver_type, timelimit)

    def _run_row_lexi_gurobi(self, solver, timelimit):
        """Gurobi path: native lexicographic multi-objective via setObjectiveN.

        Uses Gurobi's built-in hierarchical optimization — a single solve
        call that first optimises Level-1 then Level-2.
        """
        import gurobipy as gp

        model = self.model
        # Pyomo needs an objective to send the model to Gurobi; we use a
        # dummy that gets replaced by setObjectiveN below.
        model.obj_dummy = pe.Objective(expr=0, sense=pe.minimize)

        kwargs = {"report_timing": False}
        if timelimit is not None:
            kwargs["timelimit"] = int(timelimit)
        if solver.warm_start_capable():
            kwargs["warmstart"] = self.warmstart

        # Send model to Gurobi (creates the internal gurobipy model)
        solver.solve(model, load_solutions=False, **kwargs)
        grb_model = solver._solver_model

        # Map Pyomo expressions to Gurobi linear expressions
        def _pyomo_to_grb(pyomo_expr):
            """Convert a Pyomo linear expression to a gurobipy LinExpr."""
            from pyomo.repn import generate_standard_repn

            var_map = solver._pyomo_var_to_solver_var_map
            grb_expr = gp.LinExpr()
            repn = generate_standard_repn(pyomo_expr)
            if repn.constant:
                grb_expr.addConstant(repn.constant)
            for coef, var in zip(repn.linear_coefs, repn.linear_vars):
                grb_expr.add(var_map[var], coef)
            return grb_expr

        grb_z = _pyomo_to_grb(model.obj_z.expr)
        grb_d = _pyomo_to_grb(model.obj_d.expr)

        # Set lexicographic objectives: Level-1 (priority=1), Level-2 (priority=0)
        grb_model.ModelSense = gp.GRB.MINIMIZE
        grb_model.setObjectiveN(grb_z, index=0, priority=1, name="PI_violations")
        grb_model.setObjectiveN(grb_d, index=1, priority=0, name="L1_loss")
        grb_model.optimize()

        model.del_component(model.obj_dummy)

        if grb_model.Status not in (
            gp.GRB.OPTIMAL,
            gp.GRB.SUBOPTIMAL,
            gp.GRB.TIME_LIMIT,
        ):
            return None

        # Read solution back from Gurobi
        var_map = solver._pyomo_var_to_solver_var_map
        cA_row = [int(round(var_map[self._row_cA[_n]].X)) for _n in range(self.n)]
        cB_row = [int(round(var_map[self._row_cB[_n]].X)) for _n in range(self.n)]
        obj_d_val = grb_model.getObjective(1).getValue()
        return obj_d_val, cA_row, cB_row

    def _run_row_lexi_sequential(self, solver, solver_type, timelimit):
        """CBC/generic path: two sequential Pyomo solves.

        Phase 1: minimize PI violations (obj_z).
        Phase 2: fix Z = Z*, minimize L1 loss + reg (obj_d).
        """
        model = self.model

        kwargs = {"report_timing": False}
        if timelimit is not None:
            kwargs["timelimit"] = int(timelimit)
        if solver.warm_start_capable():
            kwargs["warmstart"] = self.warmstart

        # Phase 1: minimize PI violations
        model.obj_phase1 = pe.Objective(expr=model.obj_z.expr, sense=pe.minimize)
        results = solver.solve(model, **kwargs)
        if not self._check_solver_status(results):
            model.del_component(model.obj_phase1)
            return None

        z_star = int(round(pe.value(model.obj_phase1)))
        model.obj_phase1.deactivate()

        # Phase 2: fix Z = Z*, minimize L1 loss + reg
        model.z_fix = pe.Constraint(expr=model.obj_z.expr == z_star)
        model.obj_phase2 = pe.Objective(expr=model.obj_d.expr, sense=pe.minimize)
        results = solver.solve(model, **kwargs)

        result = None
        if self._check_solver_status(results):
            cA_row = [int(self._row_cA[_n].value) for _n in range(self.n)]
            cB_row = [int(self._row_cB[_n].value) for _n in range(self.n)]
            result = pe.value(model.obj_phase2), cA_row, cB_row

        # Cleanup for next row
        model.del_component(model.obj_phase1)
        model.del_component(model.obj_phase2)
        model.del_component(model.z_fix)
        return result

