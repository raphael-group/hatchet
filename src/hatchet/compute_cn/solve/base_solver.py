"""Base class for Pyomo-based copy-number solvers.

Provides shared constructor fields, constraint helper methods, and solver
invocation utilities used by ILPSubset.
"""

import numpy as np
from pyomo import environ as pe
from pyomo.opt import SolverStatus, TerminationCondition


class BaseSolver:
    """Abstract base for copy-number deconvolution solvers.

    Stores shared data (fractional CN, cluster metadata, parameters) and
    provides reusable methods for building common Pyomo constraints and
    managing solver invocation.

    Subclass: ``ILPSubset``.
    """

    def __init__(
        self,
        n,
        cn_max,
        fcn_data,
        w,
        copy_numbers,
        ampdel=True,
        base=1,
        minprop=0.01,
        zero_cn_thres=0.005,
        tol=0.001,
        balanced_clusters=None,
        mrca=False,
    ):
        f_a = fcn_data["fa"]
        f_b = fcn_data["fb"]
        assert f_a.shape == f_b.shape

        self.n = n
        self.cn_max = cn_max
        self.fcn_data = fcn_data
        self.f_a = f_a
        self.f_b = f_b
        self.m, self.k = f_a.shape
        self.cluster_ids = f_a.index
        self.sample_ids = f_a.columns
        self.w = w
        self.copy_numbers = copy_numbers
        self.ampdel = ampdel
        self.balanced_clusters = balanced_clusters or []
        self.mrca = mrca
        self._base = base
        self.zero_cn_thres = zero_cn_thres
        self.tol = tol
        self.minprop = minprop

        self.model = None
        self.warmstart = False

    @property
    def base(self):
        return self._base

    @staticmethod
    def symmCoeff(i):
        return i + 1

    # ------------------------------------------------------------------
    # Helper: per-cluster cAB bound
    # ------------------------------------------------------------------

    def cAB_bound(self, cluster_id):
        """Upper bound on cA + cB for a given cluster."""
        return max(sum(self.copy_numbers.get(cluster_id, (0, 0))), self.cn_max)

    # ------------------------------------------------------------------
    # Shared constraint builders
    # ------------------------------------------------------------------
    #
    # All constraint helpers accept callables ``get_cA(m, n)`` and
    # ``get_cB(m, n)`` that return the Pyomo variable for cluster *m*,
    # clone *n*.  This abstracts over the different variable layouts:
    #
    #   ILPSubset:  lambda m, n: self.cA[m][n]

    def _add_normal_clone_constraints(self, model, get_cA, get_cB, rows=None):
        """Fix clone 0 at (1,1) for every cluster in *rows*."""
        for _m in rows if rows is not None else range(self.m):
            model.constraints.add(get_cA(_m, 0) == 1)
            model.constraints.add(get_cB(_m, 0) == 1)

    def _add_zero_cn_constraints(self, model, get_cA, get_cB, rows=None):
        """Forbid zero total CN for significant clusters, tumor clones."""
        w_total = sum(self.w)
        for _m in rows if rows is not None else range(self.m):
            cluster_id = self.cluster_ids[_m]
            if self.w[cluster_id] / w_total >= self.zero_cn_thres:
                for _n in range(1, self.n):
                    model.constraints.add(get_cA(_m, _n) + get_cB(_m, _n) >= 1)

    def _add_ampdel_constraints(self, model, get_cA, get_cB, rows=None):
        """Add amp/del binary constraints for non-fixed clusters.

        Returns:
            Dict ``{m: (adA_var, adB_var)}`` for clusters that got constraints.
        """
        if not self.ampdel:
            return {}
        cn_max = self.cn_max
        _base = self._base
        ad_vars = {}
        for _m in rows if rows is not None else range(self.m):
            cluster_id = self.cluster_ids[_m]
            if cluster_id in self.copy_numbers:
                continue
            adA = pe.Var(bounds=(0, 1), domain=pe.Binary)
            model.add_component(f"adA_{_m}", adA)
            adB = pe.Var(bounds=(0, 1), domain=pe.Binary)
            model.add_component(f"adB_{_m}", adB)
            for _n in range(1, self.n):
                model.constraints.add(
                    get_cA(_m, _n) <= cn_max * adA + _base - _base * adA
                )
                model.constraints.add(get_cA(_m, _n) >= _base * adA)
                model.constraints.add(
                    get_cB(_m, _n) <= cn_max * adB + _base - _base * adB
                )
                model.constraints.add(get_cB(_m, _n) >= _base * adB)
            ad_vars[_m] = (adA, adB)
        return ad_vars

    def _add_cAB_upper_bound(self, model, get_cA, get_cB, rows=None):
        """cA[m][n] + cB[m][n] <= cAB_bound for all clusters and clones."""
        for _m in rows if rows is not None else range(self.m):
            cluster_id = self.cluster_ids[_m]
            bound = self.cAB_bound(cluster_id)
            for _n in range(self.n):
                model.constraints.add(get_cA(_m, _n) + get_cB(_m, _n) <= bound)

    def _add_fixed_cn_constraints(self, model, get_cA, get_cB, rows=None):
        """Fix CN state for clonal clusters from ``self.copy_numbers``."""
        for _m in rows if rows is not None else range(self.m):
            cluster_id = self.cluster_ids[_m]
            if cluster_id in self.copy_numbers:
                _cnA, _cnB = self.copy_numbers[cluster_id]
                for _n in range(1, self.n):
                    model.constraints.add(get_cA(_m, _n) == _cnA)
                    model.constraints.add(get_cB(_m, _n) == _cnB)

    def _add_balanced_constraints(self, model, get_cA, get_cB, rows=None):
        """Force cA == cB for balanced clusters, tumor clones only."""
        if not self.balanced_clusters:
            return
        bal_set = set(self.balanced_clusters)
        for _m in rows if rows is not None else range(self.m):
            cluster_id = self.cluster_ids[_m]
            if cluster_id in bal_set:
                for _n in range(1, self.n):
                    model.constraints.add(get_cA(_m, _n) == get_cB(_m, _n))

    def _add_mrca_loh_constraints(self, model, get_cA, get_cB, rows=None):
        """If clone 1 (MRCA) has lost an allele, subclones (n>=2) cannot regain it.

        Linearised as: cA[c][n] <= cA[c][1] * cn_max (big-M).
        Only active when ``self.mrca`` is True and n >= 3.
        """
        if not self.mrca or self.n < 3:
            return
        cn_max = self.cn_max
        for _m in rows if rows is not None else range(self.m):
            for _n in range(2, self.n):
                model.constraints.add(get_cA(_m, _n) <= get_cA(_m, 1) * cn_max)
                model.constraints.add(get_cB(_m, _n) <= get_cB(_m, 1) * cn_max)

    def _add_l1_constraints(self, model, get_fA, get_fB, get_yA, get_yB, rows=None):
        """L1 linearisation: yA >= |f_a_obs - fA| for all (m, k).

        ``get_fA/fB/yA/yB`` are callables ``(m, k) -> Var``.
        """
        f_a_vals = self.f_a.values
        f_b_vals = self.f_b.values
        for _m in rows if rows is not None else range(self.m):
            for _k in range(self.k):
                fa_obs = float(f_a_vals[_m, _k])
                fb_obs = float(f_b_vals[_m, _k])
                model.constraints.add(fa_obs - get_fA(_m, _k) <= get_yA(_m, _k))
                model.constraints.add(get_fA(_m, _k) - fa_obs <= get_yA(_m, _k))
                model.constraints.add(fb_obs - get_fB(_m, _k) <= get_yB(_m, _k))
                model.constraints.add(get_fB(_m, _k) - fb_obs <= get_yB(_m, _k))

    def _add_symmetry_breaking(self, model, get_cA, get_cB):
        """Clone ordering via symmCoeff-weighted CN sums."""
        for i in range(1, self.n - 1):
            sum1 = 0
            sum2 = 0
            for _m in range(self.m):
                sc = self.symmCoeff(_m)
                sum1 += get_cA(_m, i) * sc + get_cB(_m, i) * sc
                sum2 += get_cA(_m, i + 1) * sc + get_cB(_m, i + 1) * sc
            model.constraints.add(sum1 <= sum2)

    # ------------------------------------------------------------------
    # Solver invocation utilities
    # ------------------------------------------------------------------

    @staticmethod
    def _create_solver(solver_type, threads=None):
        """Create a Pyomo solver with suppressed output.

        Args:
            solver_type: ``"gurobi"`` or ``"cbc"``.
            threads: Max threads per solve (Gurobi only). ``None`` = solver default.
                Set to 1 when running many parallel workers to avoid contention.
        """
        if solver_type in ("gurobipy", "gurobi"):
            solver = pe.SolverFactory("gurobi", solver_io="python")
            solver.options["OutputFlag"] = 0
            solver.options["LogToConsole"] = 0
            solver.options["LogFile"] = ""
            if threads is not None:
                solver.options["Threads"] = threads
        else:
            solver = pe.SolverFactory(solver_type)
        return solver

    @staticmethod
    def _build_solve_kwargs(solver, warmstart, timelimit):
        """Build kwargs dict for ``solver.solve()``."""
        kwargs = {"report_timing": False}
        if timelimit is not None:
            kwargs["timelimit"] = int(timelimit)
        if solver.warm_start_capable():
            kwargs["warmstart"] = warmstart
        return kwargs

    @staticmethod
    def _check_solver_status(results):
        """Return True if the solver found a usable solution."""
        solver_ok = (
            results.solver.status == SolverStatus.ok
            and results.solver.termination_condition
            in (TerminationCondition.optimal, TerminationCondition.feasible)
        )
        time_limit_hit = (
            results.solver.status == SolverStatus.aborted
            and results.solver.termination_condition
            == TerminationCondition.maxTimeLimit
        )
        return solver_ok or time_limit_hit
