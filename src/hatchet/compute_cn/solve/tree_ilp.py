"""Monolithic tree-constrained C-step MILP for coordinate descent.

Builds a single Pyomo model with cA/cB variables for ALL clusters,
global arborescence topology variables z[i][j], nearest-neighbor
distance variables M[c][i], and a penalized objective:

    min weighted_IMF_loss + pparam * weighted_tree_edge_length

The model is built **once** via ``build_model()`` and reused across CD
iterations by calling ``update_u()``, ``update_lambda()``, ``update_degree()``.
"""

import numpy as np
from pyomo import environ as pe

from hatchet.compute_cn.solve.base_solver import BaseSolver


class TreeILP(BaseSolver):
    """Monolithic tree-constrained C-step MILP.

    All m clusters are solved jointly in a single model, linked by global
    binary tree topology variables z[i][j] that define an r-arborescence
    rooted at clone 0 (normal).
    """

    def __init__(self, n, cn_max, fcn_data, w, copy_numbers,
                 ampdel=True, base=1, minprop=0.01,
                 zero_cn_thres=0.005, tol=0.001):
        super().__init__(
            n=n, cn_max=cn_max, fcn_data=fcn_data, w=w,
            copy_numbers=copy_numbers, ampdel=ampdel, base=base,
            minprop=minprop, zero_cn_thres=zero_cn_thres, tol=tol,
        )
        self._var_cA = None
        self._var_cB = None
        self._var_z = None
        self._var_M = None
        self._param_u = None
        self._param_lam = None
        self._param_deg = None

    def build_model(self, fixed_u, max_degree=None):
        """Build the monolithic tree MILP.

        Call once per worker. Use ``update_u()``, ``update_lambda()``,
        ``update_degree()`` to change parameters between solves.
        """
        m, n, k = self.m, self.n, self.k
        cn_max = self.cn_max
        big_M = 2 * cn_max

        model = pe.ConcreteModel()
        model.constraints = pe.ConstraintList()

        # ---- Mutable Params ---------------------------------------------------
        param_u = {}
        for _n in range(n):
            for _k in range(k):
                p = pe.Param(mutable=True, initialize=float(fixed_u[_n][_k]))
                model.add_component(f"pu_{_n}_{_k}", p)
                param_u[(_n, _k)] = p

        param_lam = pe.Param(mutable=True, initialize=0.0)
        model.add_component("p_lambda", param_lam)

        param_deg = pe.Param(
            mutable=True,
            initialize=float(max_degree if max_degree is not None else n),
        )
        model.add_component("p_max_degree", param_deg)

        # ---- Decision variables: cA[m][n], cB[m][n] -------------------------
        var_cA = {}
        var_cB = {}
        for _m in range(m):
            bound = self.cAB_bound(self.cluster_ids[_m])
            for _n in range(n):
                var_cA[(_m, _n)] = pe.Var(bounds=(0, bound), domain=pe.Integers)
                model.add_component(f"cA_{_m}_{_n}", var_cA[(_m, _n)])
                var_cB[(_m, _n)] = pe.Var(bounds=(0, bound), domain=pe.Integers)
                model.add_component(f"cB_{_m}_{_n}", var_cB[(_m, _n)])

        get_cA = lambda _m, _n: var_cA[(_m, _n)]
        get_cB = lambda _m, _n: var_cB[(_m, _n)]

        # ---- Shared constraints from BaseSolver ------------------------------
        self._add_cAB_upper_bound(model, get_cA, get_cB)
        self._add_normal_clone_constraints(model, get_cA, get_cB)
        self._add_zero_cn_constraints(model, get_cA, get_cB)
        self._add_ampdel_constraints(model, get_cA, get_cB)
        self._add_fixed_cn_constraints(model, get_cA, get_cB)

        # ---- Global topology variables: z[i][j] (i >= 1, j < i) -------------
        var_z = {}
        for i in range(1, n):
            for j in range(i):
                var_z[(i, j)] = pe.Var(bounds=(0, 1), domain=pe.Binary)
                model.add_component(f"z_{i}_{j}", var_z[(i, j)])

        for i in range(1, n):
            model.constraints.add(
                sum(var_z[(i, j)] for j in range(i)) == 1
            )

        # ---- Degree constraint on all nodes (mutable) ------------------------
        out_deg_0 = sum(var_z[(i, 0)] for i in range(1, n))
        model.add_component(
            "con_degree_0", pe.Constraint(expr=out_deg_0 <= param_deg)
        )
        for c in range(1, n):
            out_deg_c = sum(var_z[(i, c)] for i in range(c + 1, n))
            model.add_component(
                f"con_degree_{c}",
                pe.Constraint(expr=1 + out_deg_c <= param_deg),
            )

        # ---- Nearest-neighbor distance variables -----------------------------
        var_M = {}
        for _m in range(m):
            for i in range(1, n):
                var_M[(_m, i)] = pe.Var(bounds=(0, None), domain=pe.Reals)
                model.add_component(f"M_{_m}_{i}", var_M[(_m, i)])
                for j in range(i):
                    dA = pe.Var(bounds=(0, None), domain=pe.Reals)
                    model.add_component(f"dA_{_m}_{i}_{j}", dA)
                    dB = pe.Var(bounds=(0, None), domain=pe.Reals)
                    model.add_component(f"dB_{_m}_{i}_{j}", dB)

                    model.constraints.add(
                        dA >= var_cA[(_m, i)] - var_cA[(_m, j)]
                        - big_M * (1 - var_z[(i, j)])
                    )
                    model.constraints.add(
                        dA >= var_cA[(_m, j)] - var_cA[(_m, i)]
                        - big_M * (1 - var_z[(i, j)])
                    )
                    model.constraints.add(
                        dB >= var_cB[(_m, i)] - var_cB[(_m, j)]
                        - big_M * (1 - var_z[(i, j)])
                    )
                    model.constraints.add(
                        dB >= var_cB[(_m, j)] - var_cB[(_m, i)]
                        - big_M * (1 - var_z[(i, j)])
                    )
                    model.constraints.add(
                        var_M[(_m, i)] >= dA + dB
                        - big_M * (1 - var_z[(i, j)])
                    )

        # ---- Dynamic LOH constraints ----------------------------------------
        for _m in range(m):
            for _n in range(n):
                lA = pe.Var(bounds=(0, 1), domain=pe.Binary)
                model.add_component(f"lostA_{_m}_{_n}", lA)
                lB = pe.Var(bounds=(0, 1), domain=pe.Binary)
                model.add_component(f"lostB_{_m}_{_n}", lB)
                model.constraints.add(var_cA[(_m, _n)] >= 1 - big_M * lA)
                model.constraints.add(var_cA[(_m, _n)] <= big_M * (1 - lA))
                model.constraints.add(var_cB[(_m, _n)] >= 1 - big_M * lB)
                model.constraints.add(var_cB[(_m, _n)] <= big_M * (1 - lB))

        for _m in range(m):
            for i in range(1, n):
                for j in range(i):
                    lA_j = model.find_component(f"lostA_{_m}_{j}")
                    lB_j = model.find_component(f"lostB_{_m}_{j}")
                    model.constraints.add(
                        var_cA[(_m, i)]
                        <= cn_max * (1 - lA_j)
                        + big_M * (1 - var_z[(i, j)])
                    )
                    model.constraints.add(
                        var_cB[(_m, i)]
                        <= cn_max * (1 - lB_j)
                        + big_M * (1 - var_z[(i, j)])
                    )

        # ---- Predicted fractional CN (mixture via mutable Param u) -----------
        var_fA = {}
        var_fB = {}
        var_yA = {}
        var_yB = {}
        for _m in range(m):
            bound = self.cAB_bound(self.cluster_ids[_m])
            for _k in range(k):
                var_fA[(_m, _k)] = pe.Var(bounds=(0, bound), domain=pe.Reals)
                model.add_component(f"fA_{_m}_{_k}", var_fA[(_m, _k)])
                var_fB[(_m, _k)] = pe.Var(bounds=(0, bound), domain=pe.Reals)
                model.add_component(f"fB_{_m}_{_k}", var_fB[(_m, _k)])
                var_yA[(_m, _k)] = pe.Var(bounds=(0, None), domain=pe.Reals)
                model.add_component(f"yA_{_m}_{_k}", var_yA[(_m, _k)])
                var_yB[(_m, _k)] = pe.Var(bounds=(0, None), domain=pe.Reals)
                model.add_component(f"yB_{_m}_{_k}", var_yB[(_m, _k)])

        # Mixture constraints
        for _m in range(m):
            for _k in range(k):
                sumA = sum(
                    var_cA[(_m, _n)] * param_u[(_n, _k)] for _n in range(n)
                )
                sumB = sum(
                    var_cB[(_m, _n)] * param_u[(_n, _k)] for _n in range(n)
                )
                model.add_component(
                    f"mix_A_{_m}_{_k}",
                    pe.Constraint(expr=var_fA[(_m, _k)] == sumA),
                )
                model.add_component(
                    f"mix_B_{_m}_{_k}",
                    pe.Constraint(expr=var_fB[(_m, _k)] == sumB),
                )

        # L1 linearisation
        self._add_l1_constraints(
            model,
            lambda _m, _k: var_fA[(_m, _k)],
            lambda _m, _k: var_fB[(_m, _k)],
            lambda _m, _k: var_yA[(_m, _k)],
            lambda _m, _k: var_yB[(_m, _k)],
        )

        # ---- Objective expressions -------------------------------------------
        obj_imf = sum(
            self.w[self.cluster_ids[_m]] * (var_yA[(_m, _k)] + var_yB[(_m, _k)])
            for _m in range(m) for _k in range(k)
        )
        model.obj_imf = pe.Expression(expr=obj_imf)

        obj_tree = sum(
            self.w[self.cluster_ids[_m]] * var_M[(_m, i)]
            for _m in range(m) for i in range(1, n)
        )
        model.obj_tree = pe.Expression(expr=obj_tree)

        # Objective: IMF + pparam * tree_edge_length
        model.obj = pe.Objective(
            expr=obj_imf + param_lam * obj_tree,
            sense=pe.minimize,
        )

        # ---- Store references ------------------------------------------------
        self.model = model
        self._var_cA = var_cA
        self._var_cB = var_cB
        self._var_z = var_z
        self._var_M = var_M
        self._param_u = param_u
        self._param_lam = param_lam
        self._param_deg = param_deg
        self.warmstart = False

    # ------------------------------------------------------------------
    # Parameter updates (no model rebuild)
    # ------------------------------------------------------------------

    def update_u(self, new_u):
        """Update the fixed u coefficients in the mixture constraints."""
        for _n in range(self.n):
            for _k in range(self.k):
                self._param_u[(_n, _k)].value = float(new_u[_n][_k])

    def update_lambda(self, lam):
        """Update the scalarization weight λ."""
        self._param_lam.value = float(lam)

    def update_degree(self, max_degree):
        """Update the max degree constraint on all nodes."""
        self._param_deg.value = float(max_degree)

    # ------------------------------------------------------------------

    def hot_start(self, cA, cB):
        """Warm-start cA/cB variables from previous iteration."""
        if self._var_cA is None:
            return
        for _m in range(self.m):
            for _n in range(self.n):
                self._var_cA[(_m, _n)].value = int(cA[_m][_n])
                self._var_cB[(_m, _n)].value = int(cB[_m][_n])
        self.warmstart = True

    def run(self, solver_type="gurobi", solver=None, timelimit=None):
        """Solve the tree MILP.

        Returns:
            ``(imf_obj, tree_edge_length, cA, cB, tree_edges)`` or ``None``.
        """
        if solver is None:
            solver = self._create_solver(solver_type)
        kwargs = self._build_solve_kwargs(solver, self.warmstart, timelimit)

        results = solver.solve(self.model, **kwargs)
        if not self._check_solver_status(results):
            return None

        cA = [
            [int(round(self._var_cA[(_m, _n)].value)) for _n in range(self.n)]
            for _m in range(self.m)
        ]
        cB = [
            [int(round(self._var_cB[(_m, _n)].value)) for _n in range(self.n)]
            for _m in range(self.m)
        ]

        tree_edges = {}
        for i in range(1, self.n):
            for j in range(i):
                z_val = self._var_z[(i, j)].value
                if z_val is not None and round(z_val) == 1:
                    tree_edges[i] = j
                    break

        imf_obj = pe.value(self.model.obj_imf)
        tree_edge_length = pe.value(self.model.obj_tree)

        return imf_obj, tree_edge_length, cA, cB, tree_edges

    @staticmethod
    def compute_edge_length_from_tree(cA, cB, tree_edges):
        """Compute total L1 edge length from CN arrays and tree topology."""
        cA = np.asarray(cA)
        cB = np.asarray(cB)
        total = 0.0
        for child, parent in tree_edges.items():
            total += np.sum(np.abs(cA[:, child] - cA[:, parent]))
            total += np.sum(np.abs(cB[:, child] - cB[:, parent]))
        return float(total)
