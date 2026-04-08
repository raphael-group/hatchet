"""Solver execution: FullILP (single-shot) and CDSolver (coordinate descent)."""

from __future__ import annotations

import logging
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd
from pyomo import environ as pe

from hatchet.compute_cn.solve.base_solver import BaseSolver
from hatchet.compute_cn.solve.model import (
    build_model,
    first_hot_start,
    hot_start,
    build_random_u,
)
from hatchet.compute_cn.solve.variables import SolverParams
from hatchet.compute_cn.solve.utils import Random


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------


def _solve_model(model, solver, warmstart, timelimit):
    """Solve a Pyomo model. Returns True on success."""
    kwargs = BaseSolver._build_solve_kwargs(solver, warmstart, timelimit)
    results = solver.solve(model, **kwargs)
    return BaseSolver._check_solver_status(results)


def extract_solution(model, p: SolverParams):
    """Extract (obj, cA, cB, u) from solved model."""
    cA = [
        [int(round(model.cA[_m, _n].value)) for _n in range(p.n)] for _m in range(p.m)
    ]
    cB = [
        [int(round(model.cB[_m, _n].value)) for _n in range(p.n)] for _m in range(p.m)
    ]
    u = [[model.u[_n, _k].value for _k in range(p.k)] for _n in range(p.n)]
    return model.obj(), cA, cB, u


def extract_tree_edges(var_z, n):
    """Extract tree edges from DRMST topology variables."""
    if var_z is None:
        return None
    tree_edges = {}
    for i in range(1, n):
        for j in range(i):
            if (i, j) not in var_z:
                continue
            z_val = var_z[(i, j)].value
            if z_val is not None and round(z_val) == 1:
                tree_edges[i] = j
                break
        if i not in tree_edges:
            logging.warning(f"DRMST: clone {i} has no parent, defaulting to root")
            tree_edges[i] = 0
    return tree_edges


# ---------------------------------------------------------------------------
# FullILP
# ---------------------------------------------------------------------------


class FullILP:
    """Single-shot ILP solve."""

    def __init__(self, params: SolverParams, penalty_param):
        self.p = params
        self.model, self.var_z = build_model(params, penalty_param)

    def solve(self, solver_type="gurobi", timelimit=None, pool_size=1, pool_gap=None):
        solver = BaseSolver._create_solver(solver_type)
        if pool_size > 1 and solver_type in ("gurobi", "gurobipy"):
            solver.options["PoolSolutions"] = pool_size
            solver.options["PoolSearchMode"] = 0
            if pool_gap is not None:
                solver.options["PoolGap"] = pool_gap
        self._solver = solver
        self._solver_type = solver_type

        if not _solve_model(self.model, solver, True, timelimit):
            return None
        return extract_solution(self.model, self.p)

    def get_tree_edges(self):
        return extract_tree_edges(self.var_z, self.p.n)

    def get_pool_solutions(self, pool_size=10):
        if not hasattr(self, "_solver") or self._solver_type not in (
            "gurobi",
            "gurobipy",
        ):
            return []
        try:
            grb_model = self._solver._solver_model
            var_map = self._solver._pyomo_var_to_solver_var_map
        except AttributeError:
            return []
        if grb_model.SolCount <= 1:
            return []

        solutions = []
        for sol_idx in range(1, min(grb_model.SolCount, pool_size)):
            grb_model.setParam("SolutionNumber", sol_idx)
            cA = [[0] * self.p.n for _ in range(self.p.m)]
            cB = [[0] * self.p.n for _ in range(self.p.m)]
            u = [[0.0] * self.p.k for _ in range(self.p.n)]
            for _m in range(self.p.m):
                for _n in range(self.p.n):
                    for arr, cX in [(cA, self.model.cA), (cB, self.model.cB)]:
                        pv = cX[_m, _n]
                        gv = var_map.get(id(pv))
                        arr[_m][_n] = int(round(gv.Xn)) if gv else int(round(pv.value))
            for _n in range(self.p.n):
                for _k in range(self.p.k):
                    pv = self.model.u[_n, _k]
                    gv = var_map.get(id(pv))
                    u[_n][_k] = gv.Xn if gv else pv.value
            solutions.append((grb_model.PoolObjVal, cA, cB, u))
        return solutions


# ---------------------------------------------------------------------------
# Coordinate Descent
# ---------------------------------------------------------------------------

_cd_global = None
_cd_solver_cache = None


def _init_cd_worker(cd, log_level):
    global _cd_global, _cd_solver_cache
    _cd_global = cd
    _cd_solver_cache = None
    logging.basicConfig(
        level=logging.ERROR,
        format="%(asctime)s.%(msecs)03d %(levelname)s [worker] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        force=True,
    )
    for name in ("pyomo", "pyomo.core"):
        logging.getLogger(name).setLevel(logging.ERROR)


def _cd_work(
    work_id, u, pparam, solver_type, max_iters, max_convergence_iters, timelimit
):
    global _cd_solver_cache
    cd = _cd_global
    p = cd.params

    if _cd_solver_cache is None:
        _cd_solver_cache = BaseSolver._create_solver(
            solver_type, threads=cd.solver_threads
        )
    solver = _cd_solver_cache

    _u, _cA, _cB = u, cd.hcA, cd.hcB
    _prev_obj_u = None
    _imf_c = _reg_c = 0.0
    _tree_edges = None
    _iters = _conv_iters = 0

    while _iters < max_iters and _conv_iters < max_convergence_iters:
        # C-step: fix u, optimize cA/cB
        p_c = SolverParams(**{**p.__dict__, "mode": "CARCH"})
        model_c, var_z_c = build_model(p_c, cd.penalty_param, fixed_u=_u)
        hot_start(model_c, p_c, _cA, _cB)
        model_c.pparam = pparam
        if not _solve_model(model_c, solver, True, timelimit):
            return None
        _, _cA, _cB, _ = extract_solution(model_c, p_c)
        _imf_c = pe.value(model_c.obj_imf)
        _reg_c = pe.value(model_c.obj_reg)
        _tree_edges = extract_tree_edges(var_z_c, p.n)

        # U-step: fix cA/cB, optimize u
        p_u = SolverParams(**{**p.__dict__, "mode": "UARCH"})
        model_u, _ = build_model(p_u, cd.penalty_param, fixed_cA=_cA, fixed_cB=_cB)
        if not _solve_model(model_u, solver, False, timelimit):
            return None
        _obj_u, _, _, _u = extract_solution(model_u, p_u)

        if _prev_obj_u is not None:
            if abs(_obj_u - _prev_obj_u) < cd.cd_tol:
                _conv_iters += 1
            else:
                _conv_iters = 0
        _prev_obj_u = _obj_u
        _iters += 1

    return _obj_u, _cA, _cB, _u, _imf_c, _reg_c, _tree_edges


class CDSolver:
    """Coordinate-descent solver with parallel restarts over a regularization path."""

    def __init__(
        self,
        params: SolverParams,
        reg_term="RAW",
        reg_steps=0,
        reg_bound=0.3,
        u_init_method="dirichlet",
        u_dir_alpha=0.3,
        solver_threads=None,
        cd_tol=0.001,
    ):
        self.params = params
        self.reg_name = reg_term if reg_term else "RAW"
        self.reg_steps = reg_steps
        self.reg_bound = reg_bound
        self.u_init_method = u_init_method
        self.u_dir_alpha = u_dir_alpha
        self.solver_threads = solver_threads
        self.cd_tol = cd_tol
        self.penalty_param = [self.reg_name, 0.0]
        self.hcA, self.hcB = first_hot_start(params)

    def run(
        self,
        solver_type="gurobi",
        max_iters=10,
        max_convergence_iters=2,
        n_seed=400,
        j=8,
        random_seed=None,
        timelimit=None,
        u0_tsv_path=None,
        **_,
    ):
        with Random(random_seed):
            seeds = [
                build_random_u(
                    self.params, method=self.u_init_method, alpha=self.u_dir_alpha
                )
                for _ in range(n_seed)
            ]

        if u0_tsv_path is not None:
            rows = []
            for restart, u in enumerate(seeds):
                for clone in range(u.shape[0]):
                    row = {"restart": restart, "clone": clone}
                    for j_idx, sid in enumerate(self.params.sample_ids):
                        row[sid] = u[clone, j_idx]
                    rows.append(row)
            pd.DataFrame(rows).to_csv(u0_tsv_path, sep="\t", index=False)

        no_effect = self.reg_name in ("DSPAN", "DRMST") and self.params.n <= 2
        if self.reg_name == "RAW" or no_effect:
            pparams = [0]
        else:
            step = self.reg_bound / max(self.reg_steps, 1)
            pparams = [round(step * i, 4) for i in range(self.reg_steps + 1)]

        n_workers = min(j, len(seeds))
        pool_instances = {}
        tree_info = {}

        executor = ProcessPoolExecutor(
            max_workers=n_workers,
            mp_context=multiprocessing.get_context("spawn"),
            initializer=_init_cd_worker,
            initargs=(self, logging.root.level),
        )
        try:
            for pparam in pparams:
                logging.info(
                    f"CD: pparam={pparam}, launching {len(seeds)} seed(s) across {n_workers} worker(s)"
                )
                futures = [
                    executor.submit(
                        _cd_work,
                        i,
                        u,
                        pparam,
                        solver_type,
                        max_iters,
                        max_convergence_iters,
                        timelimit,
                    )
                    for i, u in enumerate(seeds)
                ]
                instances = []
                n_done = 0
                for future in as_completed(futures):
                    try:
                        result = future.result()
                    except Exception as e:
                        logging.error(f"CD worker failed: {e}")
                        executor.shutdown(wait=False, cancel_futures=True)
                        raise RuntimeError(f"CD worker failed: {e}") from e
                    n_done += 1
                    if result is not None:
                        instances.append(result)
                    if n_done % 50 == 0 or n_done == len(futures):
                        logging.info(
                            f"CD: pparam={pparam}, {n_done}/{len(futures)} seeds completed"
                        )

                if not instances:
                    logging.warning(f"CD: no feasible solution at pparam={pparam}")
                    continue
                best = min(instances, key=lambda x: x[0])
                obj_u, cA, cB, u, imf_c, reg_c, t_edges = best
                logging.info(
                    f"CD: pparam={pparam}, best obj=({imf_c:.4f}, {reg_c:.1f}) from {len(instances)} feasible"
                )
                pool_instances[pparam] = [(obj_u, cA, cB, u)]
                if t_edges is not None:
                    tree_info[pparam] = {
                        "tree_edges": t_edges,
                        "total_edge_length": reg_c,
                    }
        finally:
            executor.shutdown(wait=True)

        if not pool_instances:
            raise RuntimeError("Not a single feasible solution found!")
        return pool_instances, tree_info
