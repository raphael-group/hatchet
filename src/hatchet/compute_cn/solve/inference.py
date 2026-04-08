"""Solver execution: run_full_ilp and run_coordinate_descent."""

from __future__ import annotations

import logging
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd
from pyomo import environ as pe
from pyomo.opt import SolverStatus, TerminationCondition

from hatchet.compute_cn.solve.model import (
    build_model,
    first_hot_start,
    hot_start,
    build_random_u,
)
from hatchet.compute_cn.solve.variables import SolverParams

# Stack-based random seeding for reproducibility
_random_states = []


class Random:
    """Context manager that pushes/pops numpy random state for reproducibility."""

    def __init__(self, seed=None):
        self.seed = seed

    def __enter__(self):
        if self.seed is not None:
            _random_states.append(np.random.get_state())
            np.random.set_state(np.random.RandomState(self.seed).get_state())

    def __exit__(self, *args):
        if self.seed is not None:
            np.random.set_state(_random_states.pop())


# ---------------------------------------------------------------------------
# Solver utilities
# ---------------------------------------------------------------------------


def create_solver(solver_type, threads=None):
    """Create a Pyomo solver with suppressed output."""
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


def solve_model(model, solver, warmstart, timelimit):
    """Solve a Pyomo model. Returns True on success."""
    kwargs = {"report_timing": False}
    if timelimit is not None:
        kwargs["timelimit"] = int(timelimit)
    if solver.warm_start_capable():
        kwargs["warmstart"] = warmstart
    results = solver.solve(model, **kwargs)
    ok = (
        results.solver.status == SolverStatus.ok
        and results.solver.termination_condition
        in (TerminationCondition.optimal, TerminationCondition.feasible)
    )
    time_limit = (
        results.solver.status == SolverStatus.aborted
        and results.solver.termination_condition == TerminationCondition.maxTimeLimit
    )
    return ok or time_limit


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
            if var_z[(i, j)].value is not None and round(var_z[(i, j)].value) == 1:
                tree_edges[i] = j
                break
        if i not in tree_edges:
            logging.warning(f"DRMST: clone {i} has no parent, defaulting to root")
            tree_edges[i] = 0
    return tree_edges


def extract_pool_solutions(solver, model, p: SolverParams, pool_size=10):
    """Extract additional solutions from Gurobi's solution pool."""
    try:
        grb_model = solver._solver_model
        var_map = solver._pyomo_var_to_solver_var_map
    except AttributeError:
        return []
    if grb_model.SolCount <= 1:
        return []

    solutions = []
    for sol_idx in range(1, min(grb_model.SolCount, pool_size)):
        grb_model.setParam("SolutionNumber", sol_idx)
        cA = [[0] * p.n for _ in range(p.m)]
        cB = [[0] * p.n for _ in range(p.m)]
        u = [[0.0] * p.k for _ in range(p.n)]
        for _m in range(p.m):
            for _n in range(p.n):
                for arr, cX in [(cA, model.cA), (cB, model.cB)]:
                    pv = cX[_m, _n]
                    gv = var_map.get(id(pv))
                    arr[_m][_n] = int(round(gv.Xn)) if gv else int(round(pv.value))
        for _n in range(p.n):
            for _k in range(p.k):
                pv = model.u[_n, _k]
                gv = var_map.get(id(pv))
                u[_n][_k] = gv.Xn if gv else pv.value
        solutions.append((grb_model.PoolObjVal, cA, cB, u))
    return solutions


# ---------------------------------------------------------------------------
# Full ILP run
# ---------------------------------------------------------------------------


def run_full_ilp(
    params,
    penalty_param,
    reg_steps,
    reg_bound,
    solver_type,
    timelimit,
    pool_size=1,
    pool_gap=None,
    warm_start_cA=None,
    warm_start_cB=None,
):
    """Build model once, solve across regularization path.

    Returns (pool_instances, tree_info) in the same format as run_coordinate_descent.
    """
    model, var_z = build_model(params, penalty_param)

    solver = create_solver(solver_type)
    if pool_size > 1 and solver_type in ("gurobi", "gurobipy"):
        solver.options["PoolSolutions"] = pool_size
        solver.options["PoolSearchMode"] = 0
        if pool_gap is not None:
            solver.options["PoolGap"] = pool_gap

    if warm_start_cA is not None:
        hot_start(model, params, warm_start_cA, warm_start_cB)

    pname = penalty_param[0]
    dmrca_no_effect = pname == "DMRCA_SUM" and params.n <= 2
    effective_steps = 0 if dmrca_no_effect else reg_steps

    sol_instances = {}
    pool_instances = {}
    tree_info = {}

    for i0 in range(effective_steps + 1):
        pparam = reg_bound * i0 / max(effective_steps, 1)
        model.pparam = pparam
        if i0 > 0 and 0 in sol_instances:
            hot_start(model, params, sol_instances[0][1], sol_instances[0][2])

        if not solve_model(model, solver, True, timelimit):
            raise RuntimeError(f"ILP infeasible at pparam={pparam}")
        sol = extract_solution(model, params)
        sol_instances[pparam] = sol
        pool_instances[pparam] = [sol]

        if pool_size > 1 and solver_type in ("gurobi", "gurobipy"):
            pool_sols = extract_pool_solutions(solver, model, params, pool_size)
            if pool_sols:
                pool_instances[pparam].extend(pool_sols)

    t_edges = extract_tree_edges(var_z, params.n)
    if t_edges is not None:
        for pparam in pool_instances:
            tree_info[pparam] = {"tree_edges": t_edges}

    return pool_instances, tree_info


# ---------------------------------------------------------------------------
# Coordinate Descent
# ---------------------------------------------------------------------------

_cd_global = None
_cd_solver_cache = None


def _init_cd_worker(cd_config, log_level):
    global _cd_global, _cd_solver_cache
    _cd_global = cd_config
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
    cfg = _cd_global
    p = cfg["params"]

    if _cd_solver_cache is None:
        _cd_solver_cache = create_solver(solver_type, threads=cfg["solver_threads"])
    solver = _cd_solver_cache

    _u, _cA, _cB = u, cfg["hcA"], cfg["hcB"]
    _prev_obj_u = None
    _imf_c = _reg_c = 0.0
    _tree_edges = None
    _iters = _conv_iters = 0

    while _iters < max_iters and _conv_iters < max_convergence_iters:
        # C-step
        p_c = SolverParams(**{**p.__dict__, "mode": "CARCH"})
        model_c, var_z_c = build_model(p_c, cfg["penalty_param"], fixed_u=_u)
        hot_start(model_c, p_c, _cA, _cB)
        model_c.pparam = pparam
        if not solve_model(model_c, solver, True, timelimit):
            return None
        _, _cA, _cB, _ = extract_solution(model_c, p_c)
        _imf_c = pe.value(model_c.obj_imf)
        _reg_c = pe.value(model_c.obj_reg)
        _tree_edges = extract_tree_edges(var_z_c, p.n)

        # U-step
        p_u = SolverParams(**{**p.__dict__, "mode": "UARCH"})
        model_u, _ = build_model(p_u, cfg["penalty_param"], fixed_cA=_cA, fixed_cB=_cB)
        if not solve_model(model_u, solver, False, timelimit):
            return None
        _obj_u, _, _, _u = extract_solution(model_u, p_u)

        if _prev_obj_u is not None:
            if abs(_obj_u - _prev_obj_u) < cfg["cd_tol"]:
                _conv_iters += 1
            else:
                _conv_iters = 0
        _prev_obj_u = _obj_u
        _iters += 1

    return _obj_u, _cA, _cB, _u, _imf_c, _reg_c, _tree_edges


def run_coordinate_descent(
    params,
    reg_term="RAW",
    reg_steps=0,
    reg_bound=0.3,
    u_init_method="dirichlet",
    u_dir_alpha=0.3,
    solver_threads=None,
    cd_tol=0.001,
    solver_type="gurobi",
    max_iters=10,
    max_convergence_iters=2,
    n_seed=400,
    j=8,
    random_seed=None,
    timelimit=None,
    u0_tsv_path=None,
):
    """Run coordinate descent with parallel restarts over a regularization path.

    Returns (pool_instances, tree_info).
    """
    reg_name = reg_term if reg_term else "RAW"
    penalty_param = [reg_name, 0.0]
    hcA, hcB = first_hot_start(params)

    with Random(random_seed):
        seeds = [
            build_random_u(params, method=u_init_method, alpha=u_dir_alpha)
            for _ in range(n_seed)
        ]

    if u0_tsv_path is not None:
        rows = []
        for restart, u in enumerate(seeds):
            for clone in range(u.shape[0]):
                row = {"restart": restart, "clone": clone}
                for j_idx, sid in enumerate(params.sample_ids):
                    row[sid] = u[clone, j_idx]
                rows.append(row)
        pd.DataFrame(rows).to_csv(u0_tsv_path, sep="\t", index=False)

    no_effect = reg_name in ("DSPAN", "DRMST") and params.n <= 2
    if reg_name == "RAW" or no_effect:
        pparams = [0]
    else:
        step = reg_bound / max(reg_steps, 1)
        pparams = [round(step * i, 4) for i in range(reg_steps + 1)]

    cd_config = {
        "params": params,
        "penalty_param": penalty_param,
        "hcA": hcA,
        "hcB": hcB,
        "solver_threads": solver_threads,
        "cd_tol": cd_tol,
    }

    n_workers = min(j, len(seeds))
    pool_instances = {}
    tree_info = {}

    executor = ProcessPoolExecutor(
        max_workers=n_workers,
        mp_context=multiprocessing.get_context("spawn"),
        initializer=_init_cd_worker,
        initargs=(cd_config, logging.root.level),
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
                tree_info[pparam] = {"tree_edges": t_edges, "total_edge_length": reg_c}
    finally:
        executor.shutdown(wait=True)

    if not pool_instances:
        raise RuntimeError("Not a single feasible solution found!")
    return pool_instances, tree_info
