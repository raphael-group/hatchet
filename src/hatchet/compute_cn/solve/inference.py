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
    SolverParams,
    SolverInputs,
    build_model,
    first_hot_start,
    hot_start,
    build_random_u,
    update_fixed_u,
    update_fixed_cn,
)
from hatchet.compute_cn.solve.utils import (
    dedup_pool_instances,
    split_by_chromosome,
)

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


def extract_solution(
    model,
    params: SolverParams,
    inputs: SolverInputs,
    fixed_u=None,
    fixed_cA=None,
    fixed_cB=None,
):
    """Extract solution dict from solved model.

    Returns dict with keys: imf_obj, cA, cB, u.
    """
    if fixed_cA is not None:
        cA = fixed_cA
        cB = fixed_cB
    else:
        cA = [
            [int(round(model.cA[_m, _n].value)) for _n in range(params.n)]
            for _m in range(inputs.m)
        ]
        cB = [
            [int(round(model.cB[_m, _n].value)) for _n in range(params.n)]
            for _m in range(inputs.m)
        ]
    if fixed_u is not None:
        u = fixed_u
    else:
        u = [
            [model.u[_n, _k].value for _k in range(inputs.k)] for _n in range(params.n)
        ]
    return {"imf_obj": model.obj(), "cA": cA, "cB": cB, "u": u}


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


def extract_pool_solutions(
    solver, model, params: SolverParams, inputs: SolverInputs, pool_size=10
):
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
        cA = [[0] * params.n for _ in range(inputs.m)]
        cB = [[0] * params.n for _ in range(inputs.m)]
        u = [[0.0] * inputs.k for _ in range(params.n)]
        for _m in range(inputs.m):
            for _n in range(params.n):
                for arr, cX in [(cA, model.cA), (cB, model.cB)]:
                    pv = cX[_m, _n]
                    gv = var_map.get(id(pv))
                    arr[_m][_n] = int(round(gv.Xn)) if gv else int(round(pv.value))
        for _n in range(params.n):
            for _k in range(inputs.k):
                pv = model.u[_n, _k]
                gv = var_map.get(id(pv))
                u[_n][_k] = gv.Xn if gv else pv.value
        solutions.append((grb_model.PoolObjVal, cA, cB, u))
    return solutions


def run_full_ilp(
    params,
    inputs,
    reg_steps,
    reg_bound,
    solver_type,
    timelimit,
    warm_start_cA=None,
    warm_start_cB=None,
):
    """Build model once, solve across regularization path.

    Returns pool_instances dict.
    """
    model, var_z = build_model("FULL", params, inputs)

    solver = create_solver(solver_type)

    if warm_start_cA is not None:
        hot_start(model, params, inputs, warm_start_cA, warm_start_cB)

    no_effect = params.reg_name in ("DBOX_L1", "DBOX_L0") and params.n <= 2
    effective_steps = 0 if no_effect else reg_steps

    first_sol = None
    pool_instances = {}

    for i0 in range(effective_steps + 1):
        pparam = reg_bound * i0 / max(effective_steps, 1)
        model.pparam = pparam
        if i0 > 0 and first_sol is not None:
            hot_start(model, params, inputs, first_sol["cA"], first_sol["cB"])

        if not solve_model(model, solver, True, timelimit):
            raise RuntimeError(f"ILP infeasible at pparam={pparam}")
        sol = extract_solution(model, params, inputs)
        sol["imf_obj"] = pe.value(model.obj_imf)
        sol["reg_obj"] = pe.value(model.obj_reg)
        sol_id = f"p{pparam:.4f}_s0"
        pool_instances[sol_id] = sol
        if first_sol is None:
            first_sol = sol

    pool_instances = dedup_pool_instances(pool_instances)
    return pool_instances


# ── Coordinate descent (CD and CNT-CD) ─────────────────────────────────

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


def _get_solver(solver_type):
    global _cd_solver_cache
    if _cd_solver_cache is None:
        _cd_solver_cache = create_solver(
            solver_type, threads=_cd_global["solver_threads"]
        )
    return _cd_solver_cache


def _cd_work(
    work_id, u, pparam, solver_type, max_iters, max_convergence_iters, timelimit
):
    cfg = _cd_global
    params = cfg["params"]
    inputs = cfg["inputs"]
    solver = _get_solver(solver_type)

    _u, _cA, _cB = u, cfg["hcA"], cfg["hcB"]

    model_c, var_z_c = build_model("CARCH", params, inputs, fixed_u=_u)
    model_u, _ = build_model("UARCH", params, inputs, fixed_cA=_cA, fixed_cB=_cB)

    _prev_obj_u = None
    _imf_c = _reg_c = 0.0
    _tree_edges = None
    _iters = _conv_iters = 0

    while _iters < max_iters and _conv_iters < max_convergence_iters:
        if _iters > 0:
            update_fixed_u(model_c, _u, params, inputs)
        hot_start(model_c, params, inputs, _cA, _cB)
        model_c.pparam = pparam
        if not solve_model(model_c, solver, True, timelimit):
            return None
        c_sol = extract_solution(model_c, params, inputs, fixed_u=_u)
        _cA, _cB = c_sol["cA"], c_sol["cB"]
        _imf_c = pe.value(model_c.obj_imf)
        _reg_c = pe.value(model_c.obj_reg)
        _tree_edges = extract_tree_edges(var_z_c, params.n)

        update_fixed_cn(model_u, _cA, _cB, params, inputs)
        if not solve_model(model_u, solver, False, timelimit):
            return None
        u_sol = extract_solution(model_u, params, inputs, fixed_cA=_cA, fixed_cB=_cB)
        _obj_u = u_sol["imf_obj"]
        _u = u_sol["u"]

        if _prev_obj_u is not None:
            if abs(_obj_u - _prev_obj_u) < cfg["cd_tol"]:
                _conv_iters += 1
            else:
                _conv_iters = 0
        _prev_obj_u = _obj_u
        _iters += 1

    return {
        "imf_obj": _obj_u,
        "reg_obj": _reg_c,
        "cA": _cA,
        "cB": _cB,
        "u": _u,
    }


def _cnt_cd_work(
    work_id, u_init, tree_idx, solver_type, max_iters, max_convergence_iters, timelimit
):
    from hatchet.compute_cn.solve.cnt_model import (
        build_c_step_model,
        solve_c_step,
        extract_c_step,
        build_u_step_model,
        extract_u_step,
    )

    cfg = _cd_global
    params = cfg["params"]
    inputs = cfg["inputs"]
    tree = cfg["trees"][tree_idx]
    chrom_groups = cfg["chrom_groups"]
    solver = _get_solver(solver_type)

    n = tree.n
    u = u_init.copy()
    prev_tuple = None
    conv_iters = 0

    for iteration in range(max_iters):
        all_ab = {}
        F_stage1 = F_actual = T_tree = 0.0
        failed = False

        for cg in chrom_groups:
            model, aux = build_c_step_model(tree, params, inputs, u, cg)
            result = solve_c_step(
                model, aux, solver, params.eps_fit, timelimit=timelimit
            )
            if result is None:
                failed = True
                break
            F_stage1 += result["F_star"]
            F_actual += result["F_actual"]
            T_tree += result["T_star"]
            all_ab[id(cg)] = (cg, extract_c_step(model, tree, aux))

        if failed:
            return None

        S, V = inputs.m, tree.n_nodes
        n_te = len(tree.tumor_edges)
        a_full = np.zeros((S, V + 1))
        b_full = np.zeros((S, V + 1))
        event_keys = [
            "alpha_a",
            "alpha_b",
            "delta_a",
            "delta_b",
            "abar_a",
            "abar_b",
            "dbar_a",
            "dbar_b",
        ]
        events_full = {k: np.zeros((S, n_te)) for k in event_keys}

        for cg, ab_data in all_ab.values():
            for local_s, global_s in enumerate(cg):
                a_full[global_s] = ab_data["a"][local_s]
                b_full[global_s] = ab_data["b"][local_s]
                for k in event_keys:
                    events_full[k][global_s] = ab_data[k][local_s]

        a_leaves = a_full[:, 1 : n + 1]
        b_leaves = b_full[:, 1 : n + 1]

        u_model, _ = build_u_step_model(params, inputs, a_leaves, b_leaves)
        u_result = solver.solve(u_model, tee=False)
        if u_result.solver.termination_condition != pe.TerminationCondition.optimal:
            return None
        u = extract_u_step(u_model, params, inputs)

        obj_tuple = (F_stage1, F_actual, T_tree)
        if prev_tuple is not None:
            if all(abs(a - b) < cfg["cd_tol"] for a, b in zip(obj_tuple, prev_tuple)):
                conv_iters += 1
            else:
                conv_iters = 0
        prev_tuple = obj_tuple

        if conv_iters >= max_convergence_iters:
            break

    return {
        "imf_obj": F_actual,
        "reg_obj": 0.0,
        "imf_obj_stage1": F_stage1,
        "tree_obj": T_tree,
        "cA": np.round(a_leaves).astype(int).tolist(),
        "cB": np.round(b_leaves).astype(int).tolist(),
        "u": u.tolist(),
        "a_all": a_full,
        "b_all": b_full,
        "events": events_full,
    }


def run_coordinate_descent(
    params,
    inputs,
    mode="cd",
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
    tree_file=None,
):
    """Run coordinate descent with parallel restarts.

    mode="cd": ILP-based CD over a regularization path.
    mode="cnt_cd": tree-based CNT-CD over enumerated tree shapes.

    Returns pool_instances dict.
    """
    with Random(random_seed):
        seeds = [
            build_random_u(params, inputs, method=u_init_method, alpha=u_dir_alpha)
            for _ in range(n_seed)
        ]

    if u0_tsv_path is not None:
        rows = []
        for restart, u in enumerate(seeds):
            for clone in range(u.shape[0]):
                row = {"restart": restart, "clone": clone}
                for j_idx, sid in enumerate(inputs.sample_ids):
                    row[sid] = u[clone, j_idx]
                rows.append(row)
        pd.DataFrame(rows).to_csv(u0_tsv_path, sep="\t", index=False)

    # Build work groups and config
    cd_config = {
        "params": params,
        "inputs": inputs,
        "solver_threads": solver_threads,
        "cd_tol": cd_tol,
    }

    if mode == "cnt_cd":
        from hatchet.compute_cn.solve.cnt_tree import (
            enumerate_binary_trees,
            parse_newick,
        )

        if tree_file is not None:
            with open(tree_file) as f:
                trees = [parse_newick(f.read())]
        else:
            trees = enumerate_binary_trees(params.n)

        chrom_groups = split_by_chromosome(inputs)
        cd_config["trees"] = trees
        cd_config["chrom_groups"] = chrom_groups

        work_groups = [(f"t{ti}", ti) for ti in range(len(trees))]
        worker_fn = _cnt_cd_work

        logging.info(
            f"CNT-CD: {len(trees)} tree shape(s), {n_seed} seed(s), "
            f"{len(chrom_groups)} chromosome(s), n={params.n}"
        )
    else:
        hcA, hcB = first_hot_start(params, inputs)
        cd_config["hcA"] = hcA
        cd_config["hcB"] = hcB

        no_effect = params.reg_name in ("DBOX_L1", "DBOX_L0") and params.n <= 2
        if params.reg_name == "RAW" or no_effect:
            pparams = [0]
        else:
            step = reg_bound / max(reg_steps, 1)
            pparams = [round(step * i, 4) for i in range(reg_steps + 1)]

        work_groups = [(f"p{pp:.4f}_s0", pp) for pp in pparams]
        worker_fn = _cd_work

    n_workers = min(j, len(seeds))
    pool_instances = {}

    executor = ProcessPoolExecutor(
        max_workers=n_workers,
        mp_context=multiprocessing.get_context("spawn"),
        initializer=_init_cd_worker,
        initargs=(cd_config, logging.root.level),
    )
    try:
        for sol_id, group_key in work_groups:
            logging.info(
                f"{mode.upper()}: group={sol_id}, launching {len(seeds)} seed(s) "
                f"across {n_workers} worker(s)"
            )
            futures = [
                executor.submit(
                    worker_fn,
                    i,
                    u,
                    group_key,
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
                    logging.error(f"{mode.upper()} worker failed: {e}")
                    executor.shutdown(wait=False, cancel_futures=True)
                    raise RuntimeError(f"{mode.upper()} worker failed: {e}") from e
                n_done += 1
                if result is not None:
                    instances.append(result)
                if n_done % 50 == 0 or n_done == len(futures):
                    logging.info(f"  {n_done}/{len(futures)} seeds completed")

            if not instances:
                logging.warning(f"{mode.upper()}: no feasible solution for {sol_id}")
                continue

            best = min(instances, key=lambda x: x["imf_obj"])
            sol_dict = {
                "imf_obj": best["imf_obj"],
                "reg_obj": best["reg_obj"],
                "cA": best["cA"],
                "cB": best["cB"],
                "u": best["u"],
            }

            if mode == "cnt_cd":
                tree = trees[group_key]
                sol_dict.update(
                    imf_obj_stage1=best["imf_obj_stage1"],
                    tree_obj=best["tree_obj"],
                    tree=tree.label(best["a_all"], best["b_all"], best["events"]),
                )

            logging.info(
                f"  {sol_id}: imf={sol_dict['imf_obj']:.4f} reg={sol_dict['reg_obj']:.1f} "
                f"from {len(instances)} feasible"
            )
            pool_instances[sol_id] = sol_dict
    finally:
        executor.shutdown(wait=True)

    if not pool_instances:
        raise RuntimeError(f"{mode.upper()}: no feasible solution found")

    if mode != "cnt_cd":
        pool_instances = dedup_pool_instances(pool_instances)
    return pool_instances
