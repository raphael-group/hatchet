import os
import logging
import pandas as pd
from pyomo import environ as pe

from hatchet.compute_cn.solve.utils import *
from hatchet.compute_cn.solve.ilp_subset import ILPSubset
from hatchet.compute_cn.solve.cd import CoordinateDescent


def solver_available(solver_type: str):
    found_solver = False
    if solver_type == "gurobi":
        found_solver = pe.SolverFactory("gurobi", solver_io="python").available(
            exception_flag=False
        )
    else:
        found_solver = pe.SolverFactory(solver_type).available(exception_flag=False)
    if found_solver:
        logging.info(f"solver={solver_type} is available.")
    return found_solver


def solve(
    f_a: pd.DataFrame,
    f_b: pd.DataFrame,
    n: int,
    minprop: float,
    max_ncns_seg: int,
    cn_max: int,
    weights: pd.Series,
    ampdel: bool,
    clonal: dict,
    purities: dict,
    baf: pd.DataFrame,
    copy_numbers_fixed: dict,
    reg_term: str,
    reg_steps: int,
    reg_stepsize: float,
    solver_type: str,
    solve_mode: str,
    base: int = 1,
    max_iters=10,
    max_convergence_iters=2,
    n_seed=400,
    n_worker=8,
    random_seed=42,
    timelimit=None,
    instances_dir=None,
    verbose=False,
    pool_size=1,
    pool_gap=None,
):
    cd_instances = None
    if solve_mode in ("cd", "both"):
        cd = CoordinateDescent(
            f_a=f_a,
            f_b=f_b,
            n=n,
            minprop=minprop,
            max_ncns_seg=max_ncns_seg,
            cn_max=cn_max,
            w=weights,
            ampdel=ampdel,
            cn=clonal,
            purities=purities,
            copy_numbers_fixed=copy_numbers_fixed,
            reg_term=reg_term,
            reg_steps=reg_steps,
            reg_stepsize=reg_stepsize,
            base=base,
        )

        # obj. value => (cA, cB, u) mapping
        u0_tsv_path = (
            os.path.join(instances_dir, "u0_seeds.tsv")
            if instances_dir is not None
            else None
        )
        cd_instances = cd.run(
            solver_type=solver_type,
            max_iters=max_iters,
            max_convergence_iters=max_convergence_iters,
            n_seed=n_seed,
            j=n_worker,
            random_seed=random_seed,
            timelimit=timelimit,
            u0_tsv_path=u0_tsv_path,
        )
        if instances_dir is not None:
            store_instance_tofile(
                cd_instances,
                f_a,
                f_b,
                baf,
                instances_dir,
                "cd",
                n,
            )

    sol_instances = None
    if solve_mode in ("ilp", "both"):
        sol_instances = {}
        solver = ILPSubset(
            n,
            cn_max,
            max_ncns_seg=max_ncns_seg,
            minprop=minprop,
            ampdel=ampdel,
            copy_numbers=clonal,
            f_a=f_a,
            f_b=f_b,
            w=weights,
            purities=purities,
            copy_numbers_fixed=copy_numbers_fixed,
            penalty_param=[reg_term if reg_term is not None else "RAW", 0.0],
            base=base,
        )
        solver.create_model(pprint=verbose)
        if solve_mode == "both":
            # select local-opt from coordinate-descent instances
            # TODO does the starting point be more useful to do additional model selection?
            _, [obj, cA, cB, _] = min(cd_instances.items(), key=lambda tp: tp[1][0])
            logging.info(f"use CD local opt with obj={obj} to initialize ILP model")
            solver.hot_start(cA, cB)

        pool_instances = {}
        for i0 in range(0, reg_steps + 1):
            if verbose:
                logging.info(f"running instance {i0}/{reg_steps}")
            pparam = reg_stepsize * i0
            solver.model.pparam = pparam
            if i0 > 0:
                cA, cB = sol_instances[0][1:3]
                solver.hot_start(cA, cB)
            sol_instances[pparam] = solver.run(
                solver_type=solver_type,
                timelimit=timelimit,
                pool_size=pool_size,
                pool_gap=pool_gap,
            )
            assert sol_instances[pparam] is not None, "optimization failed"

            if pool_size > 1 and solver_type in ("gurobi", "gurobipy"):
                pool_sols = solver.get_pool_solutions(pool_size=pool_size)
                if pool_sols:
                    pool_instances[pparam] = pool_sols

        if instances_dir is not None:
            store_instance_tofile(
                sol_instances,
                f_a,
                f_b,
                baf,
                instances_dir,
                solve_mode,
                n,
            )
            if pool_instances:
                store_pool_tofile(
                    pool_instances,
                    f_a,
                    f_b,
                    baf,
                    instances_dir,
                    solve_mode,
                    n,
                )

    if solve_mode == "cd":
        return cd_instances
    return sol_instances
