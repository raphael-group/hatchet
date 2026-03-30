from copy import copy
import logging
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd
from pyomo import environ as pe

from hatchet.compute_cn.solve.base_solver import BaseSolver
from hatchet.compute_cn.solve.row_ilp import RowILP
from hatchet.compute_cn.solve.tree_ilp import TreeILP
from hatchet.compute_cn.solve.utils import Random


class Worker:
    """Runs one coordinate-descent restart at a fixed pparam.

    C-step: when reg_term=DRMST uses TreeILP (monolithic tree MILP);
            otherwise uses RowILP in CARCH mode (per-row or monolithic).
    U-step: always RowILP in UARCH mode.
    """

    def __init__(self, work_id, ilp, solver_type, tree_ilp_kwargs=None,
                 solver_threads=None):
        self.work_id = work_id
        self.ilp = ilp
        self.solver_type = solver_type
        self.use_tree = tree_ilp_kwargs is not None
        self._tree_ilp_kwargs = tree_ilp_kwargs
        self._tree_ilp = None
        self._solver = BaseSolver._create_solver(solver_type, threads=solver_threads)

    def run(self, cA, cB, u, pparam, max_iters, max_convergence_iters,
            tol=0.001, timelimit=None):
        _iters = _convergence_iters = 0
        _u = u
        _cA, _cB = cA, cB
        _prev_obj_u = None
        _imf_c = 0.0
        _reg_c = 0.0

        while (_iters < max_iters) and (_convergence_iters < max_convergence_iters):
            # C-step
            if self.use_tree:
                c_result = self._c_step_tree(_cA, _cB, _u, pparam, timelimit)
            else:
                c_result = self._c_step_row(_cA, _cB, _u, pparam, timelimit)

            if c_result is None:
                logging.debug(
                    f"worker {self.work_id}: C-step infeasible at pparam={pparam}"
                )
                return None
            _imf_c, _reg_c, _cA, _cB = c_result

            # U-step: fix cA/cB, optimize u
            uarch = copy(self.ilp)
            uarch.fix_c(_cA, _cB)
            uarch.create_model()
            uarch_results = uarch.run(
                self.solver_type, timelimit=timelimit, solver=self._solver,
            )
            if uarch_results is None:
                logging.debug(f"worker {self.work_id}: U-step infeasible")
                return None
            _obj_u, _, _, _u = uarch_results

            if _prev_obj_u is not None:
                delta = abs(_obj_u - _prev_obj_u)
                if delta < tol:
                    _convergence_iters += 1
                else:
                    _convergence_iters = 0
            _prev_obj_u = _obj_u
            _iters += 1

        return _obj_u, _cA, _cB, _u, _imf_c, _reg_c

    def _c_step_row(self, cA, cB, u, pparam, timelimit):
        """C-step via RowILP (CARCH mode). Returns (imf, reg, cA, cB)."""
        carch = copy(self.ilp)
        carch.fix_u(u)
        carch.create_model()
        carch.hot_start(cA, cB)
        carch.model.pparam = pparam
        result = carch.run(
            solver_type=self.solver_type, timelimit=timelimit,
            solver=self._solver,
        )
        if result is None:
            return None
        obj_c, _cA, _cB, _ = result
        # obj_c is the composite objective; split not available from RowILP
        return obj_c, 0.0, _cA, _cB

    def _c_step_tree(self, cA, cB, u, pparam, timelimit):
        """C-step via TreeILP (monolithic tree MILP). Returns (imf, reg, cA, cB)."""
        if self._tree_ilp is None:
            self._tree_ilp = TreeILP(**self._tree_ilp_kwargs)
            self._tree_ilp.build_model(fixed_u=u)
        else:
            self._tree_ilp.update_u(u)
        self._tree_ilp.update_lambda(pparam)
        self._tree_ilp.hot_start(cA, cB)
        result = self._tree_ilp.run(
            solver_type=self.solver_type, solver=self._solver,
            timelimit=timelimit,
        )
        if result is None:
            return None
        imf_obj, tree_edge_len, _cA, _cB, tree_edges = result
        return imf_obj, tree_edge_len, _cA, _cB


# ---------------------------------------------------------------------------
# Process pool
# ---------------------------------------------------------------------------

_cd_global = None
_cd_worker = None


def _init_worker(cd, log_level):
    global _cd_global, _cd_worker
    _cd_global = cd
    _cd_worker = None
    logging.basicConfig(
        level=logging.ERROR,
        format="%(asctime)s.%(msecs)03d %(levelname)s [worker] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        force=True,
    )
    for name in ("pyomo", "pyomo.core"):
        logging.getLogger(name).setLevel(logging.ERROR)


def _work(work_id, u, pparam, solver_type, max_iters, max_convergence_iters,
          timelimit):
    """Entry point for a single process-pool worker at a fixed pparam."""
    global _cd_worker
    if _cd_worker is None:
        _cd_worker = Worker(
            work_id, _cd_global.ilp, solver_type,
            tree_ilp_kwargs=_cd_global.tree_ilp_kwargs,
            solver_threads=_cd_global.solver_threads,
        )
    return _cd_worker.run(
        _cd_global.hcA, _cd_global.hcB, u,
        pparam=pparam,
        max_iters=max_iters,
        max_convergence_iters=max_convergence_iters,
        timelimit=timelimit,
    )


class CoordinateDescent:
    """Coordinate-descent solver for copy-number deconvolution.

    Sweeps pparam over [0, step, 2*step, ..., steps*step]. For each
    pparam, runs all seed restarts in parallel and keeps the best solution.

    When reg_term=DRMST, the C-step uses a monolithic TreeILP with
    global arborescence topology. Otherwise uses per-row RowILP.
    """

    def __init__(
        self, fcn_data, n, minprop, max_ncns_seg, cn_max, cn, w, purities,
        ampdel=True, reg_term=None, reg_steps=0, reg_stepsize=0.0,
        base=1, solve_mode="cd", u_init_method="dirichlet", u_dir_alpha=0.3,
        solver_threads=None, max_degree=3,
    ):
        self.reg_name = reg_term if reg_term is not None else "RAW"
        self.reg_steps = reg_steps
        self.reg_stepsize = reg_stepsize
        self.u_init_method = u_init_method
        self.u_dir_alpha = u_dir_alpha
        self.solver_threads = solver_threads

        self.ilp = RowILP(
            n=n, cn_max=cn_max, max_ncns_seg=max_ncns_seg,
            minprop=minprop, ampdel=ampdel, copy_numbers=cn,
            fcn_data=fcn_data, w=w, purities=purities,
            penalty_param=[self.reg_name, 0.0], base=base,
        )
        self.ilp.create_model(pprint=True)
        self.hcA, self.hcB = self.ilp.first_hot_start()

        # TreeILP kwargs for DRMST mode (None otherwise)
        if self.reg_name == "DRMST":
            self.tree_ilp_kwargs = dict(
                n=n, cn_max=cn_max, fcn_data=fcn_data, w=w,
                copy_numbers=cn, ampdel=ampdel, base=base, minprop=minprop,
            )
        else:
            self.tree_ilp_kwargs = None

    def run(
        self, solver_type="gurobi", max_iters=10, max_convergence_iters=2,
        n_seed=400, j=8, random_seed=None, timelimit=None, u0_tsv_path=None,
        **_ignored,
    ):
        with Random(random_seed):
            seeds = [
                self.ilp.build_random_u(
                    method=self.u_init_method, dir_alpha=self.u_dir_alpha
                )
                for _ in range(n_seed)
            ]

        if u0_tsv_path is not None:
            sample_ids = list(self.ilp.sample_ids)
            rows = []
            for restart, u in enumerate(seeds):
                for clone in range(u.shape[0]):
                    row = {"restart": restart, "clone": clone}
                    for j_idx, sid in enumerate(sample_ids):
                        row[sid] = u[clone, j_idx]
                    rows.append(row)
            pd.DataFrame(rows).to_csv(u0_tsv_path, sep="\t", index=False)

        # Build regularization path: [0, step, 2*step, ..., steps*step]
        dmrca_no_effect = self.reg_name == "DMRCA_SUM" and self.ilp.n <= 2
        if self.reg_name == "RAW" or dmrca_no_effect:
            pparams = [0]
        else:
            pparams = [self.reg_stepsize * i for i in range(self.reg_steps + 1)]

        n_workers = min(j, len(seeds))
        pool_instances = {}

        executor = ProcessPoolExecutor(
            max_workers=n_workers,
            mp_context=multiprocessing.get_context("spawn"),
            initializer=_init_worker,
            initargs=(self, logging.root.level),
        )
        try:
            for pparam in pparams:
                logging.info(
                    f"CD: pparam={pparam}, launching {len(seeds)} seed(s) "
                    f"across {n_workers} worker(s)"
                )
                to_do = []
                for i, u in enumerate(seeds):
                    future = executor.submit(
                        _work, i, u, pparam, solver_type,
                        max_iters, max_convergence_iters, timelimit,
                    )
                    to_do.append(future)

                instances = []
                n_total = len(to_do)
                n_done = 0
                for future in as_completed(to_do):
                    try:
                        instance = future.result()
                    except Exception as e:
                        logging.error(f"CD worker failed with exception: {e}")
                        executor.shutdown(wait=False, cancel_futures=True)
                        raise RuntimeError(f"CD worker failed: {e}") from e
                    n_done += 1
                    if instance is not None:
                        instances.append(instance)
                    else:
                        logging.debug("CD: worker returned None (infeasible)")
                    if n_done % 50 == 0 or n_done == n_total:
                        logging.info(
                            f"CD: pparam={pparam}, {n_done}/{n_total} seeds completed"
                        )

                if len(instances) == 0:
                    logging.warning(
                        f"CD: no feasible solution at pparam={pparam}, skipping"
                    )
                    continue

                best = min(instances, key=lambda x: x[0])
                obj_u, cA, cB, u, imf_c, reg_c = best
                logging.info(
                    f"CD: pparam={pparam}, best obj=({imf_c:.4f}, {reg_c:.1f}) "
                    f"from {len(instances)} feasible"
                )
                pool_instances[pparam] = [(obj_u, cA, cB, u)]
        finally:
            executor.shutdown(wait=True)

        if len(pool_instances) == 0:
            raise RuntimeError("Not a single feasible solution found!")

        return pool_instances
