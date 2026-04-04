from copy import copy
import logging
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd
from pyomo import environ as pe

from hatchet.compute_cn.solve.base_solver import BaseSolver
from hatchet.compute_cn.solve.ilp_subset import ILPSubset
from hatchet.compute_cn.solve.utils import Random


class Worker:
    """Runs one coordinate-descent restart at a fixed pparam.

    C-step: ILPSubset in CARCH mode.
    U-step: ILPSubset in UARCH mode.
    """

    def __init__(self, work_id, ilp, solver_type, solver_threads=None):
        self.work_id = work_id
        self.ilp = ilp
        self.solver_type = solver_type
        self._solver = BaseSolver._create_solver(solver_type, threads=solver_threads)

    def run(
        self,
        cA,
        cB,
        u,
        pparam,
        max_iters,
        max_convergence_iters,
        tol=0.001,
        timelimit=None,
    ):
        _iters = _convergence_iters = 0
        _u = u
        _cA, _cB = cA, cB
        _prev_obj_u = None
        _imf_c = 0.0
        _reg_c = 0.0
        _tree_edges = None

        while (_iters < max_iters) and (_convergence_iters < max_convergence_iters):
            # C-step: fix u, optimize cA/cB
            carch = copy(self.ilp)
            carch.fix_u(_u)
            carch.create_model()
            carch.hot_start(_cA, _cB)
            carch.model.pparam = pparam
            result = carch.run(
                solver_type=self.solver_type,
                timelimit=timelimit,
                solver=self._solver,
            )
            if result is None:
                logging.debug(
                    f"worker {self.work_id}: C-step infeasible at pparam={pparam}"
                )
                return None
            _, _cA, _cB, _ = result
            _imf_c = pe.value(carch.model.obj_imf)
            _reg_c = pe.value(carch.model.obj_reg)
            _tree_edges = carch.get_tree_edges()

            # U-step: fix cA/cB, optimize u
            uarch = copy(self.ilp)
            uarch.fix_c(_cA, _cB)
            uarch.create_model()
            uarch_results = uarch.run(
                self.solver_type,
                timelimit=timelimit,
                solver=self._solver,
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

        return _obj_u, _cA, _cB, _u, _imf_c, _reg_c, _tree_edges


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


def _work(work_id, u, pparam, solver_type, max_iters, max_convergence_iters, timelimit):
    """Entry point for a single process-pool worker at a fixed pparam."""
    global _cd_worker
    if _cd_worker is None:
        _cd_worker = Worker(
            work_id,
            _cd_global.ilp,
            solver_type,
            solver_threads=_cd_global.solver_threads,
        )
    return _cd_worker.run(
        _cd_global.hcA,
        _cd_global.hcB,
        u,
        pparam=pparam,
        max_iters=max_iters,
        max_convergence_iters=max_convergence_iters,
        tol=_cd_global.cd_tol,
        timelimit=timelimit,
    )


class CoordinateDescent:
    """Coordinate-descent solver for copy-number deconvolution.

    Sweeps pparam over [0, step, 2*step, ..., steps*step]. For each
    pparam, runs all seed restarts in parallel and keeps the best solution.
    """

    def __init__(
        self,
        fcn_data,
        n,
        minprop,
        max_ncns_seg,
        cn_max,
        cn,
        w,
        purities,
        ampdel=True,
        reg_term=None,
        reg_steps=0,
        reg_bound=0.3,
        base=1,
        solve_mode="cd",
        u_init_method="dirichlet",
        u_dir_alpha=0.3,
        solver_threads=None,
        max_degree=3,
        balanced_clusters=None,
        mrca=False,
        zero_cn_thres=0.005,
        tol=0.001,
        cd_tol=0.001,
    ):
        self.reg_name = reg_term if reg_term is not None else "RAW"
        self.reg_steps = reg_steps
        self.reg_bound = reg_bound
        self.u_init_method = u_init_method
        self.u_dir_alpha = u_dir_alpha
        self.solver_threads = solver_threads
        self.cd_tol = cd_tol

        self.ilp = ILPSubset(
            n=n,
            cn_max=cn_max,
            max_ncns_seg=max_ncns_seg,
            minprop=minprop,
            ampdel=ampdel,
            copy_numbers=cn,
            fcn_data=fcn_data,
            w=w,
            purities=purities,
            penalty_param=[self.reg_name, 0.0],
            base=base,
            balanced_clusters=balanced_clusters,
            mrca=mrca,
            zero_cn_thres=zero_cn_thres,
            tol=tol,
            max_degree=max_degree,
        )
        self.ilp.create_model(pprint=True)
        self.hcA, self.hcB = self.ilp.first_hot_start()

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

        # Build regularization path: pparam ∈ [0, reg_bound]
        no_effect = self.reg_name in ("DSPAN", "DRMST") and self.ilp.n <= 2
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
                        _work,
                        i,
                        u,
                        pparam,
                        solver_type,
                        max_iters,
                        max_convergence_iters,
                        timelimit,
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
                obj_u, cA, cB, u, imf_c, reg_c, t_edges = best
                logging.info(
                    f"CD: pparam={pparam}, best obj=({imf_c:.4f}, {reg_c:.1f}) "
                    f"from {len(instances)} feasible"
                )
                pool_instances[pparam] = [(obj_u, cA, cB, u)]
                if t_edges is not None:
                    tree_info[pparam] = {
                        "tree_edges": t_edges,
                        "total_edge_length": reg_c,
                    }
        finally:
            executor.shutdown(wait=True)

        if len(pool_instances) == 0:
            raise RuntimeError("Not a single feasible solution found!")

        return pool_instances, tree_info
