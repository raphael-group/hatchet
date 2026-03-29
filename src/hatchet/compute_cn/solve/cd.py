from copy import copy
import logging
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, as_completed

import pandas as pd
from pyomo import environ as pe

from hatchet.compute_cn.solve.row_ilp import RowILP
from hatchet.compute_cn.solve.utils import Random


class Worker:
    """Runs one coordinate-descent restart at a fixed regularization lambda.

    Each Worker alternates between a C-step (fix u, optimise cA/cB at the
    given lambda) and a U-step (fix cA/cB, optimise u) until convergence or
    the iteration budget is exhausted.
    """

    def __init__(self, work_id: int, ilp: RowILP, solver_type: str, lexi: bool = False):
        self.work_id = work_id
        self.ilp = ilp
        self.solver_type = solver_type
        self.lexi = lexi
        self._solver = self._create_solver()

    def _create_solver(self):
        if self.solver_type in ("gurobipy", "gurobi"):
            s = pe.SolverFactory("gurobi", solver_io="python")
            s.options["OutputFlag"] = 0
            s.options["LogToConsole"] = 0
            s.options["LogFile"] = ""
        else:
            s = pe.SolverFactory(self.solver_type)
        return s

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

        while (_iters < max_iters) and (_convergence_iters < max_convergence_iters):
            # C-step: fix u, optimize cA/cB at fixed lambda
            carch = copy(self.ilp)
            carch.fix_u(_u)

            if self.lexi:
                # Per-row lexicographic: Level-1 min PI violations, Level-2 min L1
                _obj_c = 0.0
                _new_cA = [None] * carch.m
                _new_cB = [None] * carch.m
                c_infeasible = False
                for _m in range(carch.m):
                    cluster_id = carch.cluster_ids[_m]
                    if cluster_id in carch.copy_numbers:
                        # Fixed-CN cluster: solution is fully determined
                        _cnA, _cnB = carch.copy_numbers[cluster_id]
                        _new_cA[_m] = [1] + [_cnA] * (carch.n - 1)
                        _new_cB[_m] = [1] + [_cnB] * (carch.n - 1)
                        continue
                    carch.create_row_lexi_model(_m)
                    carch.hot_start_row(_cA[_m], _cB[_m])
                    row_result = carch.run_row_lexi(
                        solver_type=self.solver_type,
                        timelimit=timelimit,
                        solver=self._solver,
                        pparam=pparam,
                    )
                    if row_result is None:
                        logging.debug(
                            f"worker {self.work_id}: C-step row {_m} infeasible at lambda={pparam}"
                        )
                        c_infeasible = True
                        break
                    obj_m, cA_row, cB_row = row_result
                    _obj_c += obj_m
                    _new_cA[_m] = cA_row
                    _new_cB[_m] = cB_row
                if c_infeasible:
                    return None
                _cA, _cB = _new_cA, _new_cB
            else:
                # Monolithic C-step
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
                        f"worker {self.work_id}: C-step infeasible at lambda={pparam}"
                    )
                    return None
                _obj_c, _cA, _cB, _ = result

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

            # Convergence: compare U-step objective across iterations
            if _prev_obj_u is not None:
                delta = abs(_obj_u - _prev_obj_u)
                if delta < tol:
                    _convergence_iters += 1
                else:
                    _convergence_iters = 0
            _prev_obj_u = _obj_u
            _iters += 1

        return _obj_u, _cA, _cB, _u


# Global reference set by process pool initializer
_cd_global = None


def _init_worker(cd, log_level):
    global _cd_global
    _cd_global = cd
    logging.basicConfig(
        level=logging.ERROR,
        format="%(asctime)s.%(msecs)03d %(levelname)s [worker] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        force=True,
    )
    for name in ("pyomo", "pyomo.core"):
        logging.getLogger(name).setLevel(logging.ERROR)


def _work(work_id, u, pparam, solver_type, max_iters, max_convergence_iters, timelimit):
    """Entry point for a single process-pool worker at a fixed lambda."""
    worker = Worker(work_id, _cd_global.ilp, solver_type, lexi=_cd_global.lexi)
    return worker.run(
        _cd_global.hcA,
        _cd_global.hcB,
        u,
        pparam=pparam,
        max_iters=max_iters,
        max_convergence_iters=max_convergence_iters,
        timelimit=timelimit,
    )


class CoordinateDescent:
    """Coordinate-descent solver for copy-number deconvolution.

    For each regularization lambda in [0, delta, 2*delta, ...], runs all seed
    restarts in parallel and keeps the best solution per lambda.  Returns
    results in the ``{pparam: [best_solution]}`` format, enabling unified
    model selection across the regularization path.
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
        reg_stepsize=0.0,
        base=1,
        solve_mode="cd",
    ):
        self.reg_name = reg_term if reg_term is not None else "RAW"
        self.reg_steps = reg_steps
        self.reg_stepsize = reg_stepsize
        self.lexi = solve_mode == "cd_lexi"
        self.ilp = RowILP(
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
    ):
        with Random(random_seed):
            seeds = [self.ilp.build_random_u() for _ in range(n_seed)]

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

        # Build regularization path
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
                    f"CD: lambda={pparam}, launching {len(seeds)} seed(s) across {n_workers} worker(s)"
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
                            f"CD: lambda={pparam}, {n_done}/{n_total} seeds completed"
                        )

                if len(instances) == 0:
                    logging.warning(
                        f"CD: no feasible solution at lambda={pparam}, skipping"
                    )
                    continue

                # Keep only the best restart per lambda
                best = min(instances, key=lambda x: x[0])
                logging.info(
                    f"CD: lambda={pparam}, best obj={best[0]:.4f} from {len(instances)} feasible"
                )
                pool_instances[pparam] = [best]
        finally:
            executor.shutdown(wait=True)

        if len(pool_instances) == 0:
            raise RuntimeError("Not a single feasible solution found!")

        return pool_instances
