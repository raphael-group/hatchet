from copy import copy
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, as_completed

import pandas as pd
from pyomo import environ as pe

from hatchet.compute_cn.solve.ilp_subset import ILPSubset
from hatchet.compute_cn.solve.utils import Random, model_selection_instance


class Worker:
    """Runs one coordinate-descent restart from a fixed initial u.

    Each Worker alternates between a C-step (fix u, optimise cA/cB) and a
    U-step (fix cA/cB, optimise u) until convergence or the iteration budget
    is exhausted. When a regularization path is active, the C-step sweeps over
    λ values and uses instance-level model selection to pick the best cA/cB
    before proceeding to the U-step.
    """

    def __init__(
        self,
        work_id: int,
        ilp: ILPSubset,
        reg_name: str,
        reg_steps: int,
        reg_ssize: float,
        solver_type: str,
    ):
        self.work_id = work_id
        self.ilp = ilp
        self.reg_name = reg_name
        self.reg_steps = reg_steps
        self.reg_ssize = reg_ssize
        self.solver_type = solver_type
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
            # C-step: fix u, optimize cA/cB
            carch = copy(self.ilp)
            carch.fix_u(_u)
            carch.create_model()
            carch.hot_start(_cA, _cB)

            if self.reg_name != "RAW":
                # Regularization path: sweep λ values, warm-starting each from the previous.
                # DMRCA_SUM only penalises clones at index >= 2; with n <= 2 there are no
                # subclonal clones beyond the MRCA, so a single unregularised solve suffices.
                dmrca_no_effect = self.reg_name == "DMRCA_SUM" and self.ilp.n <= 2
                effective_reg_steps = 0 if dmrca_no_effect else self.reg_steps

                carch_instances = {}
                prev_pparam = None
                for i0 in range(0, effective_reg_steps + 1):
                    pparam = self.reg_ssize * i0
                    carch.model.pparam = pparam
                    if prev_pparam is not None:
                        # carch_instances[pparam] is a single-element list [result];
                        # index [0] unwraps it to the (obj, cA, cB, u) tuple.
                        _, prev_cA, prev_cB, _ = carch_instances[prev_pparam][0]
                        carch.hot_start(prev_cA, prev_cB)
                    result = carch.run(
                        solver_type=self.solver_type,
                        timelimit=timelimit,
                        solver=self._solver,
                    )
                    if result is None:
                        return None
                    carch_instances[pparam] = [result]
                    prev_pparam = pparam

                # Pass solve_mode="cd" so model_selection_instance uses the
                # unregularised error formula (errv = tobj - imf_obj).
                # outdir=None suppresses TSV/PNG output during the inner CD loop.
                best_result, _imf_obj, _selected_key = model_selection_instance(
                    self.ilp.f_a,
                    self.ilp.f_b,
                    self.ilp.w,
                    carch_instances,
                    self.reg_name,
                    "cd",
                    None,
                )
                _obj_c, _cA, _cB, _ = best_result
            else:
                result = carch.run(
                    self.solver_type,
                    timelimit=timelimit,
                    solver=self._solver,
                )
                if result is None:
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
                return None
            _obj_u, _, _, _u = uarch_results

            # Convergence: compare U-step objective across iterations
            # (U-step is always unregularized, so _obj_u is comparable across iters)
            if _prev_obj_u is not None:
                delta = abs(_obj_u - _prev_obj_u)
                if delta < tol:
                    _convergence_iters += 1
                else:
                    _convergence_iters = 0
            _prev_obj_u = _obj_u
            _iters += 1

        return _obj_u, _cA, _cB, _u


# Global reference set by process pool initializer — avoids pickling CoordinateDescent per task
_cd_global = None


def _init_worker(cd):
    global _cd_global
    _cd_global = cd


def _work(work_id, u, solver_type, max_iters, max_convergence_iters, timelimit):
    """Entry point for a single process-pool worker.

    Reads the shared CoordinateDescent state from the module-level
    ``_cd_global`` reference set by ``_init_worker``, constructs a Worker,
    and runs one full coordinate-descent restart.

    Args:
        work_id: Integer index identifying this restart (used for logging).
        u: Initial mixture-proportion matrix for this restart.
        solver_type: Pyomo solver name (e.g. "gurobi" or "cbc").
        max_iters: Maximum number of C/U alternation iterations.
        max_convergence_iters: Number of consecutive non-improving iterations
            before early stopping.
        timelimit: Per-solve wall-clock time limit in seconds, or None.

    Returns:
        Tuple ``(obj, cA, cB, u)`` for the converged solution, or None if
        no feasible solution was found.
    """
    worker = Worker(
        work_id,
        _cd_global.ilp,
        _cd_global.reg_name,
        _cd_global.reg_steps,
        _cd_global.reg_stepsize,
        solver_type,
    )
    return worker.run(
        _cd_global.hcA,
        _cd_global.hcB,
        u,
        max_iters=max_iters,
        max_convergence_iters=max_convergence_iters,
        timelimit=timelimit,
    )


class CoordinateDescent:
    """Coordinate-descent solver for copy-number deconvolution.

    Generates ``n_seed`` random initial mixture proportions (u), then runs
    each as an independent Worker restart in a process pool. All restarts
    share the same ILPSubset template (stored as ``self.ilp``) via the
    module-level ``_cd_global`` initializer, avoiding per-task pickling of
    the full Pyomo model. After all workers finish, solutions are collected,
    sorted by objective value, and returned as an ordered dict keyed by rank.
    """

    def __init__(
        self,
        f_a,
        f_b,
        n,
        minprop,
        max_ncns_seg,
        cn_max,
        cn,
        w,
        purities,
        ampdel=True,
        copy_numbers_fixed=None,
        reg_term=None,
        reg_steps=0,
        reg_stepsize=0.0,
        base=1,
    ):
        self.reg_name = reg_term if reg_term is not None else "RAW"
        self.reg_steps = reg_steps
        self.reg_stepsize = reg_stepsize
        # ilp attribute used here as a convenient storage container for properties
        self.ilp = ILPSubset(
            n=n,
            cn_max=cn_max,
            max_ncns_seg=max_ncns_seg,
            minprop=minprop,
            ampdel=ampdel,
            copy_numbers=cn,
            f_a=f_a,
            f_b=f_b,
            w=w,
            purities=purities,
            copy_numbers_fixed=copy_numbers_fixed,
            penalty_param=[self.reg_name, 0.0],
            base=base,
        )
        # Building the model here is not strictly necessary, as, during execution,
        #   self.carch and c.uarch will copy self.ilp and create+run those models.
        # However, we do so here simply so we can print out some diagnostic information once for the user.
        self.ilp.create_model(pprint=True)
        self.hcA, self.hcB = self.ilp.first_hot_start()

        self.seeds = None

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

        instances = []  # obj. value => (cA, cB, u) mapping
        to_do = []
        with ProcessPoolExecutor(
            max_workers=min(j, n_seed),
            mp_context=multiprocessing.get_context("spawn"),
            initializer=_init_worker,
            initargs=(self,),
        ) as executor:
            for i, u in enumerate(seeds):
                future = executor.submit(
                    _work,
                    i,
                    u,
                    solver_type,
                    max_iters,
                    max_convergence_iters,
                    timelimit,
                )
                to_do.append(future)

            for future in as_completed(to_do):
                instance = future.result()
                if instance is not None:
                    obj, cA, cB, u = instance
                    instances.append((obj, cA, cB, u))

        if len(instances) == 0:
            raise RuntimeError("Not a single feasible solution found!")

        instances_s = {}
        for idx, instance in enumerate(sorted(instances, key=lambda elem: elem[0])):
            instances_s[idx] = instance

        return instances_s
