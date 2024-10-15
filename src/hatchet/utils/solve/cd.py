from copy import copy
from concurrent.futures import ProcessPoolExecutor, as_completed

from hatchet.utils.solve.ilp_subset import ILPSubset, ILPSubsetSplit
from hatchet.utils.solve.utils import Random


class Worker:
    def __init__(self, ilp, solver):
        self.ilp = ilp
        self.solver_type = solver

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
        _cA, _cB = cA, cB  # first hot-start values

        while (_iters < max_iters) and (_convergence_iters < max_convergence_iters):
            carch = copy(self.ilp)
            carch.fix_u(_u)
            carch.create_model()
            carch.hot_start(_cA, _cB)
            carch_results = carch.run(self.solver_type, timelimit=timelimit)
            if carch_results is None:
                return None
            _obj_c, _cA, _cB, _, _, _ = carch_results

            uarch = copy(self.ilp)
            uarch.fix_c(_cA, _cB)
            uarch.create_model()
            uarch_results = uarch.run(self.solver_type, timelimit=timelimit)
            if uarch_results is None:
                return None
            _obj_u, _, _, _u, _, _ = uarch_results

            delta = abs(_obj_c - _obj_u)
            if delta < tol:
                _convergence_iters += 1
            else:
                _convergence_iters = 0

            _iters += 1

        return _obj_u, _cA, _cB, _u


# Top-level 'work' function that can be pickled for multiprocessing
def _work(cd, u, solver_type, max_iters, max_convergence_iters, timelimit):
    worker = Worker(cd.ilp, solver_type)
    return worker.run(
        cd.hcA,
        cd.hcB,
        u,
        max_iters=max_iters,
        max_convergence_iters=max_convergence_iters,
        timelimit=timelimit,
    )


class CoordinateDescent:
    def __init__(self, f_a, f_b, n, mu, d, cn_max, cn, w, purities, ampdel=True):
        # ilp attribute used here as a convenient storage container for properties
        self.ilp = ILPSubset(
            n=n,
            cn_max=cn_max,
            d=d,
            mu=mu,
            ampdel=ampdel,
            copy_numbers=cn,
            f_a=f_a,
            f_b=f_b,
            w=w,
            purities=purities,
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
    ):
        with Random(random_seed):
            seeds = [self.ilp.build_random_u() for _ in range(n_seed)]

        result = {}  # obj. value => (cA, cB, u) mapping
        to_do = []
        with ProcessPoolExecutor(max_workers=min(j, n_seed)) as executor:
            for u in seeds:
                future = executor.submit(
                    _work,
                    self,
                    u,
                    solver_type,
                    max_iters,
                    max_convergence_iters,
                    timelimit,
                )
                to_do.append(future)

            for future in as_completed(to_do):
                results = future.result()
                if results is not None:
                    obj, cA, cB, u = results
                    result[obj] = cA, cB, u

        if not result:
            raise RuntimeError("Not a single feasible solution found!")

        best = min(result)
        return (best,) + result[best] + (self.ilp.cluster_ids, self.ilp.sample_ids)


class CoordinateDescentSplit(CoordinateDescent):
    def __init__(
        self,
        f_a,
        f_b,
        n,
        mu,
        d,
        cn_max,
        cn,
        binsA,
        binsB,
        lengths,
        ampdel=True,
    ):
        # ilp attribute used here as a convenient storage container for properties
        self.ilp = ILPSubsetSplit(
            n=n,
            cn_max=cn_max,
            d=d,
            mu=mu,
            ampdel=ampdel,
            copy_numbers=cn,
            f_a=f_a,
            f_b=f_b,
            binsA=binsA,
            binsB=binsB,
            lengths=lengths,
        )
        # Building the model here is not strictly necessary, as, during execution,
        #   self.carch and c.uarch will copy self.ilp and create+run those models.
        # However, we do so here simply so we can print out some diagnostic information once for the user.
        self.ilp.create_model(pprint=True)
        self.hcA, self.hcB = self.ilp.first_hot_start()

        self.seeds = None
