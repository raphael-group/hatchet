from copy import copy
from concurrent.futures import ProcessPoolExecutor, as_completed
from .ilp_subset import ILPSubset
from .utils import Random


class Worker:

    def __init__(self, ilp, solver):
        self.ilp = ilp
        self.solver_type = solver

    def run(self, cA, cB, u, max_iters, max_convergence_iters, tol=0.001):
        _iters = _convergence_iters = 0
        _obj_u = _u = None
        _cA, _cB = cA, cB  # first hot-start values

        while (_iters < max_iters) and (_convergence_iters < max_convergence_iters):

            carch = copy(self.ilp)
            carch.fix_u(u)
            carch.create_model()
            carch.hot_start(_cA, _cB)
            _obj_c, _cA, _cB, _ = carch.run(self.solver_type)

            uarch = copy(self.ilp)
            uarch.fix_c(_cA, _cB)
            uarch.create_model()
            _obj_u, _, _, _u = uarch.run(self.solver_type)

            if abs(_obj_c - _obj_u) < tol:
                _convergence_iters += 1
            else:
                _convergence_iters = 0

            _iters += 1

        return _obj_u, _cA, _cB, _u


# Top-level 'work' function that can be pickled for multiprocessing
def _work(cd, u, solver_type, max_iters, max_convergence_iters):
    worker = Worker(cd.ilp, solver_type)
    return worker.run(cd.hcA, cd.hcB, u, max_iters=max_iters, max_convergence_iters=max_convergence_iters)


class CoordinateDescent:

    def __init__(self, f_a, f_b, n, mu, d, cn_max, cn, w, ampdel=True):
        # ilp attribute used here as a convenient storage container for properties
        self.ilp = ILPSubset(n=n, cn_max=cn_max, d=d, mu=mu, ampdel=ampdel, copy_numbers=cn, f_a=f_a, f_b=f_b, w=w)
        self.hcA, self.hcB = self.ilp.first_hot_start()

        self.seeds = None

    def run(self, solver_type='gurobipy', max_iters=10, max_convergence_iters=2, n_seed=400, j=8, random_seed=None):
        with Random(random_seed):
            seeds = [self.ilp.build_random_u() for _ in range(n_seed)]

        result = {}  # obj. value => (cA, cB, u) mapping
        to_do = []
        with ProcessPoolExecutor(max_workers=min(j, n_seed)) as executor:
            for u in seeds:
                future = executor.submit(_work, self, u, solver_type, max_iters, max_convergence_iters)
                to_do.append(future)

            for future in as_completed(to_do):
                obj, cA, cB, u = future.result()
                result[obj] = cA, cB, u

        best = min(result)
        return (best,) + result[best]
