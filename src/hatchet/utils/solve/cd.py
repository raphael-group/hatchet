from copy import copy
from concurrent.futures import ThreadPoolExecutor, as_completed
from .ilp_subset import ILPSubset
from .utils import Random


class Worker:

    def __init__(self, ilp):
        self.ilp = ilp

    def run(self, cA, cB, u, max_iters, max_convergence_iters, tol=0.001):
        _iters = _convergence_iters = 0
        _obj_u = _cA = _cB = _u = None

        while (_iters < max_iters) and (_convergence_iters < max_convergence_iters):

            carch = copy(self.ilp)
            carch.fix_u(u)
            carch.create_model()
            carch.hot_start(cA, cB)
            _obj_c, _cA, _cB, _ = carch.run()

            uarch = copy(self.ilp)
            uarch.fix_c(_cA, _cB)
            uarch.create_model()
            _obj_u, _, _, _u = uarch.run()

            if abs(_obj_c - _obj_u) < tol:
                _convergence_iters += 1
            else:
                _convergence_iters = 0

            _iters += 1

        return _obj_u, _cA, _cB, _u


class CoordinateDescent:

    def __init__(self, f_a, f_b, n, mu, d, cn_max, cn, w, ampdel=True):
        # ilp attribute used here as a convenient storage container for properties
        self.ilp = ILPSubset(n=n, cn_max=cn_max, d=d, mu=mu, ampdel=ampdel, copy_numbers=cn, f_a=f_a, f_b=f_b, w=w)
        self.hcA, self.hcB = self.ilp.first_hot_start()

        self.seeds = None

    def run(self, max_iters=10, max_convergence_iters=2, n_seed=400, j=8, random_seed=None):
        with Random(random_seed):
            seeds = [self.ilp.build_random_u() for _ in range(n_seed)]

        def _work(u):
            worker = Worker(self.ilp)
            return worker.run(self.hcA, self.hcB, u, max_iters=max_iters, max_convergence_iters=max_convergence_iters)

        result = {}  # obj. value => (cA, cB, u) mapping
        to_do = []
        with ThreadPoolExecutor(max_workers=min(j, n_seed)) as executor:
            for u in seeds:
                future = executor.submit(_work, u)
                to_do.append(future)

            for future in as_completed(to_do):
                obj, cA, cB, u = future.result()
                result[obj] = cA, cB, u

        best = min(result)
        return (best,) + result[best]
