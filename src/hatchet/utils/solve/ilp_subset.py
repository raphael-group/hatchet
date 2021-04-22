import math
import numpy as np
import pandas as pd
from pyomo import environ as pe

from .utils import Random


class ILPSubset:
    def __init__(self, n, cn_max, d, mu, ampdel, copy_numbers, f_a, f_b, w):

        # Each ILPSubset maintains its own data, so make a deep-copy of passed-in DataFrames
        f_a, f_b = f_a.copy(deep=True), f_b.copy(deep=True)

        assert f_a.shape == f_b.shape
        assert np.all(f_a.index == f_b.index)
        assert np.all(f_a.columns == f_b.columns)

        self.m, self.k = f_a.shape
        self.f_a = f_a
        self.f_b = f_b
        self.cluster_ids = f_a.index
        self.sample_ids = f_a.columns

        self.n = n
        self.cn_max = cn_max
        self.d = d
        self.mu = mu
        self.ampdel = ampdel
        self.copy_numbers = copy_numbers
        self.w = w

        self.tol = 0.001

        self.mode = 'FULL'

        # Values we want to optimize for, as dataframes
        self.cA = pd.DataFrame(index=self.cluster_ids, columns=range(n))
        self.cB = pd.DataFrame(index=self.cluster_ids, columns=range(n))
        self.u = pd.DataFrame(index=range(n), columns=self.sample_ids)

        # Fixed values of cA/cB/u
        self._fixed_cA = pd.DataFrame(index=self.cluster_ids, columns=range(n))
        self._fixed_cB = pd.DataFrame(index=self.cluster_ids, columns=range(n))
        self._fixed_u = pd.DataFrame(index=range(n), columns=self.sample_ids)

        self.model = None  # initialized on create_model()

    def __copy__(self):
        return ILPSubset(
            n=self.n, cn_max=self.cn_max, d=self.d, mu=self.mu, ampdel=self.ampdel, copy_numbers=self.copy_numbers,
            f_a=self.f_a, f_b=self.f_b, w=self.w
        )

    @property
    def M(self):
        return math.floor(math.log2(self.cn_max)) + 1

    @property
    def base(self):
        return min(2, len(self.copy_numbers))

    @staticmethod
    def symmCoeff(i):
        return math.pow(i + 1, 1)

    @property
    def optimized_cA(self):
        if self.mode == 'UARCH':
            return self._fixed_cA
        else:
            return self.cA

    @property
    def optimized_cB(self):
        if self.mode == 'UARCH':
            return self._fixed_cB
        else:
            return self.cB

    @property
    def optimized_u(self):
        if self.mode == 'CARCH':
            return self._fixed_u
        else:
            return self.u

    def create_model(self):

        m, n, k = self.m, self.n, self.k
        f_a, f_b = self.f_a, self.f_b
        cn_max = self.cn_max
        ampdel = self.ampdel
        copy_numbers = self.copy_numbers
        mode_t = self.mode
        d = self.d
        _M = self.M
        _base = self.base

        model = pe.ConcreteModel()

        fA = {}
        fB = {}
        yA = {}
        yB = {}
        adA = {}
        adB = {}

        for _m in range(m):

            cluster_id = f_a.index[_m]

            # upper bound for solver
            ub = max(sum(copy_numbers.get(cluster_id, (0, 0))), cn_max)

            for _k in range(k):
                yA[(_m, _k)] = pe.Var(bounds=(0, np.inf), domain=pe.Reals)
                model.add_component(f'yA_{_m + 1}_{_k + 1}', yA[(_m, _k)])
                yB[(_m, _k)] = pe.Var(bounds=(0, np.inf), domain=pe.Reals)
                model.add_component(f'yB_{_m + 1}_{_k + 1}', yB[(_m, _k)])
                fA[(_m, _k)] = pe.Var(bounds=(0, ub), domain=pe.Reals)
                model.add_component(f'fA_{_m + 1}_{_k + 1}', fA[(_m, _k)])
                fB[(_m, _k)] = pe.Var(bounds=(0, ub), domain=pe.Reals)
                model.add_component(f'fB_{_m + 1}_{_k + 1}', fB[(_m, _k)])

        if mode_t in ('FULL', 'CARCH'):
            for _m in range(m):

                cluster_id = f_a.index[_m]

                # upper bound for solver
                ub = max(sum(copy_numbers.get(cluster_id, (0, 0))), cn_max)

                for _n in range(n):
                    self.cA.iloc[_m][_n] = pe.Var(bounds=(0, ub), domain=pe.Integers)
                    model.add_component(f'cA_{_m + 1}_{_n + 1}', self.cA.iloc[_m][_n])
                    self.cB.iloc[_m][_n] = pe.Var(bounds=(0, ub), domain=pe.Integers)
                    model.add_component(f'cB_{_m + 1}_{_n + 1}', self.cB.iloc[_m][_n])

            if ampdel:
                for _m in range(m):
                    cluster_id = f_a.index[_m]
                    if cluster_id not in copy_numbers:
                        adA[_m] = pe.Var(bounds=(0, 1), domain=pe.Binary)
                        model.add_component(f'adA_{_m + 1}', adA[_m])
                        adB[_m] = pe.Var(bounds=(0, 1), domain=pe.Binary)
                        model.add_component(f'adB_{_m + 1}', adB[_m])

        bitcA = {}
        bitcB = {}
        if (mode_t == 'FULL') or (d > 0 and mode_t == 'CARCH'):
            for _b in range(_M):
                for _m in range(m):
                    for _n in range(n):
                        bitcA[(_b, _m, _n)] = pe.Var(bounds=(0, 1), domain=pe.Binary)
                        model.add_component(f'bitcA_{_b + 1}_{_m + 1}_{_n + 1}', bitcA[(_b, _m, _n)])
                        bitcB[(_b, _m, _n)] = pe.Var(bounds=(0, 1), domain=pe.Binary)
                        model.add_component(f'bitcB_{_b + 1}_{_m + 1}_{_n + 1}', bitcB[(_b, _m, _n)])

        if mode_t in ('FULL', 'UARCH'):
            for _n in range(n):
                for _k in range(k):
                    self.u.iloc[_n][_k] = pe.Var(bounds=(0, 1), domain=pe.Reals)
                    model.add_component(f'u_{_n + 1}_{_k + 1}', self.u.iloc[_n][_k])

        vA = {}
        vB = {}
        if mode_t == 'FULL':
            for _b in range(_M):
                for _m in range(m):
                    for _n in range(n):
                        for _k in range(k):
                            vA[(_b, _m, _n, _k)] = pe.Var(bounds=(0, 1), domain=pe.Reals)
                            model.add_component(f'vA_{_b + 1}_{_m + 1}_{_n + 1}_{_k + 1}', vA[(_b, _m, _n, _k)])
                            vB[(_b, _m, _n, _k)] = pe.Var(bounds=(0, 1), domain=pe.Reals)
                            model.add_component(f'vB_{_b + 1}_{_m + 1}_{_n + 1}_{_k + 1}', vB[(_b, _m, _n, _k)])

        x = {}
        if (mode_t in ('FULL', 'UARCH')) and (self.mu > 0):
            for _n in range(n):
                for _k in range(k):
                    x[(_n, _k)] = pe.Var(domain=pe.Binary)
                    model.add_component(f'x_{_n + 1}_{_k + 1}', x[(_n, _k)])

        # buildOptionalVariables
        z = {}
        if (mode_t in ('FULL', 'CARCH')) and d > 0:
            for _m in range(self.m):
                for _n in range(1, self.n):
                    for _d in range(d):
                        z[(_m, _n, _d)] = pe.Var(bounds=(0, 1), domain=pe.Binary)
                        model.add_component(f'z_{_m + 1}_{_n + 1}_{_d + 1}', z[(_m, _n, _d)])

        # CONSTRAINTS
        model.constraints = pe.ConstraintList()

        for _m in range(m):

            cluster_id = f_a.index[_m]
            f_a_values = f_a.loc[cluster_id].values
            f_b_values = f_b.loc[cluster_id].values

            for _k in range(k):
                _yA, _yB = yA[(_m, _k)], yB[(_m, _k)]
                _fA, _fB = fA[(_m, _k)], fB[(_m, _k)]

                model.constraints.add(float(f_a_values[_k]) - _fA <= _yA)
                model.constraints.add(_fA - float(f_a_values[_k]) <= _yA)
                model.constraints.add(float(f_b_values[_k]) - _fB <= _yB)
                model.constraints.add(_fB - float(f_b_values[_k]) <= _yB)

        if mode_t == 'FULL':
            for _m in range(m):
                for _k in range(k):

                    sum_a = 0
                    for _n in range(n):
                        for _b in range(_M):
                            sum_a += vA[(_b, _m, _n, _k)] * math.pow(2, _b)
                            model.constraints.add(vA[(_b, _m, _n, _k)] <= bitcA[(_b, _m, _n)])
                            model.constraints.add(vA[(_b, _m, _n, _k)] <= self.u.iloc[_n][_k])
                            model.constraints.add(vA[(_b, _m, _n, _k)] >= bitcA[(_b, _m, _n)] + self.u.iloc[_n][_k] - 1)

                    model.constraints.add(fA[(_m, _k)] == sum_a)

            for _m in range(m):
                for _k in range(k):

                    sum_b = 0
                    for _n in range(n):
                        for _b in range(_M):
                            sum_b += vB[(_b, _m, _n, _k)] * math.pow(2, _b)
                            model.constraints.add(vB[(_b, _m, _n, _k)] <= bitcB[(_b, _m, _n)])
                            model.constraints.add(vB[(_b, _m, _n, _k)] <= self.u.iloc[_n][_k])
                            model.constraints.add(vB[(_b, _m, _n, _k)] >= bitcB[(_b, _m, _n)] + self.u.iloc[_n][_k] - 1)

                    model.constraints.add(fB[(_m, _k)] == sum_b)

            for _n in range(n):
                for _k in range(k):
                    _sum = 0
                    for _m in range(m):
                        for _b in range(_M):
                            _sum += bitcA[(_b, _m, _n)] + bitcB[(_b, _m, _n)]
                    model.constraints.add(_sum >= self.u.iloc[_n][_k])

        if (mode_t == 'FULL') or (d > 0 and mode_t == 'CARCH'):
            for _m in range(m):
                cluster_id = f_a.index[_m]
                # upper bound for solver
                ub = max(sum(copy_numbers.get(cluster_id, (0, 0))), cn_max)

                for _n in range(n):
                    sum_a = 0
                    sum_b = 0
                    for _b in range(_M):
                        sum_a += bitcA[(_b, _m, _n)] * math.pow(2, _b)
                        sum_b += bitcB[(_b, _m, _n)] * math.pow(2, _b)

                    model.constraints.add(self.cA.iloc[_m][_n] == sum_a)
                    model.constraints.add(self.cB.iloc[_m][_n] == sum_b)
                    model.constraints.add(self.cA.iloc[_m][_n] + self.cB.iloc[_m][_n] <= ub)

        if mode_t == 'CARCH':
            # TODO: These loops can be collapsed once validation against C++ is complete
            for _m in range(m):
                for _k in range(k):
                    _sumA = 0
                    _sumB = 0
                    for _n in range(n):
                        if self._fixed_u.iloc[_n][_k] >= self.mu - self.tol:
                            _sumA += self.cA.iloc[_m][_n] * self._fixed_u.iloc[_n][_k]
                            _sumB += self.cB.iloc[_m][_n] * self._fixed_u.iloc[_n][_k]
                    model.constraints.add(fA[(_m, _k)] == _sumA)
                    model.constraints.add(fB[(_m, _k)] == _sumB)

                cluster_id = f_a.index[_m]
                # upper bound for solver
                ub = max(sum(copy_numbers.get(cluster_id, (0, 0))), cn_max)
                for _n in range(n):
                    model.constraints.add(self.cA.iloc[_m][_n] + self.cB.iloc[_m][_n] <= ub)

        if mode_t in ('FULL', 'CARCH'):
            for _m in range(m):
                model.constraints.add(self.cA.iloc[_m][0] == 1)
                model.constraints.add(self.cB.iloc[_m][0] == 1)

            if ampdel:
                for _m in range(m):
                    cluster_id = f_a.index[_m]
                    if cluster_id not in copy_numbers:
                        for _n in range(1, n):
                            model.constraints.add(self.cA.iloc[_m][_n] <= cn_max * adA[_m] + _base - _base * adA[_m])
                            model.constraints.add(self.cA.iloc[_m][_n] >= _base * adA[_m])
                            model.constraints.add(self.cB.iloc[_m][_n] <= cn_max * adB[_m] + _base - _base * adB[_m])
                            model.constraints.add(self.cB.iloc[_m][_n] >= _base * adB[_m])

        if mode_t == 'UARCH':
            # TODO: These loops can be collapsed once validation against C++ is complete
            for _m in range(m):
                for _k in range(k):
                    _sumA = 0
                    _sumB = 0
                    for _n in range(n):
                        _sumA += int(self._fixed_cA.iloc[_m][_n]) * self.u.iloc[_n][_k]
                        _sumB += int(self._fixed_cB.iloc[_m][_n]) * self.u.iloc[_n][_k]
                    model.constraints.add(fA[(_m, _k)] == _sumA)
                    model.constraints.add(fB[(_m, _k)] == _sumB)

        if mode_t in ('FULL', 'UARCH'):
            for _k in range(k):
                _sum = 0
                for _n in range(n):
                    _sum += self.u.iloc[_n][_k]
                model.constraints.add(_sum == 1)

        if (mode_t in ('FULL', 'UARCH')) and self.mu > 0:
            for _k in range(k):
                for _n in range(1, n):
                    model.constraints.add(x[(_n, _k)] >= self.u.iloc[_n][_k])
                    model.constraints.add(self.u.iloc[_n][_k] >= self.mu * x[(_n, _k)])

        # buildOptionalConstraints
        if (mode_t in ('FULL', 'CARCH')) and d > 0:
            for _m in range(self.m):
                for _n in range(1, self.n):
                    _sum = 0
                    for _d in range(self.d):
                        _sum += z[(_m, _n, _d)]
                    model.constraints.add(_sum == 1)

            for _m in range(self.m):
                for _M in range(self.M):
                    for _d in range(d):
                        for _i in range(1, self.n-1):
                            for _j in range(1, self.n):
                                model.constraints.add(bitcA[(_M, _m, _i)] - bitcA[(_M, _m, _j)] <= 2 - z[(_m, _i, _d)] - z[(_m, _j, _d)])
                                model.constraints.add(bitcA[(_M, _m, _j)] - bitcA[(_M, _m, _i)] <= 2 - z[(_m, _i, _d)] - z[(_m, _j, _d)])
                                model.constraints.add(bitcB[(_M, _m, _i)] - bitcB[(_M, _m, _j)] <= 2 - z[(_m, _i, _d)] - z[(_m, _j, _d)])
                                model.constraints.add(bitcB[(_M, _m, _j)] - bitcB[(_M, _m, _i)] <= 2 - z[(_m, _i, _d)] - z[(_m, _j, _d)])

            for _m in range(self.m):
                for _d in range(d - 1):
                    _sum_l = _sum_l1 = 0
                    for _n in range(1, self.n):
                        _sum_l += z[(_m, _n, _d)] * self.symmCoeff(_n)
                        _sum_l1 += z[(_m, _n, _d+1)] * self.symmCoeff(_n)
                        model.constraints.add(_sum_l <= _sum_l1)

        if mode_t in ('FULL', 'CARCH'):
            self.build_symmetry_breaking(model)
            self.fix_given_cn(model)

        if mode_t == 'FULL':
            self.hot_start()

        obj = 0
        for _m in range(m):
            for _k in range(k):
                cluster_id = self.cluster_ids[_m]
                obj += (yA[(_m, _k)] + yB[(_m, _k)]) * self.w[cluster_id]

        model.obj = pe.Objective(expr=obj, sense=pe.minimize)
        self.model = model

    def build_symmetry_breaking(self, model):
        for i in range(1, self.n - 1):
            _sum1 = 0
            _sum2 = 0
            for _m in range(self.m):
                _symCoeff = self.symmCoeff(_m)
                _sum1 += self.cA.iloc[_m][i] * _symCoeff + self.cB.iloc[_m][i] * _symCoeff
                _sum2 += self.cA.iloc[_m][i] * _symCoeff + self.cB.iloc[_m][i] * _symCoeff
            model.constraints.add(_sum1 <= _sum2)

    def fix_given_cn(self, model):
        for _m in range(self.m):
            cluster_id = self.f_a.index[_m]
            if cluster_id in self.copy_numbers:
                _cnA, _cnB = self.copy_numbers[cluster_id]
                for _n in range(1, self.n):
                    model.constraints.add(self.cA.iloc[_m][_n] == _cnA)
                    model.constraints.add(self.cB.iloc[_m][_n] == _cnB)

    def first_hot_start(self):

        if self.d > 0:
            targetA = np.empty((self.m, self.d))
            targetB = np.empty_like(targetA)
            for _m in range(self.m):
                targetA[_m, :] = np.random.choice(self.f_a.iloc[_m].values, self.d, replace=False)
                targetB[_m, :] = np.random.choice(self.f_b.iloc[_m].values, self.d, replace=False)
            targetA = targetA.round()
            targetB = targetB.round()
        else:
            targetA = np.round(self.f_a.values)
            targetB = np.round(self.f_b.values)

        hcA = np.zeros((self.m, self.n))
        hcB = np.zeros((self.m, self.n))

        for _m in range(self.m):
            cluster_id = self.f_a.index[_m]

            adA = 0
            adB = 0
            for _n in range(self.n):
                if _n == 0:
                    hcA[_m][0] = 1
                    hcB[_m][0] = 1
                else:
                    if cluster_id in self.copy_numbers:
                        hcA[_m][_n] = self.copy_numbers[cluster_id][0]
                        hcB[_m][_n] = self.copy_numbers[cluster_id][1]
                    else:
                        mod = min(self.n, self.k)
                        a = min(targetA[_m][_n % mod], self.cn_max)
                        b = min(targetB[_m][_n % mod], self.cn_max)

                        if self.ampdel:
                            base = self.base
                            if adA == 0 and a > base:
                                adA = 1
                            if adA == 0 and a < base:
                                adA = -1
                            if adB == 0 and b > base:
                                adB = 1
                            if adB == 0 and b < base:
                                adB = -1

                            a = max(a, base) if adA >= 0 else min(a, base)
                            b = max(b, base) if adB >= 0 else min(b, base)

                            if a + b > self.cn_max:
                                a = self.cn_max - base
                                b = self.cn_max - a

                        if a + b <= self.cn_max:
                            hcA[_m][_n] = a
                            hcB[_m][_n] = b
                        else:
                            hcA[_m][_n] = a
                            hcB[_m][_n] = max(self.cn_max - a, 0)

                        assert hcA[_m][_n] + hcB[_m][_n] <= self.cn_max
                        assert hcA[_m][_n] >= 0
                        assert hcB[_m][_n] >= 0

        cA = pd.DataFrame(index=self.cluster_ids, columns=range(self.n))
        cA[:] = hcA
        cB = pd.DataFrame(index=self.cluster_ids, columns=range(self.n))
        cB[:] = hcB

        return cA, cB

    def hot_start(self, f_a=None, f_b=None):
        def get_rank(hcA, hcB):
            m, n = hcA.shape
            rank = np.zeros(n)
            rank[0] = -1
            for _m in range(m):
                for _n in range(1, n):
                    rank[_n] += hcA.iloc[(_m, _n)] * self.symmCoeff(_m) + hcB.iloc[(_m, _n)] * self.symmCoeff(_m)

            return rank

        if f_a is None and f_b is None:
            f_a, f_b = self.first_hot_start()
        rank = get_rank(f_a, f_b)
        rank_indices = np.argsort(rank)

        for _m in range(self.m):
            for _n in range(0, self.n):
                self.cA.iloc[_m][rank_indices[_n]].start = f_a.iloc[(_m, _n)]
                self.cB.iloc[_m][rank_indices[_n]].start = f_b.iloc[(_m, _n)]

    def fix_u(self, u):
        self._fixed_u[:] = u
        self.mode = 'CARCH'
        # TODO: sanity checks in fixU

    def fix_c(self, cA, cB):
        self._fixed_cA[:] = cA
        self._fixed_cB[:] = cB
        self.mode = 'UARCH'
        # TODO: sanity checks in fixC

    def build_random_u(self, random_seed=None):
        def _calculate_size_bubbles(mu):
            if mu <= 0.1:
                return 10
            elif mu <= 0.15:
                return 6
            elif mu <= 0.2:
                return 5
            else:
                return 3

        def _build_partition_vector(n, n_parts, size_bubbles, mu=0.03):
            """
            Return fractional components of n components as an ndarray of size n
            """
            # initialize all fractions to 0
            result = np.zeros(n)
            # generate n_parts random positions from all available positions
            positions = np.random.choice(np.arange(n), n_parts, replace=False)

            bubbles = np.sort(np.random.choice(np.arange(1, size_bubbles), n_parts - 1, replace=False)) / size_bubbles
            _result = np.diff(bubbles, prepend=0, append=1)

            # set fractions at selected positions to successive difference between bubbles
            result[positions] = _result
            # for fractions that are within tol of mu, clamp them up to mu
            result[(mu - self.tol <= result) & (result < mu)] = mu

            return result

        size_bubbles = _calculate_size_bubbles(self.mu)
        U = np.empty((self.n, self.k))

        with Random(random_seed):
            for _k in range(self.k):
                # Generate 2 random ints between [1, n] and take the max - biasing towards larger numbers
                _n0 = np.random.randint(1, self.n + 1)
                _n1 = np.random.randint(1, self.n + 1)
                n_parts = min(max(_n0, _n1), size_bubbles)  # no. of clones with non-zero proportion in the mix
                v = _build_partition_vector(self.n, n_parts, size_bubbles, mu=self.mu)
                U[:, _k] = v
        return U

    def run(self, solver_type='gurobipy', write_path=None):
        if solver_type == 'gurobipy':
            solver = pe.SolverFactory('gurobi', solver_io='python')
        else:
            solver = pe.SolverFactory(solver_type)

        solver.solve(self.model)
        if write_path is not None:
            self.model.write(write_path)

        # The dataframes cA/cB/u contain Variables with optimized values as the 'value' attribute in each variable
        # We define a function to extract the optimized values out of these Dataframes
        # In cases where we have a frozen value (e.g. when cA/cB/u are fixed), we simply return the value as-is.
        def _value(var): return getattr(var, 'value', var)

        return (
            self.model.obj(),
            self.optimized_cA.applymap(_value).astype(int),
            self.optimized_cB.applymap(_value).astype(int),
            self.optimized_u.applymap(_value)
        )
