import math
import numpy as np
import pandas as pd
from pyomo import environ as pe

from .utils import Random


class ILPSubset:
    def __init__(self, n, cn_max, d, mu, ampdel, copy_numbers, f_a, f_b, w):

        # Each ILPSubset maintains its own data, so make a deep-copy of passed-in DataFrames
        f2_a, f2_b = f_a.copy(deep=True), f_b.copy(deep=True)

        assert f_a.shape == f_b.shape
        assert np.all(f_a.index == f_b.index)
        assert np.all(f_a.columns == f_b.columns)

        self.m, self.k = f_a.shape
        self.f2_a = f2_a
        self.f2_b = f2_b
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
        self.c2A = pd.DataFrame(index=self.cluster_ids, columns=range(n))
        self.c2B = pd.DataFrame(index=self.cluster_ids, columns=range(n))
        self.u2 = pd.DataFrame(index=range(n), columns=self.sample_ids)

        # Fixed values of cA/cB/u
        self._fixed2_cA = pd.DataFrame(index=self.cluster_ids, columns=range(n))
        self._fixed2_cB = pd.DataFrame(index=self.cluster_ids, columns=range(n))
        self._fixed2_u = pd.DataFrame(index=range(n), columns=self.sample_ids)

        self.model2 = None  # initialized on create_model()

    def __copy__(self):
        return ILPSubset(
            n=self.n, cn_max=self.cn_max, d=self.d, mu=self.mu, ampdel=self.ampdel, copy_numbers=self.copy_numbers,
            f_a=self.f2_a, f_b=self.f2_b, w=self.w
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
    def optimized2_cA(self):
        if self.mode == 'UARCH':
            return self._fixed2_cA
        else:
            return self.c2A

    @property
    def optimized2_cB(self):
        if self.mode == 'UARCH':
            return self._fixed2_cB
        else:
            return self.c2B

    @property
    def optimized2_u(self):
        if self.mode == 'CARCH':
            return self._fixed2_u
        else:
            return self.u2

    def create_model(self):

        m, n, k = self.m, self.n, self.k
        f2_a, f2_b = self.f2_a, self.f2_b
        cn_max = self.cn_max
        ampdel = self.ampdel
        copy_numbers = self.copy_numbers
        mode_t = self.mode
        d = self.d
        _M = self.M
        _base = self.base

        model2 = pe.ConcreteModel()

        f2A = {}
        f2B = {}
        y2A = {}
        y2B = {}
        ad2A = {}
        ad2B = {}

        for _m in range(m):

            cluster_id = f2_a.index[_m]

            # upper bound for solver
            ub = max(sum(copy_numbers.get(cluster_id, (0, 0))), cn_max)

            for _k in range(k):
                y2A[(_m, _k)] = pe.Var(bounds=(0, np.inf), domain=pe.Reals)
                model2.add_component(f'yA_{_m + 1}_{_k + 1}', y2A[(_m, _k)])
                y2B[(_m, _k)] = pe.Var(bounds=(0, np.inf), domain=pe.Reals)
                model2.add_component(f'yB_{_m + 1}_{_k + 1}', y2B[(_m, _k)])
                f2A[(_m, _k)] = pe.Var(bounds=(0, ub), domain=pe.Reals)
                model2.add_component(f'fA_{_m + 1}_{_k + 1}', f2A[(_m, _k)])
                f2B[(_m, _k)] = pe.Var(bounds=(0, ub), domain=pe.Reals)
                model2.add_component(f'fB_{_m + 1}_{_k + 1}', f2B[(_m, _k)])

        if mode_t in ('FULL', 'CARCH'):
            for _m in range(m):

                cluster_id = f2_a.index[_m]

                # upper bound for solver
                ub = max(sum(copy_numbers.get(cluster_id, (0, 0))), cn_max)

                for _n in range(n):
                    self.c2A.iloc[_m][_n] = pe.Var(bounds=(0, ub), domain=pe.Integers)
                    model2.add_component(f'cA_{_m + 1}_{_n + 1}', self.c2A.iloc[_m][_n])
                    self.c2B.iloc[_m][_n] = pe.Var(bounds=(0, ub), domain=pe.Integers)
                    model2.add_component(f'cB_{_m + 1}_{_n + 1}', self.c2B.iloc[_m][_n])

            if ampdel:
                for _m in range(m):
                    cluster_id = f2_a.index[_m]
                    if cluster_id not in copy_numbers:
                        ad2A[_m] = pe.Var(bounds=(0, 1), domain=pe.Binary)
                        model2.add_component(f'adA_{_m + 1}', ad2A[_m])
                        ad2B[_m] = pe.Var(bounds=(0, 1), domain=pe.Binary)
                        model2.add_component(f'adB_{_m + 1}', ad2B[_m])

        bitc2A = {}
        bitc2B = {}
        if (mode_t == 'FULL') or (d > 0 and mode_t == 'CARCH'):
            for _b in range(_M):
                for _m in range(m):
                    for _n in range(n):
                        bitc2A[(_b, _m, _n)] = pe.Var(bounds=(0, 1), domain=pe.Binary)
                        model2.add_component(f'bitcA_{_b + 1}_{_m + 1}_{_n + 1}', bitc2A[(_b, _m, _n)])
                        bitc2B[(_b, _m, _n)] = pe.Var(bounds=(0, 1), domain=pe.Binary)
                        model2.add_component(f'bitcB_{_b + 1}_{_m + 1}_{_n + 1}', bitc2B[(_b, _m, _n)])

        if mode_t in ('FULL', 'UARCH'):
            for _n in range(n):
                for _k in range(k):
                    self.u2.iloc[_n][_k] = pe.Var(bounds=(0, 1), domain=pe.Reals)
                    model2.add_component(f'u_{_n + 1}_{_k + 1}', self.u2.iloc[_n][_k])

        v2A = {}
        v2B = {}
        if mode_t == 'FULL':
            for _b in range(_M):
                for _m in range(m):
                    for _n in range(n):
                        for _k in range(k):
                            v2A[(_b, _m, _n, _k)] = pe.Var(bounds=(0, 1), domain=pe.Reals)
                            model2.add_component(f'vA_{_b + 1}_{_m + 1}_{_n + 1}_{_k + 1}', v2A[(_b, _m, _n, _k)])
                            v2B[(_b, _m, _n, _k)] = pe.Var(bounds=(0, 1), domain=pe.Reals)
                            model2.add_component(f'vB_{_b + 1}_{_m + 1}_{_n + 1}_{_k + 1}', v2B[(_b, _m, _n, _k)])

        x2 = {}
        if (mode_t in ('FULL', 'UARCH')) and (self.mu > 0):
            for _n in range(n):
                for _k in range(k):
                    x2[(_n, _k)] = pe.Var(domain=pe.Binary)
                    model2.add_component(f'x_{_n + 1}_{_k + 1}', x2[(_n, _k)])

        # buildOptionalVariables
        z2 = {}
        if (mode_t in ('FULL', 'CARCH')) and d > 0:
            for _m in range(self.m):
                for _n in range(1, self.n):
                    for _d in range(d):
                        z2[(_m, _n, _d)] = pe.Var(bounds=(0, 1), domain=pe.Binary)
                        model2.add_component(f'z_{_m + 1}_{_n + 1}_{_d + 1}', z2[(_m, _n, _d)])

        # CONSTRAINTS
        model2.constraints = pe.ConstraintList()

        for _m in range(m):

            cluster_id = f2_a.index[_m]
            f2_a_values = f2_a.loc[cluster_id].values
            f2_b_values = f2_b.loc[cluster_id].values

            for _k in range(k):
                _yA, _yB = y2A[(_m, _k)], y2B[(_m, _k)]
                _fA, _fB = f2A[(_m, _k)], f2B[(_m, _k)]

                model2.constraints.add(float(f2_a_values[_k]) - _fA <= _yA)
                model2.constraints.add(_fA - float(f2_a_values[_k]) <= _yA)
                model2.constraints.add(float(f2_b_values[_k]) - _fB <= _yB)
                model2.constraints.add(_fB - float(f2_b_values[_k]) <= _yB)

        if mode_t == 'FULL':
            for _m in range(m):
                for _k in range(k):

                    sum2_a = 0
                    for _n in range(n):
                        for _b in range(_M):
                            sum2_a += v2A[(_b, _m, _n, _k)] * math.pow(2, _b)
                            model2.constraints.add(v2A[(_b, _m, _n, _k)] <= bitc2A[(_b, _m, _n)])
                            model2.constraints.add(v2A[(_b, _m, _n, _k)] <= self.u2.iloc[_n][_k])
                            model2.constraints.add(v2A[(_b, _m, _n, _k)] >= bitc2A[(_b, _m, _n)] + self.u2.iloc[_n][_k] - 1)

                    model2.constraints.add(f2A[(_m, _k)] == sum2_a)

            for _m in range(m):
                for _k in range(k):

                    sum2_b = 0
                    for _n in range(n):
                        for _b in range(_M):
                            sum2_b += v2B[(_b, _m, _n, _k)] * math.pow(2, _b)
                            model2.constraints.add(v2B[(_b, _m, _n, _k)] <= bitc2B[(_b, _m, _n)])
                            model2.constraints.add(v2B[(_b, _m, _n, _k)] <= self.u2.iloc[_n][_k])
                            model2.constraints.add(v2B[(_b, _m, _n, _k)] >= bitc2B[(_b, _m, _n)] + self.u2.iloc[_n][_k] - 1)

                    model2.constraints.add(f2B[(_m, _k)] == sum2_b)

            for _n in range(n):
                for _k in range(k):
                    _sum2 = 0
                    for _m in range(m):
                        for _b in range(_M):
                            _sum2 += bitc2A[(_b, _m, _n)] + bitc2B[(_b, _m, _n)]
                    model2.constraints.add(_sum2 >= self.u2.iloc[_n][_k])

        if (mode_t == 'FULL') or (d > 0 and mode_t == 'CARCH'):
            for _m in range(m):
                cluster_id = f2_a.index[_m]
                # upper bound for solver
                ub = max(sum(copy_numbers.get(cluster_id, (0, 0))), cn_max)

                for _n in range(n):
                    sum2_a = 0
                    sum2_b = 0
                    for _b in range(_M):
                        sum2_a += bitc2A[(_b, _m, _n)] * math.pow(2, _b)
                        sum2_b += bitc2B[(_b, _m, _n)] * math.pow(2, _b)

                    model2.constraints.add(self.c2A.iloc[_m][_n] == sum2_a)
                    model2.constraints.add(self.c2B.iloc[_m][_n] == sum2_b)
                    model2.constraints.add(self.c2A.iloc[_m][_n] + self.c2B.iloc[_m][_n] <= ub)

        if mode_t == 'CARCH':
            # TODO: These loops can be collapsed once validation against C++ is complete
            for _m in range(m):
                for _k in range(k):
                    _sum2 = 0
                    for _n in range(n):
                        if self._fixed2_u.iloc[_n][_k] >= self.mu - self.tol:
                            _sum2 += self.c2A.iloc[_m][_n] * self._fixed2_u.iloc[_n][_k]
                    model2.constraints.add(f2A[(_m, _k)] == _sum2)

            for _m in range(m):
                for _k in range(k):
                    _sum2 = 0
                    for _n in range(n):
                        if self._fixed2_u.iloc[_n][_k] >= self.mu - self.tol:
                            _sum2 += self.c2B.iloc[_m][_n] * self._fixed2_u.iloc[_n][_k]
                    model2.constraints.add(f2B[(_m, _k)] == _sum2)

            for _m in range(m):
                cluster_id = f2_a.index[_m]
                # upper bound for solver
                ub = max(sum(copy_numbers.get(cluster_id, (0, 0))), cn_max)
                for _n in range(n):
                    model2.constraints.add(self.c2A.iloc[_m][_n] + self.c2B.iloc[_m][_n] <= ub)

        if mode_t in ('FULL', 'CARCH'):
            for _m in range(m):
                model2.constraints.add(self.c2A.iloc[_m][0] == 1)
                model2.constraints.add(self.c2B.iloc[_m][0] == 1)

            if ampdel:
                for _m in range(m):
                    cluster_id = f2_a.index[_m]
                    if cluster_id not in copy_numbers:
                        for _n in range(1, n):
                            model2.constraints.add(self.c2A.iloc[_m][_n] <= cn_max * ad2A[_m] + _base - _base * ad2A[_m])
                            model2.constraints.add(self.c2A.iloc[_m][_n] >= _base * ad2A[_m])
                            model2.constraints.add(self.c2B.iloc[_m][_n] <= cn_max * ad2B[_m] + _base - _base * ad2B[_m])
                            model2.constraints.add(self.c2B.iloc[_m][_n] >= _base * ad2B[_m])

        if mode_t == 'UARCH':
            # TODO: These loops can be collapsed once validation against C++ is complete
            for _m in range(m):
                for _k in range(k):
                    _sum2 = 0
                    for _n in range(n):
                        _sum2 += int(self._fixed2_cA.iloc[_m][_n]) * self.u2.iloc[_n][_k]
                    model2.constraints.add(f2A[(_m, _k)] == _sum2)

            for _m in range(m):
                for _k in range(k):
                    _sum2 = 0
                    for _n in range(n):
                        _sum2 += int(self._fixed2_cB.iloc[_m][_n]) * self.u2.iloc[_n][_k]
                    model2.constraints.add(f2B[(_m, _k)] == _sum2)

        if mode_t in ('FULL', 'UARCH'):
            for _k in range(k):
                _sum2 = 0
                for _n in range(n):
                    _sum2 += self.u2.iloc[_n][_k]
                model2.constraints.add(_sum2 == 1)

        if (mode_t in ('FULL', 'UARCH')) and self.mu > 0:
            for _k in range(k):
                for _n in range(1, n):
                    model2.constraints.add(x2[(_n, _k)] >= self.u2.iloc[_n][_k])
                    model2.constraints.add(self.u2.iloc[_n][_k] >= self.mu * x2[(_n, _k)])

        # buildOptionalConstraints
        if (mode_t in ('FULL', 'CARCH')) and d > 0:
            for _m in range(self.m):
                for _n in range(1, self.n):
                    _sum2 = 0
                    for _d in range(self.d):
                        _sum2 += z2[(_m, _n, _d)]
                    model2.constraints.add(_sum2 == 1)

            for _m in range(self.m):
                for _M in range(self.M):
                    for _d in range(d):
                        for _i in range(1, self.n-1):
                            for _j in range(1, self.n):
                                model2.constraints.add(bitc2A[(_M, _m, _i)] - bitc2A[(_M, _m, _j)] <= 2 - z2[(_m, _i, _d)] - z2[(_m, _j, _d)])
                                model2.constraints.add(bitc2A[(_M, _m, _j)] - bitc2A[(_M, _m, _i)] <= 2 - z2[(_m, _i, _d)] - z2[(_m, _j, _d)])
                                model2.constraints.add(bitc2B[(_M, _m, _i)] - bitc2B[(_M, _m, _j)] <= 2 - z2[(_m, _i, _d)] - z2[(_m, _j, _d)])
                                model2.constraints.add(bitc2B[(_M, _m, _j)] - bitc2B[(_M, _m, _i)] <= 2 - z2[(_m, _i, _d)] - z2[(_m, _j, _d)])

            for _m in range(self.m):
                for _d in range(d - 1):
                    _sum2_l = _sum2_l1 = 0
                    for _n in range(1, self.n):
                        _sum2_l += z2[(_m, _n, _d)] * self.symmCoeff(_n)
                        _sum2_l1 += z2[(_m, _n, _d+1)] * self.symmCoeff(_n)
                        model2.constraints.add(_sum2_l <= _sum2_l1)

        if mode_t in ('FULL', 'CARCH'):
            self.build_symmetry_breaking(model2)
            self.fix_given_cn(model2)

        if mode_t == 'FULL':
            self.hot_start()

        obj = 0
        for _m in range(m):
            for _k in range(k):
                cluster_id = self.cluster_ids[_m]
                # obj += (y2A[(_m, _k)] + y2B[(_m, _k)]) * self.w[cluster_id]
                obj += float(self.w[cluster_id]) * y2A[(_m, _k)] + float(self.w[cluster_id]) * y2B[(_m, _k)]

        model2.obj = pe.Objective(expr=obj, sense=pe.minimize)
        self.model2 = model2

    def build_symmetry_breaking(self, model2):
        for i in range(1, self.n - 1):
            _sum1 = 0
            _sum2 = 0
            for _m in range(self.m):
                _symCoeff = self.symmCoeff(_m)
                _sum1 += self.c2A.iloc[_m][i] * _symCoeff + self.c2B.iloc[_m][i] * _symCoeff
                _sum2 += self.c2A.iloc[_m][i] * _symCoeff + self.c2B.iloc[_m][i] * _symCoeff
            model2.constraints.add(_sum1 <= _sum2)

    def fix_given_cn(self, model2):
        for _m in range(self.m):
            cluster_id = self.f2_a.index[_m]
            if cluster_id in self.copy_numbers:
                _cnA, _cnB = self.copy_numbers[cluster_id]
                for _n in range(1, self.n):
                    model2.constraints.add(self.c2A.iloc[_m][_n] == _cnA)
                    model2.constraints.add(self.c2B.iloc[_m][_n] == _cnB)

    def first_hot_start(self):

        if self.d > 0:
            targetA = np.empty((self.m, self.d))
            targetB = np.empty_like(targetA)
            for _m in range(self.m):
                targetA[_m, :] = np.random.choice(self.f2_a.iloc[_m].values, self.d, replace=False)
                targetB[_m, :] = np.random.choice(self.f2_b.iloc[_m].values, self.d, replace=False)
            targetA = targetA.round()
            targetB = targetB.round()
        else:
            targetA = np.round(self.f2_a.values)
            targetB = np.round(self.f2_b.values)

        hcA = np.zeros((self.m, self.n))
        hcB = np.zeros((self.m, self.n))

        for _m in range(self.m):
            cluster_id = self.f2_a.index[_m]

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
                self.c2A.iloc[_m][rank_indices[_n]].start = f_a.iloc[(_m, _n)]
                self.c2B.iloc[_m][rank_indices[_n]].start = f_b.iloc[(_m, _n)]

    def fix_u(self, u):
        self._fixed2_u[:] = u
        self.mode = 'CARCH'
        # TODO: sanity checks in fixU

    def fix_c(self, cA, cB):
        self._fixed2_cA[:] = cA
        self._fixed2_cB[:] = cB
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

    def run(self, write_path=None):
        solver = pe.SolverFactory('gurobi', solver_io='python')
        solver.solve(self.model2)
        if write_path is not None:
            self.model2.write(write_path)

        # The dataframes cA/cB/u contain Variables with optimized values as the 'value' attribute in each variable
        # We define a function to extract the optimized values out of these Dataframes
        # In cases where we have a frozen value (e.g. when cA/cB/u are fixed), we simply return the value as-is.
        def _value(var): return getattr(var, 'value', var)

        return (
            self.model2.obj(),
            self.optimized2_cA.applymap(_value).astype(int),
            self.optimized2_cB.applymap(_value).astype(int),
            self.optimized2_u.applymap(_value)
        )
