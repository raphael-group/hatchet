import textwrap
import math
import numpy as np
from pyomo import environ as pe
from pyomo.opt import SolverStatus, TerminationCondition

from hatchet.utils.solve.utils import Random


class ILPSubset:
    def __init__(self, n, cn_max, d, mu, ampdel, copy_numbers,
                 f_a, f_b, w, clus_adj_counts, maes, evolcons, bp_max, uniqueclones, purities):

        # Each ILPSubset maintains its own data, so make a deep-copy of passed-in DataFrames
        f_a, f_b = f_a.copy(deep=True), f_b.copy(deep=True)

        assert f_a.shape == f_b.shape
        assert np.all(f_a.index == f_b.index)
        assert np.all(f_a.columns == f_b.columns)

        # make sure that minimum absolute error is all positive
        assert np.all(maes > 0)

        self.m, self.k = f_a.shape
        self.uniqueclones = uniqueclones

        if self.uniqueclones:
            # each sample has its own clone + there is a normal clone + there may be some extra clones
            assert n >= self.k + 1

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
        self.clus_adj_counts = clus_adj_counts
        self.maes = maes
        self.evolcons = evolcons
        self.bp_max = bp_max
        self.purities = purities

        self.tol = 0.001

        self.mode = 'FULL'

        # Values we want to optimize for, as dataframes
        self.cA = [[np.nan for _ in range(self.n)] for _ in range(self.m)]
        self.cB = [[np.nan for _ in range(self.n)] for _ in range(self.m)]
        self.u = [[np.nan for _ in range(self.k)] for _ in range(self.n)]

        # Fixed values of cA/cB/u
        self._fixed_cA = [[np.nan for _ in range(n)] for _ in range(self.m)]
        self._fixed_cB = [[np.nan for _ in range(n)] for _ in range(self.m)]
        self._fixed_u = [[np.nan for _ in range(self.k)] for _ in range(self.n)]

        self.warmstart = False  # set on hot_start()
        self.model = None  # initialized on create_model()

    def __copy__(self):
        return ILPSubset(
            n=self.n,
            cn_max=self.cn_max,
            d=self.d,
            mu=self.mu,
            ampdel=self.ampdel,
            copy_numbers=self.copy_numbers,
            f_a=self.f_a,
            f_b=self.f_b,
            w=self.w,
            clus_adj_counts=self.clus_adj_counts,
            maes=self.maes,
            evolcons=self.evolcons,
            bp_max=self.bp_max,
            uniqueclones=self.uniqueclones,
            purities=self.purities,
        )

    def __str__(self):
        # Pyomo pprint gives us too much information - too unwieldy for large models
        # This method is implemented to supply the bare-minimum but useful model information.
        if self.model is None:
            return ''
        else:
            return textwrap.dedent(
                f"""
                # ------------------------------------------
                #   Problem Information
                # ------------------------------------------
                #     Number of constraints: {self.model.nconstraints()}
                #     Number of variables: {self.model.nvariables()}
                # ------------------------------------------
            """
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

    def create_model(self, pprint=False):

        m, n, k = self.m, self.n, self.k
        f_a, f_b = self.f_a, self.f_b
        cn_max = self.cn_max
        ampdel = self.ampdel
        copy_numbers = self.copy_numbers
        mode_t = self.mode
        d = self.d
        _M = self.M
        _base = self.base
        maes = self.maes
        purities = self.purities

        model = pe.ConcreteModel()
        # CONSTRAINTS
        model.constraints = pe.ConstraintList()

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
                    self.cA[_m][_n] = pe.Var(bounds=(0, ub), domain=pe.Integers)
                    model.add_component(f'cA_{_m + 1}_{_n + 1}', self.cA[_m][_n])
                    self.cB[_m][_n] = pe.Var(bounds=(0, ub), domain=pe.Integers)
                    model.add_component(f'cB_{_m + 1}_{_n + 1}', self.cB[_m][_n])

            if ampdel:
                for _m in range(m):
                    cluster_id = f_a.index[_m]
                    if cluster_id not in copy_numbers:
                        adA[_m] = pe.Var(bounds=(0, 1), domain=pe.Binary)
                        model.add_component(f'adA_{_m + 1}', adA[_m])
                        adB[_m] = pe.Var(bounds=(0, 1), domain=pe.Binary)
                        model.add_component(f'adB_{_m + 1}', adB[_m])

        # adding evolution constraints
        if self.evolcons:
            self._add_evol_cons(model, ub)

        bitcA = {}
        bitcB = {}
        if (mode_t == 'FULL') or (d > 0 and mode_t == 'CARCH'):
            for _b in range(_M):
                for _m in range(m):
                    for _n in range(n):
                        bitcA[(_b, _m, _n)] = pe.Var(bounds=(0, 1), domain=pe.Binary)
                        model.add_component(
                            f'bitcA_{_b + 1}_{_m + 1}_{_n + 1}',
                            bitcA[(_b, _m, _n)],
                        )
                        bitcB[(_b, _m, _n)] = pe.Var(bounds=(0, 1), domain=pe.Binary)
                        model.add_component(
                            f'bitcB_{_b + 1}_{_m + 1}_{_n + 1}',
                            bitcB[(_b, _m, _n)],
                        )

        if mode_t in ('FULL', 'UARCH'):
            for _n in range(n):
                for _k in range(k):
                    self.u[_n][_k] = pe.Var(bounds=(0, 1), domain=pe.Reals)
                    model.add_component(f'u_{_n + 1}_{_k + 1}', self.u[_n][_k])

        vA = {}
        vB = {}
        if mode_t == 'FULL':
            for _b in range(_M):
                for _m in range(m):
                    for _n in range(n):
                        for _k in range(k):
                            vA[(_b, _m, _n, _k)] = pe.Var(bounds=(0, 1), domain=pe.Reals)
                            model.add_component(
                                f'vA_{_b + 1}_{_m + 1}_{_n + 1}_{_k + 1}',
                                vA[(_b, _m, _n, _k)],
                            )
                            vB[(_b, _m, _n, _k)] = pe.Var(bounds=(0, 1), domain=pe.Reals)
                            model.add_component(
                                f'vB_{_b + 1}_{_m + 1}_{_n + 1}_{_k + 1}',
                                vB[(_b, _m, _n, _k)],
                            )

        # self.mu = 0 always (future work)
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
                            model.constraints.add(vA[(_b, _m, _n, _k)] <= self.u[_n][_k])
                            model.constraints.add(vA[(_b, _m, _n, _k)] >= bitcA[(_b, _m, _n)] + self.u[_n][_k] - 1)

                    model.constraints.add(fA[(_m, _k)] == sum_a)

            for _m in range(m):
                for _k in range(k):

                    sum_b = 0
                    for _n in range(n):
                        for _b in range(_M):
                            sum_b += vB[(_b, _m, _n, _k)] * math.pow(2, _b)
                            model.constraints.add(vB[(_b, _m, _n, _k)] <= bitcB[(_b, _m, _n)])
                            model.constraints.add(vB[(_b, _m, _n, _k)] <= self.u[_n][_k])
                            model.constraints.add(vB[(_b, _m, _n, _k)] >= bitcB[(_b, _m, _n)] + self.u[_n][_k] - 1)

                    model.constraints.add(fB[(_m, _k)] == sum_b)

            for _n in range(n):
                for _k in range(k):
                    _sum = 0
                    for _m in range(m):
                        for _b in range(_M):
                            _sum += bitcA[(_b, _m, _n)] + bitcB[(_b, _m, _n)]
                    model.constraints.add(_sum >= self.u[_n][_k])

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

                    model.constraints.add(self.cA[_m][_n] == sum_a)
                    model.constraints.add(self.cB[_m][_n] == sum_b)
                    model.constraints.add(self.cA[_m][_n] + self.cB[_m][_n] <= ub)

        #  this is where we fix the purities if user provides purity values for (tumor) samples
        if mode_t in ('FULL', 'UARCH') and purities:
            lbs = purities[0]
            ubs = purities[1]
            for i, p in enumerate(lbs):
                model.constraints.add(self.u[0][i] <= 1 - p)
            for i, p in enumerate(ubs):
                model.constraints.add(self.u[0][i] >= 1 - p)

        if mode_t == 'CARCH':
            # TODO: These loops can be collapsed once validation against C++ is complete
            for _m in range(m):
                for _k in range(k):
                    _sumA = 0
                    _sumB = 0
                    for _n in range(n):
                        if self._fixed_u[_n][_k] >= self.mu - self.tol:
                            _sumA += self.cA[_m][_n] * self._fixed_u[_n][_k]
                            _sumB += self.cB[_m][_n] * self._fixed_u[_n][_k]
                    model.constraints.add(fA[(_m, _k)] == _sumA)
                    model.constraints.add(fB[(_m, _k)] == _sumB)

                cluster_id = f_a.index[_m]
                # upper bound for solver
                ub = max(sum(copy_numbers.get(cluster_id, (0, 0))), cn_max)
                for _n in range(n):
                    model.constraints.add(self.cA[_m][_n] + self.cB[_m][_n] <= ub)

        if mode_t in ('FULL', 'CARCH'):
            for _m in range(m):
                model.constraints.add(self.cA[_m][0] == 1)
                model.constraints.add(self.cB[_m][0] == 1)
                cluster_id = self.cluster_ids[_m]

                # if the fraction of the cluster size is larger than or equal to 0.02 in the whole dataset
                # then its integer copy number cannot have (0,0) copy number
                if self.w[cluster_id] / sum(self.w) >= 0.02:
                    for _n in range(1, n):
                        model.constraints.add(self.cA[_m][_n] + self.cB[_m][_n] >= 1)

            if ampdel:
                for _m in range(m):
                    cluster_id = f_a.index[_m]
                    if cluster_id not in copy_numbers:
                        for _n in range(1, n):
                            model.constraints.add(self.cA[_m][_n] <= cn_max * adA[_m] + _base - _base * adA[_m])
                            model.constraints.add(self.cA[_m][_n] >= _base * adA[_m])
                            model.constraints.add(self.cB[_m][_n] <= cn_max * adB[_m] + _base - _base * adB[_m])
                            model.constraints.add(self.cB[_m][_n] >= _base * adB[_m])

        if mode_t == 'UARCH':
            # TODO: These loops can be collapsed once validation against C++ is complete
            for _m in range(m):
                for _k in range(k):
                    _sumA = 0
                    _sumB = 0
                    for _n in range(n):
                        _sumA += int(self._fixed_cA[_m][_n]) * self.u[_n][_k]
                        _sumB += int(self._fixed_cB[_m][_n]) * self.u[_n][_k]
                    model.constraints.add(fA[(_m, _k)] == _sumA)
                    model.constraints.add(fB[(_m, _k)] == _sumB)

        if mode_t in ('FULL', 'UARCH'):
            for _k in range(k):
                _sum = 0
                for _n in range(n):
                    _sum += self.u[_n][_k]
                model.constraints.add(_sum == 1)

            if self.uniqueclones:
                for _k in range(k):
                    for _n in range(1, k + 1):
                        if _n != _k + 1:
                            model.constraints.add(self.u[_n][_k] == 0)

        # self.mu = 0 always true (future)
        if (mode_t in ('FULL', 'UARCH')) and self.mu > 0:
            for _k in range(k):
                for _n in range(1, n):
                    model.constraints.add(x[(_n, _k)] >= self.u[_n][_k])
                    model.constraints.add(self.u[_n][_k] >= self.mu * x[(_n, _k)])

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
                        for _i in range(1, self.n - 1):
                            for _j in range(1, self.n):
                                model.constraints.add(
                                    bitcA[(_M, _m, _i)] - bitcA[(_M, _m, _j)] <= 2 - z[(_m, _i, _d)] - z[(_m, _j, _d)]
                                )
                                model.constraints.add(
                                    bitcA[(_M, _m, _j)] - bitcA[(_M, _m, _i)] <= 2 - z[(_m, _i, _d)] - z[(_m, _j, _d)]
                                )
                                model.constraints.add(
                                    bitcB[(_M, _m, _i)] - bitcB[(_M, _m, _j)] <= 2 - z[(_m, _i, _d)] - z[(_m, _j, _d)]
                                )
                                model.constraints.add(
                                    bitcB[(_M, _m, _j)] - bitcB[(_M, _m, _i)] <= 2 - z[(_m, _i, _d)] - z[(_m, _j, _d)]
                                )

            for _m in range(self.m):
                for _d in range(d - 1):
                    _sum_l = _sum_l1 = 0
                    for _n in range(1, self.n):
                        _sum_l += z[(_m, _n, _d)] * self.symmCoeff(_n)
                        _sum_l1 += z[(_m, _n, _d + 1)] * self.symmCoeff(_n)
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
                sample_id = self.sample_ids[_k]
                obj += (yA[(_m, _k)] + yB[(_m, _k)]) * self.w[cluster_id] / maes.loc[cluster_id, sample_id]

        model.obj = pe.Objective(expr=obj, sense=pe.minimize)
        self.model = model

        if pprint:
            print(str(self))

    def _add_evol_cons(self, model, ub):
        n = self.n
        f_a = self.f_a
        mode_t = self.mode
        w_max = self.bp_max
        wA = {}
        wB = {}
        sAB = {}

        if mode_t in ('FULL', 'CARCH'):

            for _n1 in range(1, n):
                for _n2 in range(1, n):
                    if _n1 < _n2:
                        sAB[(_n1, _n2)] = pe.Var(bounds=(0, 2 * w_max), domain=pe.Reals)
                        model.add_component(f'sAB_{_n1 + 1}_{_n2 + 1}', sAB[(_n1, _n2)])
                        for _c1, _c2 in self.clus_adj_counts.keys():
                            wA[(_c1, _c2, _n1, _n2)] = pe.Var(bounds=(0, 2 * ub), domain=pe.Reals)
                            model.add_component(f'wA_{_c1}_{_c2}_{_n1 + 1}_{_n2 + 1}', wA[(_c1, _c2, _n1, _n2)])
                            wB[(_c1, _c2, _n1, _n2)] = pe.Var(bounds=(0, 2 * ub), domain=pe.Reals)
                            model.add_component(f'wB_{_c1}_{_c2}_{_n1 + 1}_{_n2 + 1}', wB[(_c1, _c2, _n1, _n2)])

            # breakpoint distance upperbound constraint.
            # skip the clone 1 because it's the normal clone.
            clust_to_ind_map = dict([(j, i) for i, j in list(enumerate(f_a.index))])
            # cluster id to f_a & f_b index map
            for _n1 in range(1, n):
                for _n2 in range(1, n):
                    if _n1 < _n2:
                        sum_w = 0
                        # cluster ids including the telomere "-1"
                        for adj, count in self.clus_adj_counts.items():
                            _ci, _cj = adj  # these are cluster ids
                            _wA, _wB = wA[(_ci, _cj, _n1, _n2)], wB[(_ci, _cj, _n1, _n2)]
                            sum_w += count * (_wA + _wB)
                            if _ci == -1:   # a telomere
                                _cA1i, _cB1i = 0, 0
                                _cA2i, _cB2i = 0, 0
                            else:
                                # this is cluster id index in f_a and cA matrix which is sorted alphabetically
                                _cxi = clust_to_ind_map[_ci]
                                _cA1i, _cB1i = self.cA[_cxi][_n1], self.cB[_cxi][_n1]
                                _cA2i, _cB2i = self.cA[_cxi][_n2], self.cB[_cxi][_n2]
                            # _cj cannot be -1 since during construction of clus_adj_counts _ci < _cj is warranted.
                            _cxj = clust_to_ind_map[_cj]
                            _cA1j, _cB1j = self.cA[_cxj][_n1], self.cB[_cxj][_n1]
                            _cA2j, _cB2j = self.cA[_cxj][_n2], self.cB[_cxj][_n2]
                            model.constraints.add((_cA1j - _cA1i) - (_cA2j - _cA2i) <= _wA)
                            model.constraints.add((_cA2j - _cA2i) - (_cA1j - _cA1i) <= _wA)
                            model.constraints.add((_cB1j - _cB1i) - (_cB2j - _cB2i) <= _wB)
                            model.constraints.add((_cB2j - _cB2i) - (_cB1j - _cB1i) <= _wB)
                        # for quickly accessing pairwise breakpoint distances between clones
                        model.constraints.add(sAB[(_n1, _n2)] == sum_w)
                        model.constraints.add(sAB[(_n1, _n2)] <= 2 * w_max)

    def build_symmetry_breaking(self, model):
        for i in range(1, self.n - 1):
            _sum1 = 0
            _sum2 = 0
            for _m in range(self.m):
                _symCoeff = self.symmCoeff(_m)
                _sum1 += self.cA[_m][i] * _symCoeff + self.cB[_m][i] * _symCoeff
                _sum2 += self.cA[_m][i] * _symCoeff + self.cB[_m][i] * _symCoeff
            model.constraints.add(_sum1 <= _sum2)

    def fix_given_cn(self, model):
        for _m in range(self.m):
            cluster_id = self.f_a.index[_m]
            if cluster_id in self.copy_numbers:
                _cnA, _cnB = self.copy_numbers[cluster_id]
                for _n in range(1, self.n):
                    model.constraints.add(self.cA[_m][_n] == _cnA)
                    model.constraints.add(self.cB[_m][_n] == _cnB)

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

        return hcA, hcB

    def hot_start(self, f_a=None, f_b=None):
        def get_rank(hcA, hcB):
            # m, n = hcA.shape
            m, n = len(hcA), len(hcA[0])
            rank = np.zeros(n)
            rank[0] = -1
            for _m in range(m):
                for _n in range(1, n):
                    rank[_n] += hcA[_m][_n] * self.symmCoeff(_m) + hcB[_m][_n] * self.symmCoeff(_m)

            return rank

        if f_a is None and f_b is None:
            f_a, f_b = self.first_hot_start()
        rank = get_rank(f_a, f_b)
        rank_indices = np.argsort(rank)

        for _m in range(self.m):
            for _n in range(0, self.n):
                self.cA[_m][rank_indices[_n]].value = f_a[_m][_n]
                self.cB[_m][rank_indices[_n]].value = f_b[_m][_n]
        self.warmstart = True

    def fix_u(self, u):
        self._fixed_u[:] = u
        self.mode = 'CARCH'
        # TODO: sanity checks in fixU

    def fix_c(self, cA, cB):
        self._fixed_cA = cA
        self._fixed_cB = cB
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

    def run(self, solver_type='gurobi', timelimit=None, write_path=None):
        if solver_type == 'gurobipy':
            solver = pe.SolverFactory('gurobi', solver_io='python')
        else:
            solver = pe.SolverFactory(solver_type)

        kwargs = {'report_timing': False}
        if timelimit is not None:
            kwargs['timelimit'] = int(timelimit)
        if solver.warm_start_capable():
            kwargs['warmstart'] = self.warmstart

        results = solver.solve(self.model, **kwargs)
        if (
            (results.solver.status == SolverStatus.ok)
            and (
                results.solver.termination_condition
                in (
                    TerminationCondition.optimal,
                    TerminationCondition.feasible,
                )
            )
            or (
                results.solver.status == SolverStatus.aborted
                and results.solver.termination_condition == TerminationCondition.maxTimeLimit
            )
        ):
            pass
        else:
            return None

        if write_path is not None:
            self.model.write(write_path)

        return (
            self.model.obj(),
            [[int(getattr(x, 'value', x)) for x in row] for row in self.optimized_cA],
            [[int(getattr(x, 'value', x)) for x in row] for row in self.optimized_cB],
            [[getattr(x, 'value', x) for x in row] for row in self.optimized_u],
            self.cluster_ids,
            self.sample_ids,
        )


class ILPSubsetSplit(ILPSubset):
    def __init__(
        self,
        n,
        cn_max,
        d,
        mu,
        ampdel,
        copy_numbers,
        f_a,
        f_b,
        binsA,
        binsB,
        lengths,
    ):

        # Each ILPSubset maintains its own data, so make a deep-copy of passed-in DataFrames
        f_a, f_b = f_a.copy(deep=True), f_b.copy(deep=True)
        self.binsA = {k: v.copy(deep=True) for k, v in binsA.items()}
        self.binsB = {k: v.copy(deep=True) for k, v in binsB.items()}
        self.lengths = {k: np.copy(v) for k, v in lengths.items()}

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

        self.tol = 0.001

        self.mode = 'FULL'

        # Values we want to optimize for, as dataframes
        self.cA = [[np.nan for _ in range(self.n)] for _ in range(self.m)]
        self.cB = [[np.nan for _ in range(self.n)] for _ in range(self.m)]
        self.u = [[np.nan for _ in range(self.k)] for _ in range(self.n)]

        # Fixed values of cA/cB/u
        self._fixed_cA = [[np.nan for _ in range(n)] for _ in range(self.m)]
        self._fixed_cB = [[np.nan for _ in range(n)] for _ in range(self.m)]
        self._fixed_u = [[np.nan for _ in range(self.k)] for _ in range(self.n)]

        self.warmstart = False  # set on hot_start()
        self.model = None  # initialized on create_model()

    def __copy__(self):
        return ILPSubsetSplit(
            n=self.n,
            cn_max=self.cn_max,
            d=self.d,
            mu=self.mu,
            ampdel=self.ampdel,
            copy_numbers=self.copy_numbers,
            f_a=self.f_a,
            f_b=self.f_b,
            binsA=self.binsA,
            binsB=self.binsB,
            lengths=self.lengths,
        )

    def create_model(self, pprint=False):

        m, n, k = self.m, self.n, self.k
        f_a = self.f_a
        cn_max = self.cn_max
        ampdel = self.ampdel
        copy_numbers = self.copy_numbers
        mode_t = self.mode
        d = self.d
        _M = self.M
        _base = self.base
        maes = self.maes

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

                fA[(_m, _k)] = pe.Var(bounds=(0, ub), domain=pe.Reals)
                model.add_component(f'fA_{_m + 1}_{_k + 1}', fA[(_m, _k)])
                fB[(_m, _k)] = pe.Var(bounds=(0, ub), domain=pe.Reals)
                model.add_component(f'fB_{_m + 1}_{_k + 1}', fB[(_m, _k)])

                # Need an objective function term for each bin for absolute vlaue
                for _i in range(len(self.binsA[cluster_id])):
                    yA[(_m, _k, _i)] = pe.Var(bounds=(0, np.inf), domain=pe.Reals)
                    model.add_component(f'yA_{_m + 1}_{_k + 1}_{_i + 1}', yA[(_m, _k, _i)])
                    yB[(_m, _k, _i)] = pe.Var(bounds=(0, np.inf), domain=pe.Reals)
                    model.add_component(f'yB_{_m + 1}_{_k + 1}_{_i + 1}', yB[(_m, _k, _i)])

        if mode_t in ('FULL', 'CARCH'):
            for _m in range(m):

                cluster_id = f_a.index[_m]

                # upper bound for solver
                ub = max(sum(copy_numbers.get(cluster_id, (0, 0))), cn_max)

                for _n in range(n):
                    self.cA[_m][_n] = pe.Var(bounds=(0, ub), domain=pe.Integers)
                    model.add_component(f'cA_{_m + 1}_{_n + 1}', self.cA[_m][_n])
                    self.cB[_m][_n] = pe.Var(bounds=(0, ub), domain=pe.Integers)
                    model.add_component(f'cB_{_m + 1}_{_n + 1}', self.cB[_m][_n])

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
                        model.add_component(
                            f'bitcA_{_b + 1}_{_m + 1}_{_n + 1}',
                            bitcA[(_b, _m, _n)],
                        )
                        bitcB[(_b, _m, _n)] = pe.Var(bounds=(0, 1), domain=pe.Binary)
                        model.add_component(
                            f'bitcB_{_b + 1}_{_m + 1}_{_n + 1}',
                            bitcB[(_b, _m, _n)],
                        )

        if mode_t in ('FULL', 'UARCH'):
            for _n in range(n):
                for _k in range(k):
                    self.u[_n][_k] = pe.Var(bounds=(0, 1), domain=pe.Reals)
                    model.add_component(f'u_{_n + 1}_{_k + 1}', self.u[_n][_k])

        vA = {}
        vB = {}
        if mode_t == 'FULL':
            for _b in range(_M):
                for _m in range(m):
                    for _n in range(n):
                        for _k in range(k):
                            vA[(_b, _m, _n, _k)] = pe.Var(bounds=(0, 1), domain=pe.Reals)
                            model.add_component(
                                f'vA_{_b + 1}_{_m + 1}_{_n + 1}_{_k + 1}',
                                vA[(_b, _m, _n, _k)],
                            )
                            vB[(_b, _m, _n, _k)] = pe.Var(bounds=(0, 1), domain=pe.Reals)
                            model.add_component(
                                f'vB_{_b + 1}_{_m + 1}_{_n + 1}_{_k + 1}',
                                vB[(_b, _m, _n, _k)],
                            )

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
            f_a_values = self.binsA[cluster_id].values
            f_b_values = self.binsB[cluster_id].values

            for _k in range(k):
                _fA, _fB = fA[(_m, _k)], fB[(_m, _k)]

                for _i in range(f_a_values.shape[0]):
                    _yA, _yB = yA[(_m, _k, _i)], yB[(_m, _k, _i)]

                    # Instead of 1 term per inferred CN, use all values in bin
                    model.constraints.add(f_a_values[_i, _k] - _fA <= _yA)
                    model.constraints.add(_fA - f_a_values[_i, _k] <= _yA)

                    # model.constraints.add(myA >= _yA)
                    model.constraints.add(f_b_values[_i, _k] - _fB <= _yB)
                    model.constraints.add(_fB - f_b_values[_i, _k] <= _yB)
                    # model.constraints.add(myB >= _yB)

        if mode_t == 'FULL':
            for _m in range(m):
                for _k in range(k):

                    sum_a = 0
                    for _n in range(n):
                        for _b in range(_M):
                            sum_a += vA[(_b, _m, _n, _k)] * math.pow(2, _b)
                            model.constraints.add(vA[(_b, _m, _n, _k)] <= bitcA[(_b, _m, _n)])
                            model.constraints.add(vA[(_b, _m, _n, _k)] <= self.u[_n][_k])
                            model.constraints.add(vA[(_b, _m, _n, _k)] >= bitcA[(_b, _m, _n)] + self.u[_n][_k] - 1)

                    model.constraints.add(fA[(_m, _k)] == sum_a)

            for _m in range(m):
                for _k in range(k):

                    sum_b = 0
                    for _n in range(n):
                        for _b in range(_M):
                            sum_b += vB[(_b, _m, _n, _k)] * math.pow(2, _b)
                            model.constraints.add(vB[(_b, _m, _n, _k)] <= bitcB[(_b, _m, _n)])
                            model.constraints.add(vB[(_b, _m, _n, _k)] <= self.u[_n][_k])
                            model.constraints.add(vB[(_b, _m, _n, _k)] >= bitcB[(_b, _m, _n)] + self.u[_n][_k] - 1)

                    model.constraints.add(fB[(_m, _k)] == sum_b)

            for _n in range(n):
                for _k in range(k):
                    _sum = 0
                    for _m in range(m):
                        for _b in range(_M):
                            _sum += bitcA[(_b, _m, _n)] + bitcB[(_b, _m, _n)]
                    model.constraints.add(_sum >= self.u[_n][_k])

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

                    model.constraints.add(self.cA[_m][_n] == sum_a)
                    model.constraints.add(self.cB[_m][_n] == sum_b)
                    model.constraints.add(self.cA[_m][_n] + self.cB[_m][_n] <= ub)

        if mode_t == 'CARCH':
            # TODO: These loops can be collapsed once validation against C++ is complete
            for _m in range(m):
                for _k in range(k):
                    _sumA = 0
                    _sumB = 0
                    for _n in range(n):
                        if self._fixed_u[_n][_k] >= self.mu - self.tol:
                            _sumA += self.cA[_m][_n] * self._fixed_u[_n][_k]
                            _sumB += self.cB[_m][_n] * self._fixed_u[_n][_k]
                    model.constraints.add(fA[(_m, _k)] == _sumA)
                    model.constraints.add(fB[(_m, _k)] == _sumB)

                cluster_id = f_a.index[_m]
                # upper bound for solver
                ub = max(sum(copy_numbers.get(cluster_id, (0, 0))), cn_max)
                for _n in range(n):
                    model.constraints.add(self.cA[_m][_n] + self.cB[_m][_n] <= ub)

        if mode_t in ('FULL', 'CARCH'):
            for _m in range(m):
                model.constraints.add(self.cA[_m][0] == 1)
                model.constraints.add(self.cB[_m][0] == 1)

            if ampdel:
                for _m in range(m):
                    cluster_id = f_a.index[_m]
                    if cluster_id not in copy_numbers:
                        for _n in range(1, n):
                            model.constraints.add(self.cA[_m][_n] <= cn_max * adA[_m] + _base - _base * adA[_m])
                            model.constraints.add(self.cA[_m][_n] >= _base * adA[_m])
                            model.constraints.add(self.cB[_m][_n] <= cn_max * adB[_m] + _base - _base * adB[_m])
                            model.constraints.add(self.cB[_m][_n] >= _base * adB[_m])

        if mode_t == 'UARCH':
            # TODO: These loops can be collapsed once validation against C++ is complete
            for _m in range(m):
                for _k in range(k):
                    _sumA = 0
                    _sumB = 0
                    for _n in range(n):
                        _sumA += int(self._fixed_cA[_m][_n]) * self.u[_n][_k]
                        _sumB += int(self._fixed_cB[_m][_n]) * self.u[_n][_k]
                    model.constraints.add(fA[(_m, _k)] == _sumA)
                    model.constraints.add(fB[(_m, _k)] == _sumB)

        if mode_t in ('FULL', 'UARCH'):
            for _k in range(k):
                _sum = 0
                for _n in range(n):
                    _sum += self.u[_n][_k]
                model.constraints.add(_sum == 1)

        if (mode_t in ('FULL', 'UARCH')) and self.mu > 0:
            for _k in range(k):
                for _n in range(1, n):
                    model.constraints.add(x[(_n, _k)] >= self.u[_n][_k])
                    model.constraints.add(self.u[_n][_k] >= self.mu * x[(_n, _k)])

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
                        for _i in range(1, self.n - 1):
                            for _j in range(1, self.n):
                                model.constraints.add(
                                    bitcA[(_M, _m, _i)] - bitcA[(_M, _m, _j)] <= 2 - z[(_m, _i, _d)] - z[(_m, _j, _d)]
                                )
                                model.constraints.add(
                                    bitcA[(_M, _m, _j)] - bitcA[(_M, _m, _i)] <= 2 - z[(_m, _i, _d)] - z[(_m, _j, _d)]
                                )
                                model.constraints.add(
                                    bitcB[(_M, _m, _i)] - bitcB[(_M, _m, _j)] <= 2 - z[(_m, _i, _d)] - z[(_m, _j, _d)]
                                )
                                model.constraints.add(
                                    bitcB[(_M, _m, _j)] - bitcB[(_M, _m, _i)] <= 2 - z[(_m, _i, _d)] - z[(_m, _j, _d)]
                                )

            for _m in range(self.m):
                for _d in range(d - 1):
                    _sum_l = _sum_l1 = 0
                    for _n in range(1, self.n):
                        _sum_l += z[(_m, _n, _d)] * self.symmCoeff(_n)
                        _sum_l1 += z[(_m, _n, _d + 1)] * self.symmCoeff(_n)
                        model.constraints.add(_sum_l <= _sum_l1)

        if mode_t in ('FULL', 'CARCH'):
            self.build_symmetry_breaking(model)
            self.fix_given_cn(model)

        if mode_t == 'FULL':
            self.hot_start()

        obj = 0
        for _m in range(m):
            cluster_id = self.cluster_ids[_m]
            n_bins = len(self.binsA[cluster_id])

            for _k in range(k):
                sample_id = self.sample_ids[_k]
                for _i in range(n_bins):
                    obj += (yA[(_m, _k, _i)] + yB[(_m, _k, _i)]) * self.lengths[cluster_id][_i] / maes.loc[cluster_id, sample_id]

        model.obj = pe.Objective(expr=obj, sense=pe.minimize)
        self.model = model

        if pprint:
            print(str(self))
