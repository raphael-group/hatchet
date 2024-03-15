import os
import numpy as np
import pandas as pd
from pyomo import environ as pe
from collections import OrderedDict

from hatchet.utils.solve.ilp_subset import ILPSubset, ILPSubsetSplit
from hatchet.utils.solve.cd import CoordinateDescent, CoordinateDescentSplit
from hatchet.utils.solve.utils import parse_clonal, scale_rdr
from hatchet import config
import hatchet.utils.Supporting as sp


def solver_available(solver=None):
    solver = solver or config.compute_cn.solver
    if solver == 'cpp':
        return os.getenv('GRB_LICENSE_FILE') is not None
    elif solver == 'gurobipy':
        return pe.SolverFactory('gurobi', solver_io='python').available(exception_flag=False)
    return pe.SolverFactory(solver).available(exception_flag=False)


# Compute the Mahalanobis distance for each point "a" in S with regularization
def mahalanobis_distances(S, centroid, regularization=1e-6):
    num_points = len(S)
    mahalanobis_sum = 0

    # Compute the mean and covariance matrix of S
    mean = centroid
    # we take diagonal of covariance matrix
    # gamma only scales one of the dimensions
    cov_matrix = np.cov(S, rowvar=False)
    if S.shape[1] == 1:
        diagonal_cov_matrix = cov_matrix.reshape(1, 1)
    else:
        diagonal_cov_matrix = np.diag(np.diag(cov_matrix))

    # Add regularization to the covariance matrix
    cov_matrix_regularized = diagonal_cov_matrix + regularization * np.identity(diagonal_cov_matrix.shape[0])

    for i in range(num_points):
        a = S[i]
        diff = a - mean
        mahalanobis_distance = np.sqrt(np.dot(diff, np.dot(np.linalg.inv(cov_matrix_regularized), diff)))
        mahalanobis_sum += mahalanobis_distance

    return mahalanobis_sum


def compute_mahalanobis_objective(cA, cB, U, gamma, bbc):
    fA = np.dot(cA, U)
    fB = np.dot(cB, U)
    rCent = (fA + fB) / list(gamma)[0]
    bCent = fB / (fA + fB)

    # find balanced clusters (i.e. rows with identical cA and cB)
    balanced_clusters = []
    for i in range(len(cA)):
        if np.array_equal(cA[i], cB[i]):
            balanced_clusters.append(i)

    objective = 0
    per_cluster_maha = dict()
    # groupby cluster
    for cluster_id, _df in bbc.groupby('CLUSTER'):
        rdr = _df['RD']
        baf = _df['BAF']
        S = np.array([rdr, baf]).T
        if S.shape[0] <= 1:
            continue
        # compute mahalanobis distance for each cluster
        if cluster_id not in balanced_clusters:
            increment = mahalanobis_distances(
                S, np.concatenate((rCent[cluster_id - 1], bCent[cluster_id - 1]), axis=0), regularization=1e-7
            )
        else:
            increment = mahalanobis_distances(S[:, 0].reshape(-1, 1), rCent[cluster_id - 1], regularization=1e-7)
        objective += increment
        per_cluster_maha[cluster_id] = increment
    return objective, per_cluster_maha


def solve(
    clonal,
    bbc_file,
    seg_file,
    n,
    solver='gurobi',
    solve_mode='cd',
    d=-1,
    cn_max=-1,
    mu=0.01,
    diploid_threshold=0.1,
    ampdel=True,
    n_seed=400,
    n_worker=8,
    random_seed=None,
    max_iters=None,
    timelimit=None,
    binwise=False,
    evolcons=False,
    bp_max=60,
    uniqueclones=False,
    purities=None,
):

    assert solve_mode in ('ilp', 'cd', 'both'), 'Unrecognized solve_mode'
    assert solver_available(solver), f'Solver {solver} not available or not licensed'

    if max_iters is None:
        max_iters = 10

    df = pd.read_csv(seg_file, sep='\t')
    df = df.sort_values(['#ID', 'SAMPLE'])

    bbc = pd.read_table(bbc_file)
    bbc = bbc.sort_values(by=['#CHR', 'START', 'SAMPLE'])

    # sanity-check
    sample_ids = np.sort(df['SAMPLE'].unique())
    for _cluster_id, _df in df.groupby('#ID'):
        _sample_ids = _df['SAMPLE'].values
        if not np.all(_sample_ids == sample_ids):
            raise ValueError(f'Sample IDs for cluster {_cluster_id} do not match {sample_ids}')

    rdr = df.pivot(index='#ID', columns='SAMPLE', values='RD')
    baf = df.pivot(index='#ID', columns='SAMPLE', values='BAF')

    #  compute adjacency matrix for breakpoint distance constraint

    # verts = list(df["#ID"].unique()) + [-1]
    # clus_adj_mat = {}
    # for i in verts:
    #     for j in verts:
    #         clus_adj_mat.setdefault(i, {})[j] = 0

    #  get the first sample name and extract cluster ids in the genomic order
    clust_id_ord = list(bbc[bbc['SAMPLE'] == bbc['SAMPLE'].values[0]]['CLUSTER'])

    clus_adj_list = []
    #  telomere adjacencies -1
    clus_adj_list.append((-1, clust_id_ord[0]))
    clus_adj_list.append((-1, clust_id_ord[-1]))

    for i in range(1, len(clust_id_ord)):
        c1, c2 = sorted([clust_id_ord[i - 1], clust_id_ord[i]])
        if c1 != c2:
            clus_adj_list.append((c1, c2))

    clus_adj_counts = dict()
    for i in clus_adj_list:
        clus_adj_counts[i] = clus_adj_counts.get(i, 0) + 1

    bins = {}  # cluster_id => no. of bins
    for cluster_id, _df in df.groupby('#ID'):
        if len(_df['#BINS'].unique()) != 1:
            raise ValueError(f'Bin sizes for cluster {cluster_id} across tumor samples are not identical!')
        bins[cluster_id] = _df.iloc[0]['#BINS']
    bins = pd.Series(bins)

    weights = 100 * bins / sum(bins)

    if clonal is None:
        _candidate_cluster_ids = np.all(rdr > 0.5 - diploid_threshold, axis=1)
        if not np.any(_candidate_cluster_ids):
            raise RuntimeError(f'Unable to determine cluster with diploid RDR threshold {diploid_threshold}')
        dipoid_cluster_id = (_candidate_cluster_ids * weights).idxmax()
        copy_numbers = {dipoid_cluster_id: (1, 1)}
    else:
        copy_numbers = parse_clonal(clonal)

    gamma = scale_rdr(rdr, copy_numbers)
    sp.log(
        msg='Computed scaling factor gamma = \n' + str(gamma) + '\n',
        level='INFO',
    )
    rdr = rdr * gamma
    f_b = rdr * baf
    f_a = rdr - f_b

    # compute mean absolute errors for each cluster and sample using
    # with respect to fa and fb (fractional copy numbers).
    # use rdr and baf values in the bbc file

    # bbc_rdr = bbc.pivot(index='CLUSTER', columns='SAMPLE', values='RD')
    # bbc_baf = bbc.pivot(index='CLUSTER', columns='SAMPLE', values='BAF')
    # bbc_rdr = bbc_rdr * gamma

    # bbc_fb = bbc_rdr * bbc_baf
    # bbc_fa = bbc_rdr - bbc_fb

    # # compute means of fb and fa for each cluster and sample
    # mean_fa = bbc_fa.mean(axis=2)
    # mean_fb = bbc_fb.mean(axis=2)

    # maes = np.mean(np.abs(bbc_fa - mean_fa), axis=2) + np.mean(np.abs(bbc_fb - mean_fb), axis=2)

    # create a table called maes that has the same rows and index as f_a
    maes = pd.DataFrame(index=f_a.index, columns=f_a.columns)

    for _cluster_id, _df in bbc.groupby('CLUSTER'):
        for _sample, _df2 in _df.groupby('SAMPLE'):
            _rdr = _df2['RD'].values * gamma[_sample]
            _baf = _df2['BAF'].values
            _f_b = _rdr * _baf
            _f_a = _rdr - _f_b

            mean_fa = np.mean(_f_a)
            mean_fb = np.mean(_f_b)
            # compute mean absolute error
            maes.loc[_cluster_id, _sample] = np.mean(np.abs(_f_a - mean_fa)) + np.mean(np.abs(_f_b - mean_fb))
            if maes.loc[_cluster_id, _sample] == 0:
                maes.loc[_cluster_id, _sample] = 1

    if not binwise:
        if solve_mode == 'ilp':
            ilp = ILPSubset(
                n,
                cn_max,
                d=d,
                mu=mu,
                ampdel=ampdel,
                copy_numbers=copy_numbers,
                f_a=f_a,
                f_b=f_b,
                w=weights,
                clus_adj_counts=clus_adj_counts,
                maes=maes,
                evolcons=evolcons,
                bp_max=bp_max,
                uniqueclones=uniqueclones,
                purities=purities,
            )
            ilp.create_model(pprint=True)
            solution = ilp.run(solver_type=solver, timelimit=timelimit)
            return solution
        elif solve_mode == 'cd':
            cd = CoordinateDescent(
                f_a=f_a,
                f_b=f_b,
                n=n,
                mu=mu,
                d=d,
                cn_max=cn_max,
                w=weights,
                ampdel=ampdel,
                cn=copy_numbers,
                clus_adj_counts=clus_adj_counts,
                maes=maes,
                evolcons=evolcons,
                bp_max=bp_max,
                uniqueclones=uniqueclones,
                purities=purities,
            )
            solution = cd.run(
                solver_type=solver,
                max_iters=max_iters,
                n_seed=n_seed,
                j=n_worker,
                random_seed=random_seed,
                timelimit=timelimit,
            )

            # objective, pcmaha = compute_mahalanobis_objective(solution[1], solution[2], solution[3], gamma, bbc)
            # sp.log(
            #     msg=f'\nMahanalobis Objective value: {objective} \n {pcmaha}\n',
            #     level='INFO',
            # )

            return solution
        else:
            cd = CoordinateDescent(
                f_a=f_a,
                f_b=f_b,
                n=n,
                mu=mu,
                d=d,
                cn_max=cn_max,
                w=weights,
                ampdel=ampdel,
                cn=copy_numbers,
                clus_adj_counts=clus_adj_counts,
                maes=maes,
                evolcons=evolcons,
                bp_max=bp_max,
                uniqueclones=uniqueclones,
                purities=purities,
            )
            _, cA, cB, _, _, _ = cd.run(
                solver_type=solver,
                max_iters=max_iters,
                n_seed=n_seed,
                j=n_worker,
                random_seed=random_seed,
                timelimit=timelimit,
            )

            ilp = ILPSubset(
                n,
                cn_max,
                d=d,
                mu=mu,
                ampdel=ampdel,
                copy_numbers=copy_numbers,
                f_a=f_a,
                f_b=f_b,
                w=weights,
                clus_adj_counts=clus_adj_counts,
                maes=maes,
                evolcons=evolcons,
                bp_max=bp_max,
                uniqueclones=uniqueclones,
                purities=purities,
            )
            ilp.create_model()
            ilp.hot_start(cA, cB)
            solution = ilp.run(solver_type=solver, timelimit=timelimit)
            return solution

    else:
        bins = OrderedDict()  # cluster_id => RDR for cluster
        for cluster_id, _df in df.groupby('#ID'):
            if len(_df['#BINS'].unique()) != 1:
                raise ValueError(f'Bin sizes for cluster {cluster_id} across tumor samples are not identical!')
            my_bbc = bbc[bbc.CLUSTER == cluster_id]
            if _df.iloc[0]['#BINS'] * len(_df) != len(my_bbc):   # seg file should have 1 row per sample
                raise ValueError(f'BBC and SEG files describe inconsisitent # bins for cluster {cluster_id}!')
            bins[cluster_id] = my_bbc.RD

        binned_length = (bbc[bbc.SAMPLE == bbc.iloc[0].SAMPLE].END - bbc[bbc.SAMPLE == bbc.iloc[0].SAMPLE].START).sum()

        bins_rdr = {
            k: my_df.pivot(index=['#CHR', 'START'], columns='SAMPLE', values='RD')
            for k, my_df in bbc.groupby('CLUSTER')
        }
        bins_baf = {
            k: my_df.pivot(index=['#CHR', 'START'], columns='SAMPLE', values='BAF')
            for k, my_df in bbc.groupby('CLUSTER')
        }
        bins_length = {
            k: (
                (
                    my_df.pivot(index=['#CHR', 'START'], columns='SAMPLE', values='END')
                    - my_df.pivot(
                        index=['#CHR', 'START'],
                        columns='SAMPLE',
                        values='START',
                    )
                ).values[:, 0]
                * 100
            )
            / binned_length
            for k, my_df in bbc.groupby('CLUSTER')
        }

        binsA = {}
        binsB = {}
        for k in bins_rdr.keys():
            bins_rdr[k] = bins_rdr[k] * gamma
            binsB[k] = bins_rdr[k] * bins_baf[k]
            binsA[k] = bins_rdr[k] - binsB[k]

        if solve_mode == 'ilp':
            ilp = ILPSubsetSplit(
                n,
                cn_max,
                d=d,
                mu=mu,
                ampdel=ampdel,
                copy_numbers=copy_numbers,
                f_a=f_a,
                f_b=f_b,
                binsA=binsA,
                binsB=binsB,
                lengths=bins_length,
                clus_adj_counts=clus_adj_counts,
                purities=purities,
            )
            ilp.create_model(pprint=True)
            return ilp.run(solver_type=solver, timelimit=timelimit)
        elif solve_mode == 'cd':
            cd = CoordinateDescentSplit(
                f_a=f_a,
                f_b=f_b,
                n=n,
                mu=mu,
                d=d,
                cn_max=cn_max,
                ampdel=ampdel,
                cn=copy_numbers,
                binsA=binsA,
                binsB=binsB,
                lengths=bins_length,
            )
            return cd.run(
                solver_type=solver,
                max_iters=max_iters,
                n_seed=n_seed,
                j=n_worker,
                random_seed=random_seed,
                timelimit=timelimit,
            )
        else:
            cd = CoordinateDescentSplit(
                f_a=f_a,
                f_b=f_b,
                n=n,
                mu=mu,
                d=d,
                cn_max=cn_max,
                ampdel=ampdel,
                cn=copy_numbers,
                binsA=binsA,
                binsB=binsB,
                lengths=bins_length,
            )
            _, cA, cB, _, _, _ = cd.run(
                solver_type=solver,
                max_iters=max_iters,
                n_seed=n_seed,
                j=n_worker,
                random_seed=random_seed,
                timelimit=timelimit,
            )

            ilp = ILPSubsetSplit(
                n,
                cn_max,
                d=d,
                mu=mu,
                ampdel=ampdel,
                copy_numbers=copy_numbers,
                f_a=f_a,
                f_b=f_b,
                cn=copy_numbers,
                binsA=binsA,
                binsB=binsB,
                lengths=bins_length,
            )
            ilp.create_model()
            ilp.hot_start(cA, cB)
            return ilp.run(solver_type=solver, timelimit=timelimit)
