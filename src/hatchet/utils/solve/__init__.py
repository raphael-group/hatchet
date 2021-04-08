import numpy as np
import pandas as pd

from .ilp_subset import ILPSubset
from .cd import CoordinateDescent
from .utils import parse_clonal, scale_rdr


def solve(clonal, seg_file, n, solve_mode='cd', d=-1, cn_max=-1, mu=0.01, diploid_threshold=0.1, ampdel=True, n_seed=400, n_worker=8):

    assert solve_mode in ('ilp', 'cd', 'both'), 'Unrecognized solve_mode'

    df = pd.read_csv(seg_file, sep='\t')
    df = df.sort_values(['#ID', 'SAMPLE'])

    # sanity-check
    sample_ids = np.sort(df['SAMPLE'].unique())
    for _cluster_id, _df in df.groupby('#ID'):
        _sample_ids = _df['SAMPLE'].values
        if not np.all(_sample_ids == sample_ids):
            raise ValueError(f'Sample IDs for cluster {_cluster_id} do not match {sample_ids}')

    rdr = df.pivot(index='#ID', columns='SAMPLE', values='RD')
    baf = df.pivot(index='#ID', columns='SAMPLE', values='BAF')

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

    rdr = scale_rdr(rdr, copy_numbers)
    f_b = rdr * baf
    f_a = rdr - f_b

    if solve_mode == 'ilp':
        ilp = ILPSubset(n, cn_max, d=d, mu=mu, ampdel=ampdel, copy_numbers=copy_numbers, f_a=f_a, f_b=f_b, w=weights)
        ilp.create_model()
        return ilp.run()
    elif solve_mode == 'cd':
        cd = CoordinateDescent(f_a=f_a, f_b=f_b, n=n, mu=mu, d=d, cn_max=cn_max, w=weights, ampdel=ampdel,
                               cn=copy_numbers)
        return cd.run(max_iters=10, n_seed=n_seed, j=n_worker)
    else:
        cd = CoordinateDescent(f_a=f_a, f_b=f_b, n=n, mu=mu, d=d, cn_max=cn_max, w=weights, ampdel=ampdel,
                               cn=copy_numbers)
        _, cA, cB, _ = cd.run(max_iters=10, n_seed=n_seed, j=n_worker)

        ilp = ILPSubset(n, cn_max, d=d, mu=mu, ampdel=ampdel, copy_numbers=copy_numbers, f_a=f_a, f_b=f_b, w=weights)
        ilp.create_model()
        ilp.hot_start(cA, cB)
        return ilp.run()
