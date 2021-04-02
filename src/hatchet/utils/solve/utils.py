import warnings
from collections import OrderedDict
import numpy as np
import pandas as pd


# A list of random states, used as a stack
random_states = []


class Random:
    """
    A context manager that pushes a random seed to the stack for reproducible results,
    and pops it on exit.
    """

    def __init__(self, seed=None):
        self.seed = seed

    def __enter__(self):
        if self.seed is not None:
            # Push current state on stack
            random_states.append(np.random.get_state())
            new_state = np.random.RandomState(self.seed)
            np.random.set_state(new_state.get_state())

    def __exit__(self, *args):
        if self.seed is not None:
            np.random.set_state(random_states.pop())


def parse_clonal(clonal):
    copy_numbers = OrderedDict()  # dict from cluster_id => (cn_a, cn_b) 2-tuple
    clonal_parts = clonal.split(',')
    n_clonal_parts = len(clonal_parts)
    cn_totals = set()

    for c in clonal_parts:
        cluster_id, cn_a, cn_b = [int(_c) for _c in c.split(':')]
        cn_total = cn_a + cn_b
        if (cn_total == 2) and (n_clonal_parts > 1):
            warnings.warn('Please specify a single cluster when CN_A+CN_B=2')
        if cn_total in cn_totals:
            raise ValueError('Cannot specify two clusters with same CN_A+CN_B')
        cn_totals.add(cn_total)
        if cluster_id in copy_numbers:
            raise ValueError('Already encountered cluster_id =', str(cluster_id))
        copy_numbers[cluster_id] = cn_a, cn_b

    return copy_numbers


def scale_rdr(rdr, copy_numbers, purity_tol=0.05):
    assert len(copy_numbers) >= 1
    if len(copy_numbers) == 1:
        diploid_cluster_id = list(copy_numbers.keys())[0]
        diploid_rdr = rdr.loc[diploid_cluster_id]
        scale = 2 / diploid_rdr
    else:
        # Use the first two copy-number specifications to determine scaling factors
        # Note the reversed assignment order to conform to C++ behavior - check with Simone!
        cluster_id_2, cluster_id_1 = tuple(copy_numbers.keys())[:2]  # TODO
        rdr_1, rdr_2 = rdr.loc[cluster_id_1], rdr.loc[cluster_id_2]
        cn_sum_1, cn_sum_2 = sum(copy_numbers[cluster_id_1]), sum(copy_numbers[cluster_id_2])
        purity = 2 * (rdr_1 - rdr_2) / ((2 - cn_sum_2) * rdr_1 - (2 - cn_sum_1) * rdr_2)

        purity[(1 <= purity) & (purity <= 1 + purity_tol)] = 1
        purity[(-purity_tol <= purity) & (purity <= 0)] = 0

        scale = (2 - (2 - cn_sum_1) * purity) / rdr_1
        assert np.all((0 <= purity) & (purity <= 1) & (scale >= 0)), 'scaling failed'

    rdr = scale * rdr
    return rdr


def segmentation(cA, cB, u, bbc_file, bbc_out_file=None, seg_out_file=None):
    df = pd.read_csv(bbc_file, sep='\t')
    # df = df.sort_values(['#CHR', 'START', 'END', 'SAMPLE'])

    n_cluster, n_clone = cA.shape

    # copy-numbers represented as <CN_A>|<CN_B> strings
    cN = cA.astype(str) + '|' + cB.astype(str)
    cN.columns = ['cn_normal'] + [f'cn_clone{i}' for i in range(1, n_clone)]

    # Merge in copy-number + proportion information to our original Dataframe
    df = df.merge(cN, left_on='CLUSTER', right_index=True)
    u = u.T  # Make (n_sample, n_clone) in shape; easier to merge later
    u.columns = ['u_normal'] + [f'u_clone{i}' for i in range(1, n_clone)]
    df = df.merge(u, left_on='SAMPLE', right_index=True)

    if bbc_out_file is not None:
        # rearrange columns and sort rows for easy comparison to legacy files
        # last 2*n_clone columns = [cn_normal, u_normal, cn_clone1, u_clone1, cn_clone2, ...]
        columns = df.columns.values[:-2*n_clone].tolist() + \
                  [col for sublist in zip(cN.columns, u.columns) for col in sublist]
        df = df[columns]
        df = df.sort_values(['#CHR', 'START', 'END', 'SAMPLE'])
        df.to_csv(bbc_out_file, sep='\t', index=False)

    return df


