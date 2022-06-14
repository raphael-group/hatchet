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
    copy_numbers = (
        OrderedDict()
    )  # dict from cluster_id => (cn_a, cn_b) 2-tuple
    clonal_parts = clonal.split(',')
    n_clonal_parts = len(clonal_parts)

    for i, c in enumerate(clonal_parts):
        cluster_id, cn_a, cn_b = [int(_c) for _c in c.split(':')]

        # The first two clonal clusters (used for scaling) MUST have different total copy numbers
        if i == 1:
            _first_cn_a, _first_cn_b = list(copy_numbers.values())[0]
            if _first_cn_a + _first_cn_b == cn_a + cn_b:
                raise ValueError(
                    'When >= 2 clonal copy numbers are given, the first two must be different in the two segmental clusters'
                )

        cn_total = cn_a + cn_b
        if (cn_total == 2) and (n_clonal_parts > 1):
            # warnings.warn('Please specify a single cluster when CN_A+CN_B=2')
            # TODO: The C++ implementation generates a warning corresponding to the above
            # This is suppressed by default, which is what we do here for now.
            pass
        if cluster_id in copy_numbers:
            raise ValueError(
                'Already encountered cluster_id =', str(cluster_id)
            )
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
        cn_sum_1, cn_sum_2 = sum(copy_numbers[cluster_id_1]), sum(
            copy_numbers[cluster_id_2]
        )
        purity = (
            2
            * (rdr_1 - rdr_2)
            / ((2 - cn_sum_2) * rdr_1 - (2 - cn_sum_1) * rdr_2)
        )

        purity[(1 <= purity) & (purity <= 1 + purity_tol)] = 1
        purity[(-purity_tol <= purity) & (purity <= 0)] = 0

        scale = (2 - (2 - cn_sum_1) * purity) / rdr_1
        assert np.all(
            (0 <= purity) & (purity <= 1) & (scale >= 0)
        ), 'scaling failed'

    return scale


def segmentation(
    cA,
    cB,
    u,
    cluster_ids,
    sample_ids,
    bbc_file,
    bbc_out_file=None,
    seg_out_file=None,
):
    df = pd.read_csv(bbc_file, sep='\t')
    # Chromosomes may or may not have chr notation - force a string dtype anyway
    df['#CHR'] = df['#CHR'].astype(str)
    # TODO: The legacy C++ implementation interprets the 'coverage' column as an int
    df['COV'] = df['COV'].astype(int)

    n_clone = len(cA[0])
    cA = pd.DataFrame(cA, index=cluster_ids, columns=range(n_clone))
    cB = pd.DataFrame(cB, index=cluster_ids, columns=range(n_clone))
    u = pd.DataFrame(u, index=range(n_clone), columns=sample_ids)
    # Make (n_sample, n_clone) in shape; easier to merge later
    u = u.T

    # copy-numbers represented as <CN_A>|<CN_B> strings
    cN = cA.astype(str) + '|' + cB.astype(str)
    cN.columns = ['cn_normal'] + [f'cn_clone{i}' for i in range(1, n_clone)]

    # Merge in copy-number + proportion information to our original Dataframe
    df = df.merge(cN, left_on='CLUSTER', right_index=True)
    u.columns = ['u_normal'] + [f'u_clone{i}' for i in range(1, n_clone)]
    df = df.merge(u, left_on='SAMPLE', right_index=True)

    # Sorting the values by start/end position critical for merging contiguous
    # segments with identical copy numbers later on
    df = df.sort_values(['#CHR', 'START', 'END', 'SAMPLE'])
    df = df.reset_index(drop=True)

    # last 2*n_clone columns names = [cn_normal, u_normal, cn_clone1, u_clone1, cn_clone2, ...]
    extra_columns = [
        col for sublist in zip(cN.columns, u.columns) for col in sublist
    ]
    all_columns = df.columns.values[: -2 * n_clone].tolist() + extra_columns

    if bbc_out_file is not None:
        # rearrange columns for easy comparison to legacy files
        df = df[all_columns]
        df.to_csv(bbc_out_file, sep='\t', index=False)

    if seg_out_file is not None:
        # create a new column that will use to store the contiguous segment number (1-indexed)
        df['segment'] = 0
        # all column names with cnA|cnB information (normal + clones)
        cN_column_names = cN.columns.tolist()
        # create a new column with all cnA|cnB strings joined as a single column
        df['all_copy_numbers'] = df[cN_column_names].apply(
            lambda x: ','.join(x), axis=1
        )

        _first_sample_name = df['SAMPLE'].iloc[0]

        # Grouping by consecutive identical values
        # See https://towardsdatascience.com/pandas-dataframe-group-by-consecutive-same-values-128913875dba
        group_name_to_indices = df.groupby(
            (
                # Find indices where we see the first sample name AND
                (df['SAMPLE'] == _first_sample_name)
                & (
                    # The chromosome changed values from the previous row OR
                    # any of the copy-numbers changed from the previous row OR
                    # the START changed from the END in the previous row
                    (df['#CHR'] != df['#CHR'].shift())
                    | (
                        df['all_copy_numbers']
                        != df['all_copy_numbers'].shift()
                    )
                    | (df['START'] != df['END'].shift())
                )
            ).cumsum(),
            # cumulative sum increments whenever a True is encountered, thus creating a series of monotonically
            # increasing values we can use as segment numbers
            sort=False,
        ).indices
        # 'indices' of a Pandas GroupBy object gives us a mapping from the group 'name'
        # (numbers starting from 1) -> indices in the Dataframe

        for group_name, indices in group_name_to_indices.items():
            df.loc[indices, 'segment'] = group_name

        aggregation_rules = {
            '#CHR': 'first',
            'START': 'min',
            'END': 'max',
            'SAMPLE': 'first',
        }
        aggregation_rules.update({c: 'first' for c in extra_columns})
        df = df.groupby(['segment', 'SAMPLE']).agg(aggregation_rules)

        df.to_csv(seg_out_file, sep='\t', index=False)
