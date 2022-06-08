import sys
from io import StringIO
import pytest
import shutil
import os
from mock import patch
import pandas as pd
from pandas.testing import assert_frame_equal

from hatchet import config
from hatchet.utils.cluster_bins_loc import main as cluster_bins_loc


this_dir = os.path.dirname(__file__)

@pytest.fixture(scope='module')
def output_folder():
    out = os.path.join(this_dir, 'out_cluster_loc')
    shutil.rmtree(out, ignore_errors=True)
    for sub_folder in ['bbc']:
        os.makedirs(os.path.join(out, sub_folder))
    return out

def test_cluster_bins_loc(output_folder):
    cluster_bins_loc(args=[
        f'{this_dir}/data/fl/bb/bulk.bb',
        '-o', f'{output_folder}/bbc/bulk.seg',
        '-O', f'{output_folder}/bbc/bulk.bbc',
        '--seed', '11111',
        '--transmat', 'diag',
        '--covar', 'diag',
        '--minK', '2',
        '--maxK', '10'
    ])

    seg1 = pd.read_csv(os.path.join(output_folder, 'bbc', 'bulk.seg'), sep = '\t')
    bbc1 = pd.read_csv(os.path.join(output_folder, 'bbc', 'bulk.bbc'), sep = '\t')

    seg2 = pd.read_csv(os.path.join(this_dir, 'data', 'vl', 'bbc', 'bulk.seg'), sep = '\t')
    bbc2 = pd.read_csv(os.path.join(this_dir, 'data', 'vl', 'bbc', 'bulk.bbc'), sep = '\t')
    assert_frame_equal(seg1, seg2)
    assert_frame_equal(bbc1, bbc2)


def test_cluster_bins_loc_singleton(output_folder):
    cluster_bins_loc(args=[
        f'{this_dir}/data/fl/bb/bulk.bb',
        '-o', f'{output_folder}/bbc/bulk_onecl.seg',
        '-O', f'{output_folder}/bbc/bulk_onecl.bbc',
        '--seed', '11111',
        '--transmat', 'diag',
        '--covar', 'diag',
        '--exactK', '1',
    ])

    seg1 = pd.read_csv(os.path.join(output_folder, 'bbc', 'bulk_onecl.seg'), sep = '\t')
    bbc1 = pd.read_csv(os.path.join(output_folder, 'bbc', 'bulk_onecl.bbc'), sep = '\t')

    seg2 = pd.read_csv(os.path.join(this_dir, 'data', 'vl', 'bbc', 'bulk_onecl.seg'), sep = '\t')
    bbc2 = pd.read_csv(os.path.join(this_dir, 'data', 'vl', 'bbc', 'bulk_onecl.bbc'), sep = '\t')
    assert_frame_equal(seg1, seg2)
    assert_frame_equal(bbc1, bbc2)