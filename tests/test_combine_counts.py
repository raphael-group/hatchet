import pytest
import sys
import os
import glob
from io import StringIO
from mock import patch
import shutil
import pandas as pd
from pandas.testing import assert_frame_equal

import hatchet
from hatchet import config
from hatchet.utils.combine_counts import main as combine_counts


this_dir = os.path.dirname(__file__)
SOLVE = os.path.join(os.path.dirname(hatchet.__file__), 'solve')

@pytest.fixture(scope='module')
def output_folder():
    out = os.path.join(this_dir, 'out_abin_ph')
    shutil.rmtree(out, ignore_errors=True)
    for sub_folder in ['bb']:
        os.makedirs(os.path.join(out, sub_folder))
    return out

@pytest.mark.skipif(not config.paths.reference, reason='paths.reference not set')
@patch('hatchet.utils.ArgParsing.extractChromosomes', return_value=['chr22'])
def test_script(_, output_folder):
    # Test without phasing
    combine_counts(args=[
        '-A', os.path.join(this_dir, 'data', 'vl', 'rdr'),
        '-t',  os.path.join(this_dir, 'data', 'vl', 'rdr', 'total.tsv'),
        '-j', '1',
        '-b',  os.path.join(this_dir, 'data', 'vl', 'baf', 'bulk.1bed'),
        '-o',  os.path.join(output_folder, 'bb', 'bulk_nophase.bb'),
        '-V', 'hg19'
    ])
    df1 = pd.read_csv(os.path.join(output_folder, 'bb', 'bulk_nophase.bb'), sep='\t')
    df2 = pd.read_csv(os.path.join(this_dir, 'data', 'vl', 'bb', 'bulk_nophase.bb'), sep = '\t')
    assert_frame_equal(df1, df2)
    
    # Test with phasing
    combine_counts(args=[
        '-A', os.path.join(this_dir, 'data', 'vl', 'rdr'),
        '-t',  os.path.join(this_dir, 'data', 'vl', 'rdr', 'total.tsv'),
        '-j', '1',
        '-b',  os.path.join(this_dir, 'data', 'vl', 'baf', 'bulk.1bed'),
        '-o',  os.path.join(output_folder, 'bb', 'bulk_yesphase.bb'),
        '-V', 'hg19', 
        '-p', os.path.join(this_dir, 'data', 'vl', 'phase', 'phased.vcf.gz')
    ])

    df3 = pd.read_csv(os.path.join(output_folder, 'bb', 'bulk_yesphase.bb'), sep='\t')
    df4 = pd.read_csv(os.path.join(this_dir, 'data', 'vl', 'bb', 'bulk_yesphase.bb'), sep = '\t')
    assert_frame_equal(df3, df4)
    


