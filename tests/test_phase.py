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
from hatchet.utils.phase_snps import main as phase_snps
from hatchet.utils.download_panel import main as download_panel

this_dir = os.path.dirname(__file__)


@pytest.fixture(scope='module')
def output_folder():
    out = os.path.join(this_dir, 'out', 'vl', 'phase')
    shutil.rmtree(out, ignore_errors=True)
    os.makedirs(out)
    return out


@patch('hatchet.utils.ArgParsing.extractChromosomes', return_value=['chr22'])
def test_script(_, output_folder):
    download_panel(
        args=[
            '-R', '1000GP_Phase3',
        ]
    )
    
    phase_snps(
        args=[
            '-g', config.paths.reference,
            '-V', 'hg19',
            '-N',
            '-o', os.path.join(output_folder, 'phase'),
            '-L', os.path.join(this_dir, 'data', 'vl', 'snps', 'chr22.vcf.gz')
            ]
    )

    df1 = pd.read_table(os.path.join(output_folder, 'phase', 'phased.vcf.gz'), comment = '#')
    df2 = pd.read_table(os.path.join(this_dir, 'data', 'vl', 'phase', 'phased.vcf.gz'), comment = '#')
    assert_frame_equal(df1, df2)

