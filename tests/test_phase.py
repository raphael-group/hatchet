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
    out = os.path.join(this_dir, 'parts')
    shutil.rmtree(out, ignore_errors=True)
    for sub_folder in ['phase', 'panel']:
        os.makedirs(os.path.join(out, sub_folder))
    return out

<<<<<<< HEAD
@patch('hatchet.utils.ArgParsing.extractChromosomes', return_value=['chr22'])
def test_script(_, output_folder):    
=======

@pytest.mark.skipif(shutil.which('picard') is None, reason='picard is not available on PATH')
@pytest.mark.skipif(shutil.which('shapeit') is None, reason='shapeit is not available on PATH')
@patch('hatchet.utils.ArgParsing.extractChromosomes', return_value=['chr22'])
def test_script(_, output_folder):
>>>>>>> 2834f935f1745ae61ffeeaa017921add2cdceb2b
    download_panel(
        args=[
            '-D', os.path.join(output_folder, 'panel'),
            '-R', '1000GP_Phase3',
            '-V', 'hg19',
            '-N', 'true'
        ]
    )
    # TODO: test that panel was downloaded successfully
    
    phase_snps(
        args=[
            '-D', os.path.join(output_folder, 'panel'),
            '-g', config.paths.reference,
<<<<<<< HEAD
            '-si', config.paths.shapeit,
            '-pc', config.paths.picard,
            '-bg', config.paths.bgzip,
=======
>>>>>>> 2834f935f1745ae61ffeeaa017921add2cdceb2b
            '-V', 'hg19',
            '-N', 'true',
            '-o', os.path.join(output_folder, 'phase'),
            '-L', os.path.join(this_dir, 'data', 'test_parts', 'snps', 'chr22.vcf.gz')
            ]
    )
    
    
    df1 = pd.read_table(os.path.join(output_folder, 'phase', 'phased.vcf.gz'), comment = '#')
    df2 = pd.read_table(os.path.join(this_dir, 'data', 'test_parts', 'phase', 'phased.vcf.gz'), comment = '#')
    assert_frame_equal(df1, df2)
    
