import pytest
import sys
import os
import glob
from io import StringIO
from mock import patch
import shutil
import pandas as pd
from pandas.testing import assert_frame_equal
from numpy.testing import assert_array_equal
import numpy as np

from hatchet import config
from hatchet.utils.count_reads import main as count_reads

this_dir = os.path.dirname(__file__)

@pytest.fixture(scope='module')
def bams():
    bam_directory = config.tests.bam_directory
    normal_bam = os.path.join(bam_directory, 'normal.bam')
    if not os.path.exists(normal_bam):
        pytest.skip('File not found: {}/{}'.format(bam_directory, normal_bam))
    tumor_bams = sorted([f for f in glob.glob(bam_directory + '/*.bam') if os.path.basename(f) != 'normal.bam'])
    if not tumor_bams:
        pytest.skip('No tumor bams found in {}'.format(bam_directory))

    return normal_bam, tumor_bams

@pytest.fixture(scope='module')
def output_folder():
    out = os.path.join(this_dir, 'parts')
    shutil.rmtree(out, ignore_errors=True)
    for sub_folder in ['rdr']:
        os.makedirs(os.path.join(out, sub_folder))
    return out

@patch('hatchet.utils.ArgParsing.extractChromosomes', return_value=['chr22'])
def test_script(_, bams, output_folder):
    normal_bam, tumor_bams = bams
    
    count_reads(
        args=[
            '-N', normal_bam,
            '-T'
        ] + tumor_bams + [
            '-st', config.paths.samtools,
            '-S', 'Normal', 'Tumor1', 'Tumor2', 'Tumor3',
            '-j', '12',
            '-b', os.path.join(this_dir, 'data', 'vl', 'baf', 'bulk.1bed'),
            '-O', os.path.join(output_folder, 'rdr'),
            '-V', 'hg19'
        ]
    )

    arr1 = np.loadtxt(os.path.join(output_folder, 'rdr', 'chr22.thresholds.gz'))
    arr2 = np.loadtxt(os.path.join(output_folder, 'rdr', 'chr22.total.gz'))
    arr3 = [l for l in open(os.path.join(output_folder, 'rdr', 'samples.txt'), 'r')]
    arr4 = pd.read_table(os.path.join(output_folder, 'rdr', 'total.tsv'))
    
    truth1 = np.loadtxt(os.path.join(this_dir, 'data', 'vl', 'rdr', 'chr22.thresholds.gz'))
    truth2 = np.loadtxt(os.path.join(this_dir, 'data', 'vl', 'rdr', 'chr22.total.gz'))
    truth3 = [l for l in open(os.path.join(this_dir, 'data', 'vl', 'rdr', 'samples.txt'), 'r')]
    truth4 = pd.read_table(os.path.join(this_dir, 'data', 'vl', 'rdr', 'total.tsv'))
  
    assert_array_equal(arr1, truth1)
    assert_array_equal(arr2, truth2)
    assert len(arr3) == len(truth3) and all([arr3[i] == truth3[i] for i in range(len(truth3))])
    assert_frame_equal(arr4, truth4)
    
    
                       

