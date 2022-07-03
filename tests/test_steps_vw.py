import os
import glob
from mock import patch
import shutil
import pandas as pd
from pandas.testing import assert_frame_equal
from numpy.testing import assert_array_equal
import numpy as np
import pytest

from hatchet import config
from hatchet.utils.count_reads import main as count_reads
from hatchet.utils.cluster_bins import main as cluster_bins
from hatchet.utils.combine_counts import main as combine_counts

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
    out = os.path.join(this_dir, 'out', 'vw')
    shutil.rmtree(out, ignore_errors=True)
    for sub_folder in ['rdr', 'bb', 'bbc']:
        os.makedirs(os.path.join(out, sub_folder))
    return out


@patch('hatchet.utils.ArgParsing.extractChromosomes', return_value=['chr22'])
def test_count_reads(_, bams, output_folder):
    normal_bam, tumor_bams = bams

    count_reads(
        args=['-N', normal_bam, '-T']
        + tumor_bams
        + [
            '-st',
            config.paths.samtools,
            '-S',
            'Normal',
            'Tumor1',
            'Tumor2',
            'Tumor3',
            '-j',
            '1',
            '-b',
            os.path.join(this_dir, 'data', 'vw', 'baf', 'bulk.1bed'),
            '-O',
            os.path.join(output_folder, 'rdr'),
            '-V',
            'hg19',
        ]
    )

    arr1 = np.loadtxt(os.path.join(output_folder, 'rdr', 'chr22.thresholds.gz'))
    arr2 = np.loadtxt(os.path.join(output_folder, 'rdr', 'chr22.total.gz'))
    arr3 = [l for l in open(os.path.join(output_folder, 'rdr', 'samples.txt'), 'r')]
    arr4 = pd.read_table(os.path.join(output_folder, 'rdr', 'total.tsv'))

    truth1 = np.loadtxt(os.path.join(this_dir, 'data', 'vw', 'rdr', 'chr22.thresholds.gz'))
    truth2 = np.loadtxt(os.path.join(this_dir, 'data', 'vw', 'rdr', 'chr22.total.gz'))
    truth3 = [l for l in open(os.path.join(this_dir, 'data', 'vw', 'rdr', 'samples.txt'), 'r')]
    truth4 = pd.read_table(os.path.join(this_dir, 'data', 'vw', 'rdr', 'total.tsv'))

    assert_array_equal(arr1, truth1)
    assert_array_equal(arr2, truth2)
    assert len(arr3) == len(truth3) and all([arr3[i] == truth3[i] for i in range(len(truth3))])
    assert_frame_equal(arr4, truth4)


@pytest.mark.skipif(not config.paths.reference, reason='paths.reference not set')
@patch('hatchet.utils.ArgParsing.extractChromosomes', return_value=['chr22'])
def test_combine_counts(_, output_folder):
    # Test without phasing
    combine_counts(
        args=[
            '-A',
            os.path.join(this_dir, 'data', 'vw', 'rdr'),
            '-t',
            os.path.join(this_dir, 'data', 'vw', 'rdr', 'total.tsv'),
            '-j',
            '1',
            '-b',
            os.path.join(this_dir, 'data', 'vw', 'baf', 'bulk.1bed'),
            '-o',
            os.path.join(output_folder, 'bb', 'bulk_nophase.bb'),
            '-V',
            'hg19',
        ]
    )
    df1 = pd.read_csv(os.path.join(output_folder, 'bb', 'bulk_nophase.bb'), sep='\t')
    df2 = pd.read_csv(os.path.join(this_dir, 'data', 'vw', 'bb', 'bulk_nophase.bb'), sep='\t')
    assert_frame_equal(df1, df2)

    # Test with phasing
    combine_counts(
        args=[
            '-A',
            os.path.join(this_dir, 'data', 'vw', 'rdr'),
            '-t',
            os.path.join(this_dir, 'data', 'vw', 'rdr', 'total.tsv'),
            '-j',
            '1',
            '-b',
            os.path.join(this_dir, 'data', 'vw', 'baf', 'bulk.1bed'),
            '-o',
            os.path.join(output_folder, 'bb', 'bulk_yesphase.bb'),
            '-V',
            'hg19',
            '-p',
            os.path.join(this_dir, 'data', 'vw', 'phase', 'phased.vcf.gz'),
        ]
    )

    df3 = pd.read_csv(os.path.join(output_folder, 'bb', 'bulk_yesphase.bb'), sep='\t')
    df4 = pd.read_csv(
        os.path.join(this_dir, 'data', 'vw', 'bb', 'bulk_yesphase.bb'),
        sep='\t',
    )
    assert_frame_equal(df3, df4)


def test_cluster_bins(output_folder):
    cluster_bins(
        args=[
            f'{this_dir}/data/fw/bb/bulk.bb',
            '-o',
            f'{output_folder}/bbc/bulk.seg',
            '-O',
            f'{output_folder}/bbc/bulk.bbc',
            '--seed',
            '11111',
            '--transmat',
            'diag',
            '--covar',
            'diag',
            '--minK',
            '2',
            '--maxK',
            '10',
        ]
    )

    seg1 = pd.read_csv(os.path.join(output_folder, 'bbc', 'bulk.seg'), sep='\t')
    bbc1 = pd.read_csv(os.path.join(output_folder, 'bbc', 'bulk.bbc'), sep='\t')

    seg2 = pd.read_csv(os.path.join(this_dir, 'data', 'vw', 'bbc', 'bulk.seg'), sep='\t')
    bbc2 = pd.read_csv(os.path.join(this_dir, 'data', 'vw', 'bbc', 'bulk.bbc'), sep='\t')
    assert_frame_equal(seg1, seg2)
    assert_frame_equal(bbc1, bbc2)


def test_cluster_bins_singleton(output_folder):
    cluster_bins(
        args=[
            f'{this_dir}/data/fw/bb/bulk.bb',
            '-o',
            f'{output_folder}/bbc/bulk_onecl.seg',
            '-O',
            f'{output_folder}/bbc/bulk_onecl.bbc',
            '--seed',
            '11111',
            '--transmat',
            'diag',
            '--covar',
            'diag',
            '--exactK',
            '1',
        ]
    )

    seg1 = pd.read_csv(os.path.join(output_folder, 'bbc', 'bulk_onecl.seg'), sep='\t')
    bbc1 = pd.read_csv(os.path.join(output_folder, 'bbc', 'bulk_onecl.bbc'), sep='\t')

    seg2 = pd.read_csv(os.path.join(this_dir, 'data', 'vw', 'bbc', 'bulk_onecl.seg'), sep='\t')
    bbc2 = pd.read_csv(os.path.join(this_dir, 'data', 'vw', 'bbc', 'bulk_onecl.bbc'), sep='\t')
    assert_frame_equal(seg1, seg2)
    assert_frame_equal(bbc1, bbc2)
