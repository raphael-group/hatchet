import sys
from io import StringIO
import pytest
import os
import glob
import gzip
import hashlib
from mock import patch
import shutil
import pickle
import pandas as pd
from pandas.testing import assert_frame_equal

from hatchet import config
from hatchet.utils import ArgParsing
from hatchet.utils.count_reads_fw import main as count_reads
from hatchet.utils.genotype_snps import main as genotype_snps
from hatchet.utils.count_alleles import counting
from hatchet.utils.combine_counts_fw import main as combine_counts
from hatchet.utils.cluster_bins import main as cluster_bins
from hatchet.utils.plot_bins import main as plot_bins
from hatchet.bin.HATCHet import main as main
from hatchet.utils.plot_cn import main as plot_cn
from hatchet.utils.solve import solver_available

this_dir = os.path.dirname(__file__)


@pytest.fixture(scope='module')
def bams():
    bam_directory = config.tests.bam_directory
    normal_bam = os.path.join(bam_directory, 'normal.bam')
    if not os.path.exists(normal_bam):
        pytest.skip('File not found: {}'.format(os.path.join(bam_directory, 'normal.bam')))
    tumor_bams = sorted([f for f in glob.glob(bam_directory + '/*.bam') if os.path.basename(f) != 'normal.bam'])
    if not tumor_bams:
        pytest.skip('No tumor bams found in {}'.format(bam_directory))

    return normal_bam, tumor_bams


@pytest.fixture(scope='module')
def normal_snps():
    with open(f'{this_dir}/data/fl/snps/normal_snps.pik', 'rb') as f:
        return pickle.load(f)


@pytest.fixture(scope='module')
def output_folder():
    out = os.path.join(this_dir, 'out')
    shutil.rmtree(out, ignore_errors=True)
    for sub_folder in (
        'bin',
        'snps',
        'baf',
        'bb',
        'bbc',
        'plots',
        'results',
        'evaluation',
        'analysis',
    ):
        os.makedirs(os.path.join(out, sub_folder))
    return out


@pytest.mark.skipif(not config.paths.reference, reason='paths.reference not set')
@patch('hatchet.utils.ArgParsing.extractChromosomes', return_value=['chr22'])
@patch(
    'hatchet.utils.count_reads_fw.knownRegions',
    return_value={'chr22': [(30256931, 32622323)]},
)
def test_count_reads(_mock0, _mock1, bams, output_folder):
    normal_bam, tumor_bams = bams

    count_reads(
        args=['-N', normal_bam, '-T']
        + tumor_bams
        + [
            '-b',
            '50kb',
            '-st',
            config.paths.samtools,
            '-S',
            'Normal',
            'Tumor1',
            'Tumor2',
            'Tumor3',
            '-g',
            config.paths.reference,
            '-j',
            '12',
            '-q',
            '11',
            '-O',
            os.path.join(output_folder, 'bin/normal.1bed'),
            '-o',
            os.path.join(output_folder, 'bin/bulk.1bed'),
            '-v',
        ]
    )

    df1 = pd.read_csv(f'{output_folder}/bin/normal.1bed', sep='\t')
    df2 = pd.read_csv(f'{this_dir}/data/fl/bin/normal.1bed', sep='\t')
    assert_frame_equal(df1, df2)

    df1 = pd.read_csv(f'{output_folder}/bin/bulk.1bed', sep='\t')
    df2 = pd.read_csv(f'{this_dir}/data/fl/bin/bulk.1bed', sep='\t')
    assert_frame_equal(df1, df2)


@pytest.mark.skipif(not config.paths.reference, reason='paths.reference not set')
@patch(
    'hatchet.utils.ArgParsing.extractChromosomes',
    return_value=['chr22:30256931-32622323'],
)
def test_genotype_snps(_mock1, bams, output_folder):
    normal_bam, _ = bams

    genotype_snps(
        args=[
            '-N',
            normal_bam,
            '-r',
            config.paths.reference,
            '-c',
            '290',  # min reads
            '-C',
            '300',  # max reads
            '-o',
            f'{output_folder}/snps',
            '-st',
            config.paths.samtools,
            '-bt',
            config.paths.bcftools,
            '-j',
            '12',
        ]
    )

    with gzip.open(f'{output_folder}/snps/chr22:30256931-32622323.vcf.gz', 'rb') as f:
        # ignore commented lines since these may have timestamps and software version numbers etc.
        lines = '\n'.join([line for line in f.read().decode('utf8').split('\n') if not line.startswith('#')])
        assert hashlib.md5(lines.encode('utf8')).hexdigest() == '3d81c51d21c22334ce1fc069cb005328'


@patch(
    'hatchet.utils.ArgParsing.extractChromosomes',
    return_value=['chr22:30256931-32622323'],
)
def test_count_alleles_normal_snps(_mock1, bams, normal_snps, output_folder):
    normal_bam, tumor_bams = bams

    args = (
        ['-N', normal_bam, '-T']
        + tumor_bams
        + [
            '-S',
            'Normal',
            'Tumor1',
            'Tumor2',
            'Tumor3',
            '-r',
            config.paths.reference,
            '-j',
            '12',
            '-q',
            '3',
            '-Q',
            '3',
            '-U',
            '3',
            '-c',
            '8',
            '-C',
            '300',
            '-O',
            f'{output_folder}/baf/normal.1bed',
            '-o',
            f'{output_folder}/baf/bulk.1bed',
            '-L',
            f'{this_dir}/data/fl/snps/chr22:30256931-32622323.vcf.gz',
        ]
    )

    args = ArgParsing.parse_count_alleles_arguments(args)

    snps = counting(
        bcftools=args['bcftools'],
        reference=args['reference'],
        samples=[args['normal']],
        chromosomes=args['chromosomes'],
        num_workers=args['j'],
        snplist=args['snps'],
        q=args['q'],
        Q=args['Q'],
        mincov=args['mincov'],
        dp=args['maxcov'],
        E=args['E'],
        verbose=False,
        outdir=args['outputSnps'],
    )

    assert snps == normal_snps


def test_combine_counts(output_folder):
    _stdout = sys.stdout
    sys.stdout = StringIO()

    combine_counts(
        args=[
            '-c',
            f'{this_dir}/data/fl/bin/normal.1bed',
            '-C',
            f'{this_dir}/data/fl/bin/bulk.1bed',
            '-B',
            f'{this_dir}/data/fl/baf/bulk.1bed',
            '-e',
            '12',
        ]
    )

    out = sys.stdout.getvalue()
    sys.stdout.close()
    sys.stdout = _stdout

    with open(f'{this_dir}/data/fl/bb/bulk.bb', 'r') as f:
        assert out == f.read()


def test_cluster_bins(output_folder):
    cluster_bins(
        args=[
            f'{this_dir}/data/fl/bb/bulk.bb',
            '-o',
            f'{output_folder}/bbc/bulk.seg',
            '-O',
            f'{output_folder}/bbc/bulk.bbc',
            '-e',
            '22171',  # random seed
            '-tB',
            '0.04',
            '-tR',
            '0.15',
            '-d',
            '0.4',
            '-K',
            '20',
        ]
    )

    df1 = pd.read_csv(f'{output_folder}/bbc/bulk.seg', sep='\t')
    df2 = pd.read_csv(f'{this_dir}/data/fl/bbc/bulk.seg', sep='\t')
    assert_frame_equal(df1, df2)


def test_plot_bins(output_folder):
    # We simply check if we're able to run plot_bins without exceptions
    plot_bins(
        args=[
            os.path.join(f'{this_dir}/data/fl/bbc/bulk.bbc'),
            '--rundir',
            os.path.join(output_folder, 'plots'),
        ]
    )


def test_compute_cn(output_folder):
    if solver_available():
        main(
            args=[
                '-x',
                f'{output_folder}/results',
                '-i',
                f'{this_dir}/data/fl/bbc/bulk',
                '-n2',
                '-p',
                '100',
                '-v',
                '3',
                '-u',
                '0.03',
                '-r',
                '6700',  # random seed
                '-j',
                '8',
                '-eD',
                '6',
                '-eT',
                '12',
                '-g',
                '0.35',
                '-l',
                '0.6',
            ]
        )

        df1 = pd.read_csv(f'{output_folder}/results/best.bbc.ucn', sep='\t')
        df2 = pd.read_csv(f'{this_dir}/data/fl/results/best.bbc.ucn', sep='\t')
        assert_frame_equal(df1, df2)


def test_plot_cn(output_folder):
    # We simply check if we're able to run plot_cn without exceptions
    plot_cn(
        args=[
            f'{this_dir}/data/fl/results/best.bbc.ucn',
            '--rundir',
            f'{output_folder}/evaluation',
        ]
    )
