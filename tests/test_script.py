import pytest
import sys
import os
from io import BytesIO as StringIO
from mock import patch
import hashlib
import shutil
import bnpy

import hatchet
from hatchet.utils.binBAM import main as binBAM
from hatchet.utils.deBAF import main as deBAF
from hatchet.utils.comBBo import main as comBBo
from hatchet.utils.cluBB import main as cluBB
from hatchet.bin.HATCHet import main as main


this_dir = os.path.dirname(__file__)
DATA_FOLDER = os.path.join(this_dir, 'data')
SOLVE = os.path.join(os.path.dirname(hatchet.__file__), 'solve')

# The cluBB.py script wants a path to the parent of the bnpy folder, which is now listed as a dependency for our code
BNPY_FOLDER = os.path.abspath(os.path.join(os.path.dirname(bnpy.__file__), '..'))


@pytest.fixture(scope='module')
def output_folder():
    out = os.path.join(this_dir, 'out')
    shutil.rmtree(out, ignore_errors=True)
    for sub_folder in ('bin', 'baf', 'bb', 'bbc', 'results', 'evaluation', 'analysis'):
        os.makedirs(os.path.join(out, sub_folder))
    return out


@pytest.mark.skipif(not os.getenv('HG19_FILE'), reason='HG19_FILE not set')
@pytest.mark.skipif(not os.getenv('SAMTOOLS_PATH'), reason='SAMTOOLS_PATH not set')
@pytest.mark.skipif(not os.getenv('BCFTOOLS_PATH'), reason='BCFTOOLS_PATH not set')
@patch('hatchet.utils.ArgParsing.extractChromosomes', return_value=['chr22'])
def test_script(_, output_folder):

    hg19_file = os.getenv('HG19_FILE')
    samtools_path = os.getenv('SAMTOOLS_PATH')
    bcftools_path = os.getenv('BCFTOOLS_PATH')

    binBAM(args=[
        '-N', os.path.join(DATA_FOLDER, 'SRR5906250-chr22.sorted.bam'),
        '-T', os.path.join(DATA_FOLDER, 'SRR5906251-chr22.sorted.bam'),
              os.path.join(DATA_FOLDER, 'SRR5906253-chr22.sorted.bam'),
        '-b', '50kb',
        '-st', samtools_path,
        '-S', 'Normal', 'TumorOP', 'Tumor2',
        '-g', hg19_file,
        '-j', '12',
        '-q', '11',
        '-O', os.path.join(output_folder, 'bin/normal.bin'),
        '-o', os.path.join(output_folder, 'bin/bulk.bin'),
        '-v'
    ])

    deBAF(args=[
        '-bt', bcftools_path,
        '-st', samtools_path,
        '-N', os.path.join(DATA_FOLDER, 'SRR5906250-chr22.sorted.bam'),
        '-T', os.path.join(DATA_FOLDER, 'SRR5906251-chr22.sorted.bam'),
              os.path.join(DATA_FOLDER, 'SRR5906253-chr22.sorted.bam'),
        '-S', 'Normal', 'TumorOP', 'Tumor2',
        '-r', hg19_file,
        '-j', '12',
        '-q', '11',
        '-Q', '11',
        '-U', '11',
        '-c', '8',
        '-C', '300',
        '-O', os.path.join(output_folder, 'baf/normal.baf'),
        '-o', os.path.join(output_folder, 'baf/bulk.baf'),
        '-v'
    ])

    _stdout = sys.stdout
    sys.stdout = StringIO()

    comBBo(args=[
        '-c', os.path.join(output_folder, 'bin/normal.bin'),
        '-C', os.path.join(output_folder, 'bin/bulk.bin'),
        '-B', os.path.join(output_folder, 'baf/bulk.baf'),
        '-m', 'MIRROR',
        '-e', '12'
    ])

    out = sys.stdout.getvalue()
    sys.stdout.close()
    sys.stdout = _stdout

    with open(os.path.join(output_folder, 'bb/bulk.bb'), 'w') as f:
        f.write(out)

    cluBB(args=[
        os.path.join(output_folder, 'bb/bulk.bb'),
        # Note: Since bnpy is installed as a dependency for Hatchet
        # The cluBB script should be modified to not worry about this argument at all
        '-by', BNPY_FOLDER,
        '-o', os.path.join(output_folder, 'bbc/bulk.seg'),
        '-O', os.path.join(output_folder, 'bbc/bulk.bbc'),
        '-e', '22171',  # random seed
        '-tB', '0.04',
        '-tR', '0.15',
        '-d', '0.4'  # 0.08 in script
    ])

    main(args=[
        SOLVE,
        '-x', os.path.join(output_folder, 'results'),
        '-i', os.path.join(output_folder, 'bbc/bulk'),
        '-n2',  # -n2,8 in script
        '-p', '400',
        '-v', '3',
        '-u', '0.03',
        '-r', '6700',  # random seed
        '-j', '8',
        '-eD', '6',
        '-eT', '12',
        '-g', '0.35',
        '-l', '0.6'
    ])

    assert hashlib.md5(os.path.join(output_folder, 'results/best.bbc.ucn')).hexdigest() == \
           'da20c70188fa1fa9be430b7182c301e5'
