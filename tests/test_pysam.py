import os
import hashlib
import pysam


this_dir = os.path.dirname(__file__)
DATA_FOLDER = os.path.join(this_dir, 'data')


def test_pysam():
    samfile = os.path.join(DATA_FOLDER, 'SRR5906250-chr22.sorted.bam')
    assert hashlib.md5(open(samfile, 'rb').read()).hexdigest() == 'a0f00352db430e606c601b850dc80abb'
    # sam = pysam.AlignmentFile(samfile, 'rb')