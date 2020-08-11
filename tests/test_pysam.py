import os
import hashlib
import pysam


this_dir = os.path.dirname(__file__)
DATA_FOLDER = os.path.join(this_dir, 'data')


def test_pysam():
    samfile = os.path.join(DATA_FOLDER, 'SRR5906250-chr22.sorted.bam')
    assert hashlib.md5(samfile).hexdigest() == '5a7beb6deb0beabfcf64108ee10fa690'
    # sam = pysam.AlignmentFile(samfile, 'rb')