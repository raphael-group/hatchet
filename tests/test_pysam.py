import os
import pysam


this_dir = os.path.dirname(__file__)
DATA_FOLDER = os.path.join(this_dir, 'data')


def test_pysam():
    samfile = os.path.join(DATA_FOLDER, 'SRR5906250-chr22.sorted.bam')
    sam = pysam.AlignmentFile(samfile, 'rb')