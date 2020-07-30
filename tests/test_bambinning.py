from hatchet.utils.BAMBinning2 import bin


def test_bins():
    bins = bin(
        samtools='/opt/raphael-group/samtools/bin/samtools',
        samples=[('/media/vineetb/t5-vineetb/raphael-group/data/hatchet/SRR5906250.sorted.bam', 'Normal')],
        chromosomes=['chr19', 'chr22'],
        num_workers=2,
        q=11,
        size=50000,
        regions={'chr19': [(0, 59128983)], 'chr22': [(0, 51304566)]},
        verbose=True
    )

    assert ('Normal', 'chr19') in bins
    bins19 = bins[('Normal', 'chr19')]
    assert len(bins19) == 1183

    assert bins19[0] == ('Normal', 'chr19', 0, 50000, '0')
    assert bins19[1] == ('Normal', 'chr19', 50000, 100000, '81')

    assert ('Normal', 'chr22') in bins
    bins22 = bins[('Normal', 'chr22')]
    assert len(bins22) == 1027

    assert bins22[0] == ('Normal', 'chr22', 0, 50000, '0')
    assert bins22[1] == ('Normal', 'chr22', 50000, 100000, '0')
