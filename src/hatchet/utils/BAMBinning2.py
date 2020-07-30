import pysam
from concurrent import futures


def bin(samtools, samples, chromosomes, num_workers, q, size, regions, verbose=False):
    with futures.ProcessPoolExecutor(max_workers=num_workers) as executor:
        results = {}
        to_do_map = {}
        for sample in samples:
            samfile, samplename = sample
            for chromosome in chromosomes:
                future = executor.submit(work, samfile, samplename, chromosome, q, regions[chromosome], size)
                to_do_map[future] = samplename, chromosome

        for future in futures.as_completed(to_do_map):
            results[to_do_map[future]] = future.result()

    return results


def work(samfile, name, chrom, min_quality, regions=None, size=None):

    def region_string(chrom, start=None, end=None):
        if start is None and end is None:
            return chrom
        if start in (None, 0):
            start = 1
        end = '' if end is None else ('-' + str(end))

        return chrom + ':' + str(start) + str(end)

    sam = pysam.AlignmentFile(samfile, 'rb')

    if regions is None:
        reads = sam.fetch(region=region_string(chrom))
        n = 0
        for read in reads:
            if read.mapping_quality >= min_quality:
                n += 1
        results = [(name, chrom, None, None, str(n))]
    else:
        results = []
        for region in regions:
            start, end = region
            while start < end:
                _end = min(end, start + size)
                reads = sam.fetch(region=region_string(chrom, start, _end))
                n = 0
                for read in reads:
                    if read.mapping_quality >= min_quality:
                        n += 1

                results.append((name, chrom, start, _end, str(n)))
                start += size

    return results
