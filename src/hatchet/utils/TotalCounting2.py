from concurrent import futures
from BAMBinning2 import work


def tcount(samtools, samples, chromosomes, num_workers, q, verbose=False):
    with futures.ProcessPoolExecutor(max_workers=num_workers) as executor:
        results = {}
        to_do_map = {}
        for sample in samples:
            samfile, samplename = sample
            for chromosome in chromosomes:
                future = executor.submit(work, samfile, samplename, chromosome, q, None, None)
                to_do_map[future] = samplename, chromosome

        for future in futures.as_completed(to_do_map):
            results[to_do_map[future]] = int(future.result()[0][-1])

    return results
