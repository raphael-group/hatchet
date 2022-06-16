from multiprocessing import Process, Queue, JoinableQueue, Lock, Value
from collections import deque
import pysam

import hatchet.utils.ProgressBar as pb


def bin(
    samtools,
    samples,
    chromosomes,
    num_workers,
    q,
    size,
    regions,
    verbose=False,
):
    # Define a Lock and a shared value for log printing through ProgressBar
    err_lock = Lock()
    counter = Value('i', 0)
    progress_bar = pb.ProgressBar(
        total=len(samples) * len(chromosomes),
        length=40,
        lock=err_lock,
        counter=counter,
        verbose=verbose,
    )

    # Establish communication queues
    tasks = JoinableQueue()
    results = Queue()

    # Enqueue jobs
    jobs_count = 0
    for bam in samples:
        for chro in chromosomes:
            tasks.put((bam[0], bam[1], chro))
            jobs_count += 1

    # Setting up the workers
    workers = [
        Binner(tasks, results, progress_bar, samtools, q, size, regions, verbose)
        for i in range(min(num_workers, jobs_count))
    ]

    # Add a poison pill for each worker
    for i in range(len(workers)):
        tasks.put(None)

    # Start the workers
    for w in workers:
        w.start()

    # Wait for all of the tasks to finish
    tasks.join()

    # Get the results
    sorted_results = {}
    for i in range(jobs_count):
        res = results.get()
        sorted_results[res[0][0], res[0][1]] = res

    # Close Queues
    tasks.close()
    results.close()

    # Ensure each worker terminates
    for w in workers:
        w.terminate()
        w.join()

    return sorted_results


class Binner(Process):
    def __init__(
        self,
        task_queue,
        result_queue,
        progress_bar,
        samtools,
        q,
        size,
        regions,
        verbose,
    ):
        Process.__init__(self)
        self.task_queue = task_queue
        self.result_queue = result_queue
        self.progress_bar = progress_bar
        self.samtools = samtools
        self.q = q
        self.size = size
        self.regions = regions
        self.verbose = verbose

    def run(self):
        while True:
            next_task = self.task_queue.get()
            if next_task is None:
                # Poison pill means shutdown
                self.task_queue.task_done()
                break

            self.progress_bar.progress(
                advance=False,
                msg='{} starts on {} for {})'.format(self.name, next_task[1], next_task[2]),
            )
            bins = self.binChr(
                bamfile=next_task[0],
                samplename=next_task[1],
                chromosome=next_task[2],
            )
            self.progress_bar.progress(
                advance=True,
                msg='{} ends on {} for {})'.format(self.name, next_task[1], next_task[2]),
            )
            self.task_queue.task_done()
            self.result_queue.put(bins)
        return

    def binChr(self, bamfile, samplename, chromosome):
        read_callback = lambda r: r.mapping_quality >= self.q
        sam = pysam.AlignmentFile(bamfile, 'rb')
        bins = deque()
        for start, end in self.regions[chromosome]:
            while start < end:
                stop = min(start + self.size, end)
                if start > 0:
                    n_reads = sam.count(
                        chromosome,
                        start=start - 1,
                        stop=stop,
                        read_callback=read_callback,
                    )
                else:
                    n_reads = sam.count(
                        chromosome,
                        start=0,
                        stop=stop,
                        read_callback=read_callback,
                    )
                bins.append((samplename, chromosome, max(start, 1), stop, str(n_reads)))
                start += self.size
        return list(bins)
