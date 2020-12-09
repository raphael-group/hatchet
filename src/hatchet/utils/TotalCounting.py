import os
import sys
import shlex
import subprocess
from multiprocessing import Process, Queue, JoinableQueue, Lock, Value

from . import ProgressBar as pb
from . import Supporting as sp


def tcount(samtools, samples, chromosomes, num_workers, q, verbose=False):
    # Define a Lock and a shared value for log printing through ProgressBar
    err_lock = Lock()
    counter = Value('i', 0)
    progress_bar = pb.ProgressBar(total=len(samples)*len(chromosomes), length=40, lock=err_lock, counter=counter, verbose=verbose)

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
    workers = [TotalCounter(tasks, results, progress_bar, samtools, q, verbose) for i in range(min(num_workers, jobs_count))]

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
        sorted_results[res[0], res[1]] = res[2]

    # Close Queues
    tasks.close()
    results.close()

    # Ensure each worker terminates
    for w in workers:
        w.terminate()
        w.join()

    return sorted_results



class TotalCounter(Process):

    def __init__(self, task_queue, result_queue, progress_bar, samtools, q, verbose):
        Process.__init__(self)
        self.task_queue = task_queue
        self.result_queue = result_queue
        self.progress_bar = progress_bar
        self.samtools = samtools
        self.q = q
        self.verbose = verbose

    def run(self):
        while True:
            next_task = self.task_queue.get()
            if next_task is None:
                # Poison pill means shutdown
                self.task_queue.task_done()
                break

            self.progress_bar.progress(advance=False, msg="{} starts total counting on (sample={}, chromosome={})".format(self.name, next_task[1], next_task[2]))
            count = self.binChr(bamfile=next_task[0], samplename=next_task[1], chromosome=next_task[2])
            self.progress_bar.progress(advance=True, msg="{} ends total counting on (sample={}, chromosome={})".format(self.name, next_task[1], next_task[2]))
            self.task_queue.task_done()
            self.result_queue.put(count)
        return

    def binChr(self, bamfile, samplename, chromosome):
        popen = subprocess.Popen
        pipe = subprocess.PIPE
        split = shlex.split
        cmd = "{} view {} -c -q {} {}".format(self.samtools, bamfile, self.q, chromosome)
        stdout, stderr = popen(split(cmd), stdout=pipe, stderr=pipe, text=True).communicate()
        if stderr != "":
            self.progress_bar.progress(advance=False, msg="{}{}: samtools warns \"{}\"on (sample={}, chromosome={}){}".format(sp.bcolors.WARNING, self.name, stderr, samplename, chromosome, sp.bcolors.ENDC))
        return (samplename, chromosome, int(stdout.strip()))
