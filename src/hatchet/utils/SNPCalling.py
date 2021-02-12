import os
import sys
import shlex
import subprocess
from multiprocessing import Process, Queue, JoinableQueue, Lock, Value

import hatchet.utils.ProgressBar as pb



def call(samtools, bcftools, reference, samples, chromosomes, num_workers, q, Q, qual, mincov, dp, E, regions, snplist=None, verbose=False):
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
    workers = [SNPCaller(tasks, results, progress_bar, samtools, bcftools, reference, q, Q, qual, mincov, dp, E, snplist, regions, verbose) for i in range(min(num_workers, jobs_count))]

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
        if len(res) > 0:
            sorted_results[res[0][0], res[0][1]] = res

    # Close Queues
    tasks.close()
    results.close()

    # Ensure each worker terminates
    for w in workers:
        w.terminate()
        w.join()

    return sorted_results


class SNPCaller(Process):

    def __init__(self, task_queue, result_queue, progress_bar, samtools, bcftools, reference, q, Q, qual, mincov, dp, E, snplist=None, regions=None, verbose=False):
        Process.__init__(self)
        self.task_queue = task_queue
        self.result_queue = result_queue
        self.progress_bar = progress_bar
        self.samtools = samtools
        self.bcftools = bcftools
        self.reference = reference
        self.q = q
        self.Q = Q
        self.qual = qual
        self.mincov = mincov
        self.dp = dp
        self.E = E
        self.snplist = snplist
        self.regions = regions
        self.verbose = verbose

    def run(self):
        while True:
            next_task = self.task_queue.get()
            if next_task is None:
                # Poison pill means shutdown
                self.task_queue.task_done()
                break

            self.progress_bar.progress(advance=False, msg="{} starts SNP calling on (sample={}, chromosome={})".format(self.name, next_task[1], next_task[2]))
            snps = self.callSNPs(bamfile=next_task[0], samplename=next_task[1], chromosome=next_task[2])
            self.progress_bar.progress(advance=True, msg="{} ends SNP calling on (sample={}, chromosome={})".format(self.name, next_task[1], next_task[2]))
            self.task_queue.task_done()
            self.result_queue.put(snps)
        return

    def callSNPs(self, bamfile, samplename, chromosome):
        cmd_mpileup = "{} mpileup {} -ugf {} -q {} -Q {} --skip-indels -t INFO/AD,AD,DP".format(self.samtools, bamfile, self.reference, self.q, self.Q)
        cmd_call = "{} call -vmO u --skip-variants indels --variants-only".format(self.bcftools)
        # cmd_filter = "{} filter -s LowQual -e '%QUAL<{} || DP>{}'".format(self.bcftools, self.qual, self.dp)
        cmd_query = "{} query -f '%CHROM,%POS,%AD{},%AD{}\n' -i 'AN==2 & AC==1 & %QUAL>={} & DP<={} & SUM(AD)>={}'".format(self.bcftools,"{0}","{1}", self.qual, self.dp, self.mincov)
        if chromosome is not None:
            cmd_mpileup += " -r {}".format(chromosome)
        assert not (self.snplist != None and self.regions != None), "SNP list and regions cannot be both provided!"
        if self.snplist is not None:
            cmd_mpileup += " -l {}".format(self.snplist)
        if self.regions is not None:
            cmd_mpileup += " -l {}".format(self.regions)
        if self.E:
            cmd_mpileup += " -E"
        if self.verbose:
            err = open("samtools.log", 'a')
        else:
            err = open(os.devnull, 'w')
        mpileup = subprocess.Popen(shlex.split(cmd_mpileup), stdout=subprocess.PIPE, stderr=err, shell=False)
        call = subprocess.Popen(shlex.split(cmd_call), stdin=mpileup.stdout, stdout=subprocess.PIPE, stderr=err, shell=False)
        query = subprocess.Popen(shlex.split(cmd_query), stdin=call.stdout, stdout=subprocess.PIPE, stderr=err, shell=False)
        stdout, stderr = query.communicate()
        err.close()

        return [(samplename, el[0], el[1], int(el[2]), int(el[3])) for el in (tuple(line.split(',')) for line in stdout.strip().split('\n') if line != "")]
