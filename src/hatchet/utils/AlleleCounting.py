import os
import sys
import shlex
import subprocess
from multiprocessing import Process, Queue, JoinableQueue, Lock, Value

from . import ProgressBar as pb



def count(samtools, bcftools, reference, samples, chromosomes, num_workers, q, Q, mincov, dp, E, snplist=None, verbose=False):
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
    workers = [AlleleCounter(tasks, results, progress_bar, samtools, bcftools, reference, q, Q, mincov, dp, E, snplist, verbose) for i in range(min(num_workers, jobs_count))]

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


class AlleleCounter(Process):

    def __init__(self, task_queue, result_queue, progress_bar, samtools, bcftools, reference, q, Q, mincov, dp, E, snplist=None, verbose=False):
        Process.__init__(self)
        self.task_queue = task_queue
        self.result_queue = result_queue
        self.progress_bar = progress_bar
        self.samtools = samtools
        self.bcftools = bcftools
        self.reference = reference
        self.q = q
        self.Q = Q
        self.mincov = mincov
        self.dp = dp
        self.E = E
        self.snplist = snplist
        self.verbose = verbose

    def run(self):
        while True:
            next_task = self.task_queue.get()
            if next_task is None:
                # Poison pill means shutdown
                self.task_queue.task_done()
                break

            self.progress_bar.progress(advance=False, msg="{} starts allele counting on (sample={}, chromosome={})".format(self.name, next_task[1], next_task[2]))
            snps = self.countAlleles(bamfile=next_task[0], samplename=next_task[1], chromosome=next_task[2])
            self.progress_bar.progress(advance=True, msg="{} ends allele counting on (sample={}, chromosome={})".format(self.name, next_task[1], next_task[2]))
            self.task_queue.task_done()
            self.result_queue.put(snps)
        return

    def countAlleles(self, bamfile, samplename, chromosome):
        cmd_mpileup = "{} mpileup {} -ugf {} -q {} -Q {} --skip-indels -t INFO/AD,AD,DP".format(self.samtools, bamfile, self.reference, self.q, self.Q)
        cmd_query = "{} query -f '%CHROM,%POS,%AD{},%AD{}\n' -i 'DP<={} & SUM(AD)>={}'".format(self.bcftools, "{0}", "{1}", self.dp, self.mincov)
        if chromosome is not None:
            cmd_mpileup += " -r {}".format(chromosome)
        if self.snplist is not None:
            cmd_mpileup += " -l {}".format(self.snplist)
        if self.E:
            cmd_mpileup += " -E"
        if self.verbose:
            err = open("samtools.log", 'a')
        else:
            err = open(os.devnull, 'w')
        mpileup = subprocess.Popen(shlex.split(cmd_mpileup), stdout=subprocess.PIPE, stderr=err, shell=False, text=True)
        query = subprocess.Popen(shlex.split(cmd_query), stdin=mpileup.stdout, stdout=subprocess.PIPE, stderr=err, shell=False, text=True)
        stdout, stderr = query.communicate()
        err.close()

        return [(samplename, el[0], el[1], int(el[2]), int(el[3])) for el in (tuple(line.split(',')) for line in stdout.strip().split('\n') if line != "")]



def naiveCount(samtools, samples, chromosomes, num_workers, q, Q, E, snplist=None, verbose=False):
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
    workers = [NaiveAlleleCounter(tasks, results, progress_bar, samtools, q, Q, E, snplist, verbose) for i in range(min(num_workers, jobs_count))]

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



class NaiveAlleleCounter(Process):

    def __init__(self, task_queue, result_queue, progress_bar, samtools, q, Q, E, snplist=None, verbose=False):
        Process.__init__(self)
        self.task_queue = task_queue
        self.result_queue = result_queue
        self.progress_bar = progress_bar
        self.samtools = samtools
        self.q = q
        self.Q = Q
        self.E = E
        self.snplist = snplist
        self.verbose = verbose

    def run(self):
        while True:
            next_task = self.task_queue.get()
            if next_task is None:
                # Poison pill means shutdown
                self.task_queue.task_done()
                break
            self.progress_bar.progress(advance=False, msg="{} starts allele counting on (sample={}, chromosome={})".format(self.name, next_task[1], next_task[2]))
            allelecount = self.parsePileup(bamfile=next_task[0], samplename=next_task[1], chromosome=next_task[2])
            self.progress_bar.progress(msg="{} ends allele counting on (sample={}, chromosome={})".format(self.name, next_task[1], next_task[2]))
            self.task_queue.task_done()
            self.result_queue.put(allelecount)
        return


    def parsePileup(self, bamfile, samplename, chromosome):
        # Creating pileup
        cmd_mpileup = "{} mpileup {} -q {} -Q {} --skip-indels -t INFO/AD,AD,DP".format(self.samtools, bamfile, self.q, self.Q)
        if chromosome is not None:
            cmd_mpileup += " -r {}".format(chromosome)
        if self.snplist is not None:
            cmd_mpileup += " -l {}".format(self.snplist)
        if self.E:
            cmd_mpileup += " -E"
        if self.verbose:
            err = open("samtools.log", 'a')
        else:
            err = open(os.devnull, 'w')
        mpileup = subprocess.Popen(shlex.split(cmd_mpileup), stdout=subprocess.PIPE, stderr=err, shell=False, text=True)
        stdout, stderr = mpileup.communicate()

        # Counting alleles
        return [self.parsePileupRecord(record, samplename) for record in stdout.strip().split('\n') if record != ""]


    def parsePileupRecord(self, record, sample):## , minQuality=0):
        fields = record.strip().split()
        chromosome = fields[0]
        position = fields[1]
        coverage = int(fields[3])
        if coverage == 0:
            return (sample, chromosome, position, 0, 0, 0, 0, 0)

        bases = fields[4].replace(",", fields[2].lower()).replace(".", fields[2])
        ##qualities = [x-33 for x in map(ord, fields[5])]
        index = 0
        skip = 0
        counts = {"A" : 0, "T" : 0, "C" : 0, "G" : 0, "+" : 0, "-" : 0}

        for ichar in range(len(bases)):
            if skip == 0:
                char = bases[ichar].upper()
                if char in ("A", "T", "C", "G"):
                    ##if (qualities[index] > minQuality): counts[char.upper()] += 1
                    counts[char] += 1
                    index += 1
                elif char in ("<", ">", "*"):
                    index += 1
                    counts["-"] += 1
                elif char in ("^"):
                    skip = 1
                elif char in ("$"):
                    skip = 0
                elif char in ("+", "-"):
                    jump = ""
                    offset = 0
                    while bases[ichar+1].isdigit():
                        ichar += 1
                        offset += 1
                        jump += bases[ichar]
                    jump = int(jump)
                    skip = jump + offset
                    counts["+"] += jump
                else:
                    raise ValueError("Unexpected element \"{}\" found in the following pileup record: ".format(char))
            else:
                skip -= 1

        assert index == coverage
        assert (counts["A"] + counts["T"] + counts["C"] + counts["G"] + counts["-"]) == coverage
        return (sample, chromosome, position, coverage, counts["A"], counts["T"], counts["C"], counts["G"])
