#!/usr/bin/python2

import os, sys
import os.path
import argparse
import shlex
import subprocess as pr
from multiprocessing import Process, Queue, JoinableQueue, Lock, Value
from scipy.stats import beta
# from statsmodels.stats.proportion import *

import ArgParsing as ap
from Supporting import *
import Supporting as sp
import ProgressBar as pb



def main(args=None):
    log(msg="# Parsing the input arguments, checking the consistency of given files, and extracting required information\n", level="STEP")
    args = ap.parse_baf_arguments(args)
    logArgs(args, 80)

    log(msg="# Inferring SNPs from the normal sample\n", level="STEP")
    snps = call(bcftools=args["bcftools"], reference=args["reference"], samples=[args["normal"]], chromosomes=args["chromosomes"], 
                num_workers=args["j"], q=args["q"], Q=args["Q"], mincov=args["mincov"], dp=args["maxcov"], 
                E=args["E"], snplist=args["snps"], outdir=args['outputsnps'], verbose=args["verbose"])
    
    log(msg="# Counting number of identified SNPs\n", level="STEP")
    def count(f):
        cmd_bcf = 'bcftools query -f \'%CHROM,%POS\n\' {}'.format(f)
        cmf_wcl = 'wc -l'
        bcf = pr.Popen(shlex.split(cmd_bcf), stdout=pr.PIPE, stderr=pr.PIPE)
        return pr.Popen(shlex.split(cmd_wcl), stdin=bcf.stdout(), stdout=pr.PIPE, stderr=pr.PIPE).communicate()[0]
    number_snps = sum(count(f) for f in snps)
    if number_snps == 0:
        raise ValueError(sp.error("No SNPs found in the normal!\n"))
    else:
        log(msg="{} SNPs have been identified in total\n".format(number_snps), level="INFO")
    
    
def call(bcftools, reference, samples, chromosomes, num_workers, q, Q, mincov, dp, E, outdir, snplist=None, verbose=False):
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
    workers = [Caller(tasks, results, progress_bar, bcftools, reference, q, Q, mincov, dp, E, outdir, snplist, verbose) for i in range(min(num_workers, jobs_count))]

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

    # Close Queues
    tasks.close()
    results.close()

    # Ensure each worker terminates
    for w in workers:
        w.terminate()
        w.join()

    return sorted_results


class Caller(Process):

    def __init__(self, task_queue, result_queue, progress_bar, bcftools, reference, q, Q, mincov, dp, E, outdir, snplist, verbose):
        Process.__init__(self)
        self.task_queue = task_queue
        self.result_queue = result_queue
        self.progress_bar = progress_bar
        self.bcftools = bcftools
        self.reference = reference
        self.q = q
        self.Q = Q
        self.mincov = mincov
        self.dp = dp
        self.E = E
        self.outdir = outdir
        self.snplist = snplist
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
        cmd_mpileup = "{} mpileup {} -Ou -f {} --skip-indels -a INFO/AD,AD,DP -q {} -Q {} -d {}".format(self.bcftools, bamfile, self.reference, self.q, self.Q, self.dp)
        cmd_call = "{} call -mv -Oz -o {}".format(self.bcftools, os.path.join(self.outdir, '{}.vcf.gz'.format(chromosome)))
        if chromosome is not None:
            cmd_mpileup += " -r {}".format(chromosome)
        if self.snplist is not None:
            cmd_mpileup += " -R {}".format(self.snplist)
        if self.E:
            cmd_mpileup += " -E"
        errname = os.path.join(outdir, "{}_{}_bcftools.log".format(samplename, chromosome))
        with open(errname, 'w') as err:
            mpileup = pr.Popen(shlex.split(cmd_mpileup), stdout=pr.PIPE, stderr=err)
            call = pr.Popen(shlex.split(cmd_call), stdin=mpileup.stdout, stdout=pr.PIPE, stderr=err)
            stdout, stderr = call.communicate()
        if mpileup.returncode != 0 or call.returncode != 0:
            raise ValueError(sp.error('SNP Calling failed on {} of {}, please check errors in {}!').format(chromosome, samplename, errname))
        else:
            os.remove(errname)
        return os.path.join(self.outdir, '{}.vcf.gz'.format(chromosome))


if __name__ == '__main__':
    main()
