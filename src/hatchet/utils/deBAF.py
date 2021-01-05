#!/usr/bin/python2

import sys, os
import os.path
import argparse
import shlex
import subprocess as pr
from multiprocessing import Process, Queue, JoinableQueue, Lock, Value
from scipy.stats import beta
# from statsmodels.stats.proportion import *


import ProgressBar as pb
import ArgParsing as ap
from Supporting import *
import Supporting as sp


def main(args=None):
    log(msg="# Parsing the input arguments, checking the consistency of given files, and extracting required information\n", level="STEP")
    args = ap.parse_baf_arguments(args)
    logArgs(args, 80)

    log(msg="# Counting SNPs alleles from the matched-normal sample\n", level="STEP")
    snps = counting(bcftools=args["bcftools"], reference=args["reference"], samples=[args["normal"]], chromosomes=args["chromosomes"], num_workers=args["j"], 
                    snplist=args["snps"], q=args["q"], Q=args["Q"], mincov=args["mincov"], dp=args["maxcov"], E=args["E"], 
                    verbose=args["verbose"], outdir=args['outputSnps'])

    log(msg="# Selecting heterozygous SNPs\n", level="STEP")
    hetSNPs = selectHetSNPs(counts=snps, gamma=args["gamma"], maxshift=args["maxshift"])
    if not hetSNPs: sp.close("No heterozygous SNPs found in the selected regions of the normal!\n")

    log(msg="# Writing the list of selected SNPs, covered and heterozygous in the normal sample\n", level="STEP")
    hetsnpsfiles = {}
    for chro in args["chromosomes"]:
        hetsnpsfiles[chro] = os.path.join(args["outputSnps"], 'TMP_{}.tsv'.format(chro))
        with open(hetsnpsfiles[chro], 'w') as f:
            if (args["normal"][1], chro) in hetSNPs:
                for snp in sorted(hetSNPs[args["normal"][1], chro]):
                    f.write("{}\t{}\n".format(chro, snp))

    log(msg="# Writing the allele counts of the normal sample for selected SNPs\n", level="STEP")
    if args["outputNormal"] is not None:
        with open(args["outputNormal"], 'w') as f:
            for chro in args["chromosomes"]:
                if (args["normal"][1], chro) in hetSNPs:
                    for snp in sorted(hetSNPs[args["normal"][1], chro]):
                        count = hetSNPs[args["normal"][1], chro][snp]
                        f.write("{}\t{}\t{}\t{}\t{}\n".format(args["normal"][1], chro, snp, count[0][1], count[1][1]))
    else:
        for chro in args["chromosomes"]:
            if (args["normal"][1], chro) in hetSNPs:
                for snp in sorted(hetSNPs[args["normal"][1], chro]):
                    count = hetSNPs[args["normal"][1], chro][snp]
                    sys.stdout.write("{}\t{}\t{}\t{}\t{}\n".format(args["normal"][1], chro, snp, count[0][1], count[1][1]))

    log(msg="# Counting SNPs alleles from tumour samples\n", level="STEP")
    rcounts = counting(bcftools=args["bcftools"], reference=args["reference"], samples=args["samples"], chromosomes=args["chromosomes"], num_workers=args["j"], 
                       snplist=hetsnpsfiles, q=args["q"], Q=args["Q"], mincov=args["mincov"], dp=args["maxcov"], E=args["E"], 
                       verbose=args["verbose"], outdir=args['outputSnps'])
    if not rcounts: sp.close("The selected SNPs are not covered in the tumors!\n")
    rcounts = {c : dict(map(lambda r : (int(r[2]), dict(r[3])), rcounts[c])) for c in rcounts}
    het = (lambda chro : hetSNPs[args["normal"][1], chro])
    form = (lambda REF, ALT, T : ((REF, T[REF] if REF in T else 0), (ALT, T[ALT] if ALT in T else 0)))
    counts = {c : {o : form(het(c[1])[o][0][0], het(c[1])[o][1][0], rcounts[c][o]) for o in rcounts[c]} for c in rcounts}

    log(msg="# Writing the allele counts of tumor samples for selected SNPs\n", level="STEP")
    map(lambda f : os.remove(hetsnpsfiles[f]), hetsnpsfiles)
    if args["outputTumors"] is not None:
        with open(args["outputTumors"], 'w') as f:
            for sample in args["samples"]:
                for chro in args["chromosomes"]:
                    if (sample[1], chro) in counts:
                        for snp in counts[sample[1], chro]:
                            count = counts[sample[1], chro][snp]
                            f.write("{}\t{}\t{}\t{}\t{}\n".format(sample[1], chro, snp, count[0][1], count[1][1]))
    else:
        for sample in args["samples"]:
            for chro in args["chromosomes"]:
                if (sample[1], chro) in counts:
                    for snp in counts[sample[1], chro]:
                        count = counts[sample[1], chro][snp]
                        sys.stdout.write("{}\t{}\t{}\t{}\t{}\n".format(sample[1], chro, snp, count[0][1], count[1][1]))


def selectHetSNPs(counts, gamma, maxshift):
    getma = (lambda d : max(d.items(), key=(lambda x : x[1])))
    hetSNPs = {c : [(r[0], r[1], r[2], r[3][0], max(r[3][1:], key=(lambda x : x[1]))) for r in counts[c] if len(r[3]) >= 2] for c in counts}
    check = (lambda r : isHet(r[3][1], r[4][1], gamma) and checkShift(r[3][1], r[4][1], maxshift))
    hetSNPs = {c : filter(check, hetSNPs[c]) for c in hetSNPs if len(hetSNPs[c]) > 0}
    return {c : {int(r[2]) : (r[3], r[4]) for r in reversed(hetSNPs[c])} for c in hetSNPs if len(hetSNPs[c]) > 0}


def isHet(countA, countB, gamma):
    p_lower = gamma / 2.0
    p_upper = 1.0 - p_lower
    [c_lower, c_upper] = beta.ppf([p_lower, p_upper], countA + 1, countB + 1)
    return c_lower <= 0.5 <= c_upper

# def isHet(countA, countB, gamma):
#     lb, ub = proportion_confint(min(countA, countB), countA+countB, alpha=gamma, method='beta')
#     return lb <= 0.5 <= ub

# def isHet(countA, countB, gamma):
#     lb, ub = proportion_confint(min(countA, countB), countA+countB, alpha=gamma, method='jeffreys')
#     return lb <= 0.5 <= ub


def checkShift(countA, countB, maxshift):
    return (0.5 - (float(min(countA, countB)) / float(countA + countB)) ) <= maxshift


def counting(bcftools, reference, samples, chromosomes, num_workers, snplist, q, Q, mincov, dp, E, verbose, outdir):
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
    workers = [AlleleCounter(tasks, results, progress_bar, bcftools, reference, q, Q, mincov, dp, E, snplist, verbose, outdir) for i in range(min(num_workers, jobs_count))]

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

    def __init__(self, task_queue, result_queue, progress_bar, bcftools, reference, q, Q, mincov, dp, E, snplist, verbose, outdir):
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
        self.snplist = snplist
        self.verbose = verbose
        self.outdir = outdir

    def run(self):
        while True:
            next_task = self.task_queue.get()
            if next_task is None:
                # Poison pill means shutdown
                self.task_queue.task_done()
                break

            self.progress_bar.progress(advance=False, msg="{} starts on {} for {}".format(self.name, next_task[1], next_task[2]))
            snps = self.countAlleles(bamfile=next_task[0], samplename=next_task[1], chromosome=next_task[2])
            self.progress_bar.progress(advance=True, msg="{} ends on {} for {}".format(self.name, next_task[1], next_task[2]))
            self.task_queue.task_done()
            self.result_queue.put(snps)
        return

    def countAlleles(self, bamfile, samplename, chromosome):
        cmd_mpileup = "{} mpileup {} -Ou -f {} --skip-indels -a INFO/AD -q {} -Q {} -d {}".format(self.bcftools, bamfile, self.reference, self.q, self.Q, self.dp)
        cmd_mpileup = "{} -r {} -R {}".format(cmd_mpileup, chromosome, self.snplist[chromosome])
        cmd_query = "{} query -f '%CHROM\\t%POS\\t%REF,%ALT\\t%AD\\n' -i 'SUM(AD)<={} & SUM(AD)>={}'".format(self.bcftools, self.dp, self.mincov)
        if self.E:
            cmd_mpileup += " -E"
        errname = os.path.join(self.outdir, "{}_{}_bcftools.log".format(samplename, chromosome))
        with open(errname, 'w') as err:
            mpileup = pr.Popen(shlex.split(cmd_mpileup), stdout=pr.PIPE, stderr=err)
            query = pr.Popen(shlex.split(cmd_query), stdin=mpileup.stdout, stdout=pr.PIPE, stderr=err)
            stdout, _ = query.communicate()
            codes = map(lambda p : p.wait(), [mpileup, query])
        if any(c != 0 for c in codes):
            raise ValueError(sp.error('Allele counting failed on {} of {}, please check errors in {}!').format(chromosome, samplename, errname))
        else:
            os.remove(errname)
        alleles = {'A', 'C', 'G', 'T'}
        mkcounts = (lambda p, q : map(lambda y : (y[0], int(y[1])), filter(lambda x : x[0] in alleles, zip(p, q))))
        form = (lambda p : (samplename, p[0], p[1], mkcounts(p[2].split(','), p[3].split(','))))
        return [form(line.strip().split()) for line in stdout.strip().split('\n') if line != ""]


if __name__ == '__main__':
    main()
