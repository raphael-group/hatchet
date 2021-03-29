#!/usr/bin/python3

import os, sys
import os.path
import argparse
import shlex
import subprocess as pr
from multiprocessing import Process, Queue, JoinableQueue, Lock, Value
import requests
import tarfile

from . import ArgParsing as ap
from .Supporting import *
from . import Supporting as sp
from . import ProgressBar as pb



def main(args=None):
    log(msg="# log notes\n", level="STEP")
    args = ap.parse_phase_arguments(args)
    logArgs(args, 80)

    for i in args:
        print(i, args[i])

    if args["refpanel"] == "1000GP_Phase3":
        """
        url = "https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.tgz"
        out = args["outputphase"] + "/1000GP_Phase3.tgz"
        with requests.get(url, stream=True) as r:
            r.raise_for_status()
            with open(out, 'wb') as f:
                for chunk in r.iter_content(chunk_size=8192):
                    f.write(chunk)
            f.write(r.content)
        # extract tar file, remove
        t = tarfile.open(out)
        t.extractall(args["outputphase"])
        t.close()
        os.remove(out)
        """
        panel = '{}/1000GP_Phase3'.format(args["outputphase"])
    else:
        panel = args["refpanel"]

    vcfs = phase(panel, snplist=args["snps"], outdir=args["outputphase"], chromosomes=args["chromosomes"], num_workers=args["j"], verbose=False) 

    print("SORTED RESULTS!")
    print(vcfs)

def phase(panel, snplist, outdir, chromosomes, num_workers, verbose=False):
    # Define a Lock and a shared value for log printing through ProgressBar
    err_lock = Lock()
    counter = Value('i', 0)
    progress_bar = pb.ProgressBar(total=len(chromosomes), length=40, lock=err_lock, counter=counter, verbose=verbose)

    # Establish communication queues
    tasks = JoinableQueue()
    results = Queue()

    # Enqueue jobs
    jobs_count = 0
    for chro in chromosomes:
        tasks.put( (snplist[chro], chro) )
        jobs_count += 1

    # Setting up the workers
    workers = [Phaser(tasks, results, progress_bar, panel, outdir, snplist, verbose) for i in range(min(num_workers, jobs_count))]

    # Add a poison pill for each worker
    for i in range(len(workers)):
        tasks.put(None)

    # Start the workers
    for w in workers:
        w.start()

    # Wait for all of the tasks to finish
    tasks.join()

    # Get the results
    sorted_results = sorted([results.get() for i in range(jobs_count)])

    # Close Queues
    tasks.close()
    results.close()

    # Ensure each worker terminates
    for w in workers:
        w.terminate()
        w.join()

    return sorted_results

class Phaser(Process):

    def __init__(self, task_queue, result_queue, progress_bar, panel, outdir, snplist, verbose):
        Process.__init__(self)
        self.task_queue = task_queue
        self.result_queue = result_queue
        self.progress_bar = progress_bar
        self.panel = panel 
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

            #self.progress_bar.progress(advance=False, msg="{} starts on {} for {})".format(self.name, next_task[1], next_task[2]))
            self.progress_bar.progress(advance=False, msg="{} starts on vcf {} for chromosome {}".format(self.name, next_task[0], next_task[1]))
            phased = self.phasevcfs(vcffile=next_task[0], chromosome=next_task[1])
            self.progress_bar.progress(advance=True, msg="{} ends on vcf {} for chromosome {})".format(self.name, next_task[0], next_task[1]))
            self.task_queue.task_done()
            self.result_queue.put(phased)
        return

    def phasevcfs(self, vcffile, chromosome):
        errname = os.path.join(self.outdir, "{}_phase.log".format(chromosome))
        cmd_bcf = "bcftools view --max-alleles 2 --exclude-types indels --output-type z "
        cmd_bcf += "--output-file {} {}".format( os.path.join(self.outdir, "{}_filtered.vcf.gz".format(chromosome)), self.snplist[chromosome])

        with open(errname, 'w') as err:
            biallelic = pr.Popen(shlex.split(cmd_bcf), stdout=pr.PIPE, stderr=err, universal_newlines=True)
            # biallilic in list to demonstrate you can pass multiple
            codes = map(lambda p : p.wait(), [biallelic])
        if any(c != 0 for c in codes):
            raise ValueError(sp.error('Phasing failed on {}, please check errors in {}!').format(self.snplist[chromosome], errname))
        else:
            os.remove(errname)
        return os.path.join(self.outdir, '{}_filtered.vcf.gz'.format(chromosome))



if __name__ == '__main__':
    main()
