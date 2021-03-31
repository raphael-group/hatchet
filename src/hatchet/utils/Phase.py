#!/usr/bin/python3

import os, sys
import os.path
import argparse
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
    concat_vcf = concat(vcfs, outdir=args["outputphase"])

    # read shapeit output, get fraction of phased snps
    out = open( os.path.join(args["outputphase"], "phased.log"), 'w')
    print("chrom", "phased_snps", "original_snps", "proportion", file=out, sep="\t")
    for c in args["chromosomes"]:
        fn = os.path.join(args["outputphase"], f"{c}_alignments.log")
        for l in open (fn, 'r'):
            if "SNPs included" in l:
                snps = int(l.split()[1])
            elif "reference panel sites included" in l:
                phased_snps = int(l.split()[1])
        print(c, phased_snps, snps, float(phased_snps/snps), file=out, sep="\t")
    out.close()
            
    print(concat_vcf)

def concat(vcfs, outdir):
    errname = os.path.join(outdir, f"concat.log")
    infiles = ' '.join(vcfs)
    outfile = os.path.join(outdir, f"phased.vcf.gz")
    cmd_bcf = f"bcftools concat --output-type z --output {outfile} "
    cmd_bcf += f"{infiles}"
    with open(errname, 'w') as err:
        run = pr.run(cmd_bcf, stdout=err, stderr=err, shell=True, universal_newlines=True)
    if run.returncode != 0:
        raise ValueError(sp.error(f"bcftools concat failed, please check errors in {errname}!"))
    else:
        os.remove(errname)
    return(outfile)

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

            self.progress_bar.progress(advance=False, msg=f"{self.name} starts on vcf {next_task[0]} for chromosome {next_task[1]}")
            bi = self.biallelic(infile=next_task[0], chromosome=next_task[1])
            phased = self.shapeit(infile=next_task[0], chromosome=next_task[1])
            self.progress_bar.progress(advance=True, msg=f"{self.name} ends on vcf {next_task[0]} for chromosome {next_task[1]})")
            self.task_queue.task_done()
            self.result_queue.put(phased)
        return

    def biallelic(self, infile, chromosome):
        # use bcftools to discard multi-allelic sites and indels
        errname = os.path.join(self.outdir, f"{chromosome}_bcftools.log")
        outfile = os.path.join(self.outdir, f"{chromosome}_filtered.vcf.gz")
        cmd_bcf = f"bcftools view --max-alleles 2 --exclude-types indels --output-type z "
        cmd_bcf += f"--output-file {outfile} {infile}"
        with open(errname, 'w') as err:
            run = pr.run(cmd_bcf, stdout=err, stderr=err, shell=True, universal_newlines=True)
        if run.returncode != 0:
            raise ValueError(sp.error(f"Biallelic sites filtering failed on {infile}, please check errors in {errname}!"))
        else:
            os.remove(errname)
        return os.path.join(outfile)

    def shapeit(self, infile, chromosome): 
        # use shapeit with reference panel to phase vcf files
        errname = os.path.join(self.outdir, f"{chromosome}_shapeit.log")

        # define params used across shapeit functions
        inmap = f"{self.panel}/genetic_map_chr{chromosome}_combined_b37.txt"
        inref = f"{self.panel}/1000GP_Phase3_chr{chromosome}.hap.gz "
        inref += f"{self.panel}/1000GP_Phase3_chr{chromosome}.legend.gz "
        inref += f"{self.panel}/1000GP_Phase3.sample"

        # check data with shapeit -check; get list of sites to exclude, such as sites in target VCF that are not present in reference panel
        cmd1 = f"shapeit -check --input-vcf {self.outdir}/{chromosome}_filtered.vcf.gz "
        cmd1 += f"--input-map {inmap} "
        cmd1 += f"--input-ref {inref} "
        cmd1 += f"--output-log {self.outdir}/{chromosome}_alignments"
        cmd1 += "\n"
        # phase
        cmd2 = f"shapeit --input-vcf {self.outdir}/{chromosome}_filtered.vcf.gz "
        cmd2 += f"--input-map {inmap} "
        cmd2 += f"--input-ref {inref} "
        cmd2 += f"--exclude-snp {self.outdir}/{chromosome}_alignments.snp.strand.exclude "
        cmd2 += f"--output-max {self.outdir}/{chromosome}.haps {self.outdir}/{chromosome}.sample "
        cmd2 += f"--chrX --no-mcmc "
        # convert output file to vcf
        cmd3 = f"shapeit -convert --input-haps {self.outdir}/{chromosome} "
        cmd3 += f"--output-vcf {self.outdir}/{chromosome}_phased.vcf"
        # compress vcf 
        cmd4 = f"bgzip {self.outdir}/{chromosome}_phased.vcf"
        # index vcf
        cmd5 = f"bcftools index {self.outdir}/{chromosome}_phased.vcf.gz"
        with open(errname, 'w') as err:
            run1 = pr.run(cmd1, stdout=err, stderr=err, shell=True, universal_newlines=True)
            run2 = pr.run(cmd2, stdout=err, stderr=err, shell=True, universal_newlines=True)
            run3 = pr.run(cmd3, stdout=err, stderr=err, shell=True, universal_newlines=True)
            run4 = pr.run(cmd4, stdout=err, stderr=err, shell=True, universal_newlines=True)
            run5 = pr.run(cmd5, stdout=err, stderr=err, shell=True, universal_newlines=True)
            codes = [run2.returncode, run3.returncode, run4.returncode, run5.returncode] # not collecting error codes from run1 at moment; has nonzero exit status even when it works
        if any(c != 0 for c in codes):
            raise ValueError(sp.error(f"Phasing failed on {infile}, please check errors in {errname}!"))
        else:
            os.remove(errname)
            os.remove( os.path.join(self.outdir, f"{chromosome}.haps") )
            os.remove( os.path.join(self.outdir, f"{chromosome}.sample") )
            os.remove( os.path.join(self.outdir, f"{chromosome}_alignments.snp.strand") )
        return os.path.join(self.outdir, f"{chromosome}_phased.vcf.gz")

if __name__ == '__main__':
    main()
