#!/usr/bin/python3

import os, sys
import os.path
import argparse
import subprocess as pr
from multiprocessing import Process, Queue, JoinableQueue, Lock, Value
import glob
import shlex
import shutil
from . import ArgParsing as ap
from .Supporting import *
from . import Supporting as sp
from . import ProgressBar as pb

#import resource
#import tracemalloc

def main(args=None):
    log(msg="# log notes\n", level="STEP")
    args = ap.parse_phase_snps_arguments(args)
    logArgs(args, 80)

    #tracemalloc.start()

    if args["refvers"] not in ["hg19", "hg38"]:
        raise ValueError(sp.error("The reference genome version of your samples is not \"hg19\" or \"hg38\", please specify one of these two options!\n"))

    rpd = args["refpaneldir"]
    try:
        open( os.path.join(rpd, "1000GP_Phase3", "1000GP_Phase3.sample"), 'r' ) 
    except FileNotFoundError:
        print("Please download the 1000GP reference panel before proceeding")
    panel = os.path.join(rpd, "1000GP_Phase3" )
            
    hg19_path = ""      # path to hg19, 1000GP in hg19 coords, potentially needed for liftover
    chains = ""         # chain files for liftover, chains["hg38_hg19"]=path, chains["hg19_hg38"]=path
    rename_files = ""   # file for renaming chrs with bcftools, rename_files[0] for removing "chr, rename_files[1] for adding "chr" 
    if args["refvers"] == "hg38":
        #dwnld reference panel genome and chain files
        hg19_path = os.path.join(rpd, "hg19_no_chr.fa")
        if args["chrnot"] == "true":
            chains = {"hg38_hg19" : os.path.join(rpd, "hg38ToHg19.chr.chain"), 
                    "hg19_hg38" : os.path.join(rpd, "hg19ToHg38.chr.chain")}
        elif args["chrnot"] == "false":
            chains = {"hg38_hg19" : os.path.join(rpd, "hg38ToHg19.no_chr.chain"), 
                    "hg19_hg38" : os.path.join(rpd, "hg19ToHg38.no_chr.chain")}
        if not os.path.isfile(chains["hg38_hg19"]) or not os.path.isfile(chains["hg19_hg38"]):
            raise ValueError(sp.error("The appropriate liftover chain files could not be located! Please run the PhasePrep.py command that downloads these\n"))

    elif args["refvers"] == "hg19" and args["chrnot"] == "true":
        rename_files = [ os.path.join(rpd, f"rename_chrs{i}.txt") for i in range(1,3) ] 

    # liftover VCFs, phase, liftover again to original coordinates 
    if not os.path.exists(args["outdir"]):
        os.makedirs(args["outdir"])
    
    chromosomes = []
    for chro in args["chromosomes"]:
        if chro.endswith('X') or chro.endswith('Y'):
            sp.log(msg=f'Skipping chromosome {chro} (because it ends with X or Y)\n', level = 'WARN')
        else:
            chromosomes.append(chro)
    
    vcfs = phase(panel, snplist=args["snps"], outdir=args["outdir"], chromosomes=chromosomes, 
                hg19=hg19_path, ref=args["refgenome"], chains=chains, rename=rename_files, refvers=args["refvers"], chrnot=args["chrnot"], 
                num_workers=args["j"], verbose=False) 
    concat_vcf = concat(vcfs, outdir=args["outdir"], chromosomes=chromosomes)


    # read shapeit output, print fraction of phased snps per chromosome
    print_log(path=args["outdir"], chromosomes=chromosomes)
    cleanup(args["outdir"])
            
    """
    current, peak = tracemalloc.get_traced_memory()
    print(f"Current memory usage is {current / 10**6}MB; Peak was {peak / 10**6}MB")
    print(concat_vcf)
    print(resource.getrusage(resource.RUSAGE_SELF))
    print()
    print(resource.getrusage(resource.RUSAGE_CHILDREN))
    """
    
def cleanup(outdir):
    f = []
    # shapeit logs
    [ f.extend( glob.glob(f"shapeit*{ext}"))  for ext in [".log", ".mm", ".hh"]]
    # intermediate files
    exts = ["_phased.vcf.gz", "_phased.vcf.gz.csi", "_filtered.vcf.gz", 
            "_rejected.vcf.gz", "_lifted.vcf.gz", "_lifted.vcf.gz.tbi",
            "_toFilter.vcf.gz", "_toConcat.vcf.gz", "_toConcat.vcf.gz.csi", 
            ".haps", ".sample", "_alignments.snp.strand", "_alignments.snp.strand.exclude"]
    [ f.append( os.path.join(outdir, f"{c}{e}") ) for c in range(1,23) for e in exts ]
    [ os.remove(i) for i in f if os.path.isfile(i) ]

def print_log(path, chromosomes):
    out = open( os.path.join(path, "phased.log"), 'w')
    print("chrom", "phased_snps", "original_snps", "proportion", file=out, sep="\t")
    for c in chromosomes:
        for l in open( os.path.join(path, f"{c}_alignments.log"), 'r' ):
            if "SNPs included" in l:
                snps = int(l.split()[1])
            elif "reference panel sites included" in l:
                phased_snps = int(l.split()[1])
        print(c, phased_snps, snps, float(phased_snps/snps), file=out, sep="\t")

def concat(vcfs, outdir, chromosomes):
    errname = os.path.join(outdir, f"concat.log")
    infiles = ' '.join(vcfs)
    outfile = os.path.join(outdir, f"phased.vcf.gz")
    cmd_bcf = f"bcftools concat --output-type z --output {outfile} {infiles}"
    with open(errname, 'w') as err:
        run = pr.run(cmd_bcf, stdout=err, stderr=err, shell=True, universal_newlines=True)
    if run.returncode != 0:
        raise ValueError(sp.error(f"bcftools concat failed, please check errors in {errname}!"))
    else:
        os.remove(errname)
    return(outfile)

def phase(panel, snplist, outdir, chromosomes, hg19, ref, chains, rename, refvers, chrnot, num_workers, verbose=False):
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
    workers = [ Phaser(tasks, results, progress_bar, panel, outdir, snplist, hg19, ref, chains, rename, refvers, chrnot, verbose) 
                for i in range(min(num_workers, jobs_count)) ]

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

    def __init__(self, task_queue, result_queue, progress_bar, panel, outdir, snplist, hg19, ref, chains, rename, refvers, chrnot, verbose):
        Process.__init__(self)
        self.task_queue = task_queue
        self.result_queue = result_queue
        self.progress_bar = progress_bar
        self.panel = panel 
        self.outdir = outdir
        self.snplist = snplist
        self.hg19 = hg19
        self.ref = ref
        self.chains = chains
        self.rename = rename
        self.refvers = refvers
        self.chrnot = chrnot
        self.verbose = verbose

    def run(self):
        while True:
            next_task = self.task_queue.get() # tuple (snplist[chro], chro)
            if next_task is None:
                # Poison pill means shutdown
                self.task_queue.task_done()
                break

            # begin work
            self.progress_bar.progress(advance=False, msg=f"{self.name} starts on vcf {next_task[0]} for chromosome {next_task[1]}")
            # (1) PREPROCESS
            if self.refvers == "hg19":
                # no need for liftover, just deal with chr annotation
                if self.chrnot == "true":
                    vcf_toFilter = self.change_chr(infile=next_task[0], chromosome=next_task[1], outname="toFilter", rename=self.rename[0])
                else:
                    vcf_toFilter = self.stage_vcfs(infile=next_task[0], chromosome=next_task[1]) # just copy files, vcfs already appropriately formatted
            else:
                # liftover
                vcf_toFilter = self.liftover(infile=next_task[0], chromosome=next_task[1], outname="toFilter", chain=self.chains["hg38_hg19"], refgen=self.hg19, ch="false")

            # (2) FILTERING AND PHASING
            vcf_filtered = self.biallelic(infile=vcf_toFilter, chromosome=next_task[1]) # filter out multi-allelic sites and indels
            vcf_phased = self.shapeit(infile=vcf_filtered, chromosome=next_task[1]) # phase

            # (3) POSTPROCESS
            if self.refvers == "hg19":
                if self.chrnot == "true":
                    vcf_toConcat = self.change_chr(infile=vcf_phased, chromosome=next_task[1], outname="toConcat", rename=self.rename[1])
                    self.index(infile=vcf_toConcat, chromosome=next_task[1]) # re-index with bcftools after renaming
                else:
                    vcf_toConcat = vcf_phased # do nothing; vcfs already in original format
            else:
                vcf_toConcat = self.liftover(infile=vcf_phased, chromosome=next_task[1], outname="toConcat", chain=self.chains["hg19_hg38"], refgen=self.ref, ch=self.chrnot)

            self.progress_bar.progress(advance=True, msg=f"{self.name} ends on vcf {next_task[0]} for chromosome {next_task[1]})")
            self.task_queue.task_done()
            self.result_queue.put(vcf_toConcat)
        return

    def liftover(self, infile, chromosome, outname, chain, refgen, ch):
        errname = os.path.join(self.outdir, f"{chromosome}_picard.log")
        tmpfile = os.path.join(self.outdir, f"{chromosome}_lifted.vcf.gz") # output from picard liftover, to be filtered
        outfile = os.path.join(self.outdir, f"{chromosome}_{outname}.vcf.gz") # filtered with bcftools
        rejfile = os.path.join(self.outdir, f"{chromosome}_rejected.vcf.gz") # required file of SNPs that didn't liftover
        cmd1 = f"picard LiftoverVcf -Xmx4g I={infile} O={tmpfile} CHAIN={chain} R={refgen} REJECT={rejfile}"
        c = chromosome if ch == "false" else f"chr{chromosome}" # need to change "chr" notation depending on liftover direction
        cmd2 = f"bcftools filter --output-type z --regions {c} {tmpfile}" # filter out mapping to other chromosomes/contigs!
        cmd3 = f"bcftools norm --remove-duplicates --output {outfile}" # remove duplicate sites from liftover
        with open(errname, 'w') as err:
            pic = pr.run(cmd1, stdout=err, stderr=err, shell=True, universal_newlines=True)
            filt = pr.Popen(shlex.split(cmd2), stdout=pr.PIPE, stderr=err, universal_newlines=True)
            norm = pr.Popen(shlex.split(cmd3), stdin=filt.stdout, stdout=err, stderr=err, universal_newlines=True)
            codes = map(lambda p : p.wait(), [filt, norm])
        if any(c != 0 for c in codes) or pic.returncode != 0:
            raise ValueError(sp.error(f"Failed to liftover chromosomes with picard on {infile}, please check errors in {errname}!"))
        else:
            os.remove(errname)
        return outfile

    def change_chr(self, infile, chromosome, outname, rename):
        # use bcftools to rename chromosomes
        errname = os.path.join(self.outdir, f"{chromosome}_bcftools.log")
        outfile = os.path.join(self.outdir, f"{chromosome}_{outname}.vcf.gz")
        cmd_bcf = f"bcftools annotate --rename-chrs {rename} --output-type z --output {outfile} {infile}"
        print(cmd_bcf)
        with open(errname, 'w') as err:
            run = pr.run(cmd_bcf, stdout=err, stderr=err, shell=True, universal_newlines=True)
        if run.returncode != 0:
            raise ValueError(sp.error(f"Failed to reannotate chromosomes with bcftools on {infile}, please check errors in {errname}!"))
        else:
            os.remove(errname)
        return outfile

    def stage_vcfs(self, infile, chromosome):
        # copy file, so that phasing takes same input file name regardless of conditions, unzip
        outfile = os.path.join(self.outdir, f"{chromosome}_toFilter.vcf.gz") 
        shutil.copyfile(infile, outfile)
        return outfile

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
        return outfile

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
            codes = map(lambda p : p.returncode, [run2, run3, run4, run5])
            #codes = [run2.returncode, run3.returncode, run4.returncode, run5.returncode] # not collecting error codes from run1 at moment; has nonzero exit status even when it works
        if any(c != 0 for c in codes):
            raise ValueError(sp.error(f"Phasing failed on {infile}, please check errors in {errname}!"))
        else:
            os.remove(errname)
        return os.path.join(self.outdir, f"{chromosome}_phased.vcf.gz")

    def index(self, infile, chromosome): 
        # use bcftools to rename chromosomes
        errname = os.path.join(self.outdir, f"{chromosome}_bcftools.log")
        cmd = f"bcftools index {infile}"
        with open(errname, 'w') as err:
            run = pr.run(cmd, stdout=err, stderr=err, shell=True, universal_newlines=True)
        if run.returncode != 0:
            raise ValueError(sp.error(f"Failed to index {infile}with bcftools, please check errors in {errname}!"))
        else:
            os.remove(errname)
        return 

if __name__ == '__main__':
    main()
