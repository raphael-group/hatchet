#!/usr/bin/python3

import os, sys
import os.path
import argparse
import subprocess as pr
from multiprocessing import Process, Queue, JoinableQueue, Lock, Value
import requests
import tarfile
import glob
import gzip
import shutil
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

    if args["refvers"] not in ["hg19", "hg38"]:
        raise ValueError(sp.error("The reference genome version of your samples is not \"hg19\" or \"hg38\", please specify one of these two options!\n"))

    # download reference panel, prepare files for liftover
    hg19_path = ""      # path to hg19, 1000GP in hg19 coords, potentially needed for liftover
    chains = ""         # chain files for liftover, chains["hg38_hg19"]=path, chains["hg19_hg38"]=path
    rename_files = ""   # file for renaming chrs with bcftools, rename_files[0] for removing "chr, rename_files[1] for adding "chr" 
    if args["refpanel"] == "1000GP_Phase3":
        # download 1000GP ref panel
        if not os.path.isfile( os.path.join(args["outdir"], "1000GP_Phase3", "1000GP_Phase3.sample") ):
            dwnld_1kgp(path = args["outdir"])
        panel = os.path.join( args["outdir"], "1000GP_Phase3" )
        # download necessary liftover files; 1000GP in hg19 coordinates
        if args["refvers"] == "hg38":
            #dwnld reference panel genome and chain files
            hg19_path = dwnld_refpanel_genome(path=args["outdir"])
            chains = dwnld_chains(path=args["outdir"], refpanel_chr="false", sample_chr=args["chrnot"] )
        elif args["refvers"] == "hg19" and args["chrnot"] == "true":
            rename_files = mk_rename_file(path = args["outdir"]) 
    else:
        raise ValueError(sp.error("Currently, only the 1000 genome panel aligned to GRCh37 without \"chr\" prefix is supported\n")) 
        #panel = args["refpanel"]       # for future: include custom panel




    # DO THIS IN ARG PARSING FUNCTION????
    if args["chrnot"] == "true":
        chromosomes = [ i.replace("chr","") for i in args["chromosomes"] ] # rename chromosomes; used to locate ref panel files!
        snplist = {k.replace("chr","") : v for k,v in args["snps"].items()} # keeps values -> filenames don't change
    else:
        chromosomes = args["chromosomes"]
        snplist = args["snps"]





    # liftover VCFs, phase, liftover again to original coordinates 
    vcfs = phase(panel, snplist=snplist, outdir=args["outdir"], chromosomes=chromosomes, 
                hg19=hg19_path, chains=chains, rename=rename_files, refvers=args["refvers"], chrnot=args["chrnot"], 
                num_workers=args["j"], verbose=False) 
    concat_vcf = concat(vcfs, outdir=args["outdir"], chromosomes=chromosomes)

    # read shapeit output, print fraction of phased snps per chromosome
    print_log(path=args["outdir"], chromosomes = chromosomes)
    cleanup(args["outdir"])
            
    print(concat_vcf)

def cleanup(outdir):
    f = []
    # shapeit logs
    [ f.extend( glob.glob(f"shapeit*{ext}"))  for ext in [".log", ".mm", ".hh"]]
    # intermediate files
    exts = ["_phased.vcf.gz", "_phased.vcf.gz.csi", "_filtered.vcf.gz", 
            "_toFilter.vcf.gz", "_toConcat.vcf.gz", "_toConcat.vcf.gz.csi", 
            ".haps", ".sample", "_alignments.snp.strand", "_alignments.snp.strand.exclude"]
    [ f.append( os.path.join(outdir, f"{c}{e}") ) for c in range(1,23) for e in exts]
    [os.remove(i) for i in f]

def dwnld_chains(path, refpanel_chr, sample_chr):
    # order of urls important! [0] chain for sample -> ref panel, [1] chain for ref panel -> sample
    urls = ("https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz", "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz")
    paths = [os.path.join(path, os.path.basename(i)) for i in urls] # paths for url downloads
    for i, url in enumerate(urls):
        r = requests.get(url, allow_redirects=True)
        open(paths[i], 'wb').write(r.content)
    c1 = mod_chain(paths[0], refpanel_chr, sample_chr, refpanel_index=7, sample_index=2) # modify chr notation of hg38ToHg19, ref panel chr in 7th field, sample chr in 2nd field
    c2 = mod_chain(paths[1], refpanel_chr, sample_chr, refpanel_index=2, sample_index=7) # modify chr notation of hg19ToHg38, ref panel chr in 2th field, sample chr in 7nd field
    [ os.remove(p) for p in paths ] # remove original chain files
    return {"hg38_hg19" : c1, "hg19_hg38" : c2}

def mod_chain(infile, refpanel_chr, sample_chr, refpanel_index, sample_index):
    name = infile.strip(".gz").replace("over", "renamed")
    with open(name, 'w') as new:
        with gzip.open(infile, 'rt') as f:
            for l in f:
                if l.startswith("chain"):
                    l = l.split()
                    if refpanel_chr == "false": l[refpanel_index] = l[refpanel_index].replace("chr","") 
                    if sample_chr == "false": l[sample_index] = l[sample_index].replace("chr","") 
                    new.write(" ".join(l) + "\n")
                else:
                    new.write(l)
    return name

def dwnld_refpanel_genome(path):
    url = "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz"
    out = os.path.join(path, "hg19.fa.gz")
    # download hg19
    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        with open(out, 'wb') as f:                                                                                                                    
            for chunk in r.iter_content(chunk_size=8192):                                                                                             
                f.write(chunk)      
    # change chr notation
    newref = os.path.join(path, "hg19_renamed.fa")
    with open(newref, 'w') as new:
        with gzip.open(out, 'rt') as f:
            [ new.write(l.replace("chr","")) if l.startswith(">") else new.write(l) for l in f ]
    os.remove(out)
    return newref

def mk_rename_file(path):
    # makes rename_chrs1.txt for removing "chr", rename_chrs2.txt for adding "chr"
    names = [ os.path.join(path, f"rename_chrs{i}.txt") for i in range(1,3) ]
    for i, n in enumerate(names):
        with open(n, 'w') as f:
            for j in range(1,23):
                f.write(f"chr{j} {j}\n") if i == 0 else f.write(f"{j} chr{j}\n") 
    return names

def dwnld_1kgp(path):
    url = "https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.tgz"
    out = os.path.join(path + "1000GP_Phase3.tgz")
    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        with open(out, 'wb') as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)
        #f.write(r.content)
    # extract tar file, remove
    t = tarfile.open(out)
    t.extractall(path)
    t.close()
    os.remove(out)

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

def phase(panel, snplist, outdir, chromosomes, hg19, chains, rename, refvers, chrnot, num_workers, verbose=False):
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
    workers = [ Phaser(tasks, results, progress_bar, panel, outdir, snplist, hg19, chains, rename, refvers, chrnot, verbose) 
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

    def __init__(self, task_queue, result_queue, progress_bar, panel, outdir, snplist, hg19, chains, rename, refvers, chrnot, verbose):
        Process.__init__(self)
        self.task_queue = task_queue
        self.result_queue = result_queue
        self.progress_bar = progress_bar
        self.panel = panel 
        self.outdir = outdir
        self.snplist = snplist
        self.hg19 = hg19
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
                    vcf_toFilter = self.stage_vcfs(infile=next_task[0], chromosome=next_task[1])
            else:
                print("Do nothing for now!")
                # liftover

            # (2) FILTERING AND PHASING
            vcf_toPhase = self.biallelic(infile=vcf_toFilter, chromosome=next_task[1]) # filter out multi-allelic sites and indels
            vcf_phased = self.shapeit(infile=vcf_toPhase, chromosome=next_task[1]) # phase

            # (3) POSTPROCESS
            if self.refvers == "hg19":
                if self.chrnot == "true":
                    vcf_toConcat = self.change_chr(infile=vcf_phased, chromosome=next_task[1], outname="toConcat", rename=self.rename[1])
                    self.index(infile=vcf_toConcat, chromosome=next_task[1]) # re-index with bcftools after renaming
                else:
                    vcf_toConcat = vcf_phased
            else:
                print("Do nothing for now!")

            self.progress_bar.progress(advance=True, msg=f"{self.name} ends on vcf {next_task[0]} for chromosome {next_task[1]})")
            self.task_queue.task_done()
            #self.result_queue.put(phased)
            self.result_queue.put(vcf_toConcat)
        return

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
            codes = [run2.returncode, run3.returncode, run4.returncode, run5.returncode] # not collecting error codes from run1 at moment; has nonzero exit status even when it works
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
