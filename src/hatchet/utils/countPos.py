import multiprocessing as mp
import os
import subprocess as sp
import numpy as np
from shutil import which
import sys

from .ArgParsing import parse_count_arguments
from .Supporting import log, logArgs, error

def main(args=None):
    log(msg="# Parsing and checking input arguments\n", level="STEP")
    args = parse_count_arguments(args)
    logArgs(args, 80)
    
    bams = args["bams"]
    names = args["names"]
    chromosomes = args["chromosomes"]
    samtools = args["samtools"]
    jobs = args["j"]
    outdir = args["outdir"]
    compression = args['compression']
    
    if which("mosdepth") is None:
        raise ValueError(error("The 'mosdepth' executable was not found on PATH."))

    
    params = zip(np.repeat(chromosomes,  len(bams)), 
        [outdir] * len(bams) * len(chromosomes), 
        [samtools] * len(bams) * len(chromosomes), 
        bams * len(chromosomes),
        names * len(chromosomes),
        [compression] * len(bams) * len(chromosomes))

    n_workers_samtools = min(int(np.round(jobs / 2)), len(bams) * len(chromosomes))
    with mp.Pool(n_workers_samtools) as p: # divide by 2 because each worker starts 2 processes
        p.map(process_chromosome_wrapper, params)
       
    n_workers_mosdepth = min(jobs, len(bams))

    # compute number of decompression threads to use for each call to mosdepth
    if jobs > len(bams):
        base = min(int(jobs / n_workers_mosdepth), 4)
        threads_per_worker = [base] * n_workers_mosdepth

        if base < 4:
            remainder = jobs % n_workers_mosdepth
            i = 0
            while remainder > 0:
                threads_per_worker[i] += 1
                i += 1
                remainder -= 1
    else:
        threads_per_worker = [1] * n_workers_mosdepth

    mosdepth_params = [(outdir, names[i], bams[i], threads_per_worker[i]) for i in range(len(bams))]
    with mp.Pool(n_workers_mosdepth) as p:
        p.map(mosdepth_wrapper, mosdepth_params)
        
def mosdepth_wrapper(params):
    run_mosdepth(*params)
    
def run_mosdepth(outdir, sample_name, bam, threads):
    try:
        last_file = os.path.join(outdir, sample_name + '.per-base.bed.gz.csi')
        if os.path.exists(last_file):
            log(f"Skipping mosdepth on sample {sample_name} (output file {last_file} exists)\n", level = "STEP")
            return
        
        log(f"Starting mosdepth on sample {sample_name} with {threads} threads\n", level = "STEP")
        sys.stderr.flush()
        
        mosdepth = sp.run(['mosdepth', '-t', str(threads), os.path.join(outdir, sample_name), bam])
        mosdepth.check_returncode()
        log(f"Done mosdepth on sample {sample_name}\n", level = "STEP")

    except Exception as e:
        log("Exception in countPos: {}\n".format(e), level = "ERROR")
        raise e  

def process_chromosome(ch, outdir, samtools, bam, sample_name, compression):
    try:    
        outfile = os.path.join(outdir, f'{sample_name}.{ch}.starts')
        if os.path.exists(outfile):
            log(f"Skipping sample {sample_name} chromosome {ch} (output file exists)\n", level = "STEP")
            return
        
        log(f"Sample {sample_name} -- Starting chromosome {ch}\n", level = "STEP")
        sys.stderr.flush()
        
        # Get start positions
        st = sp.Popen((samtools, 'view', bam,  ch), stdout=sp.PIPE)
        cut = sp.Popen(('cut', '-f', '4'), stdin=st.stdout, stdout=sp.PIPE)
        gzip = sp.Popen(('gzip', '-{}'.format(compression)),  stdin = cut.stdout, stdout=open(outfile + ".gz", 'w'))
        st.wait()
        cut.wait()
        gzip.wait()
        if st.returncode != 0:
            raise ValueError("samtools subprocess returned nonzero value: {}".format(st.returncode))
        if cut.returncode != 0:
            raise ValueError("cut subprocess returned nonzero value: {}".format(cut.returncode))

        log(f"Sample {sample_name} -- Done chromosome {ch}\n", level = "STEP")
        
    except Exception as e:
        log("Exception in countPos: {}\n".format(e), level = "ERROR")
        raise e
    
def process_chromosome_wrapper(param):
    process_chromosome(*param)

if __name__ == '__main__':
    main()

