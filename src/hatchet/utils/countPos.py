import multiprocessing as mp
import os
from resource import RUSAGE_SELF
import sys
import subprocess as sp
from datetime import datetime
from collections import Counter, defaultdict
import numpy as np
import argparse
import gzip
import logging
import pysam

from .ArgParsing import parse_count_arguments
from .Supporting import log, logArgs
from hatchet import config

import tracemalloc

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
    
    logger = mp.log_to_stderr(logging.INFO)
    params = zip(np.repeat(chromosomes,  len(bams)), 
        [outdir] * len(bams) * len(chromosomes), 
        [samtools] * len(bams) * len(chromosomes), 
        bams * len(chromosomes),
        names * len(chromosomes), [logger] * len(bams) * len(chromosomes))

    with mp.Pool(np.round(jobs / 2)) as p: # divide by 2 because each worker starts 2 processes
        p.map(process_chromosome_wrapper, params)
        
    mosdepth = sp.run()

def process_chromosome(ch, outdir, samtools, bam, sample_name, logger):
    try:    
        tracemalloc.start()
        outfile = os.path.join(outdir, f'{sample_name}.{ch}.starts')
        if os.path.exists(outfile):
            logger.info(f"Skipping sample {sample_name}-chromosome {ch} (output file exists)")
            return
        
        logger.info(f"{datetime.now()} Sample {sample_name} -- Starting chromosome {ch}")
        
        # Get start positions
        st = sp.Popen((samtools, 'view', bam,  ch), stdout=sp.PIPE)
        cut = sp.Popen(('cut', '-f', '4'), stdin=st.stdout, stdout=open(outfile, 'w'))
        st.wait()
        cut.wait()
        if st.returncode != 0:
            raise ValueError("samtools subprocess returned nonzero value: {}".format(st.returncode))
        if cut.returncode != 0:
            raise ValueError("cut subprocess returned nonzero value: {}".format(cut.returncode))

        logger.info(f"{datetime.now()} Sample {sample_name} -- Done chromosome {ch}")
        current, peak = tracemalloc.get_traced_memory()
        logger.info(f"Current memory usage is {current / 10**6}MB; Peak was {peak / 10**6}MB")
        
        tracemalloc.stop()
        
    except Exception as e:
        logger.exception(e)
        pass
    
def process_chromosome_wrapper(param):
    process_chromosome(*param)

if __name__ == '__main__':
    main()

