from multiprocessing import Pool
import os
import sys
import subprocess as sp
from datetime import datetime
from collections import Counter, defaultdict
import numpy as np
import argparse
import gzip

from .ArgParsing import parse_count_arguments
from .Supporting import log, logArgs
from hatchet import config

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
    
    params = zip(np.repeat(chromosomes,  len(bams)), 
        [outdir] * len(bams) * len(chromosomes), 
        [samtools] * len(bams) * len(chromosomes), 
        bams * len(chromosomes),
        names * len(chromosomes))

    with Pool(jobs) as p:
        # TODO: incorporate ProgressBar, logging
        p.map(process_chromosome_wrapper, params)

def process_chromosome(ch, outdir, samtools, bam, sample_name):
    # TODO: replace print statements with sp.log
    
    outfile = os.path.join(outdir, f'{sample_name}.{ch}.gz')
    if os.path.exists(outfile):
        print(f"Skipping sample {sample_name}-chromosome {ch} (output file exists)")
        return
    
    print(datetime.now(), f"Sample {sample_name} -- Starting chromosome {ch}")
    
    # Get start positions
    cmd = " ".join([samtools, 'view', bam,  ch, '|',  'awk', "-F'\t'", "'{print $4}'"])
    #print(cmd)
    result = sp.run(cmd, shell = True, capture_output = True)
    result.check_returncode()


    starts = Counter(result.stdout.decode('utf-8').strip().split('\n'))
    
    # Count reads per position
    read_length = 151

    startpos = min([int(a) for a in starts.keys()])
    endpos = max([int(a) for a in starts.keys()])

    cov_per_pos = np.zeros(endpos - startpos + read_length, dtype = np.int32)

    for pos, n in starts.items():
        pos = int(pos)
        cov_per_pos[pos - startpos:pos - startpos + read_length - 1] += n


    starts_dict = {int(k):v for k,v in starts.items()}
    with gzip.open(outfile, 'wb') as f:
        f.write(f'{startpos}\n'.encode())
        [f.write(f'{int(starts_dict[i + startpos]) if (i + startpos) in starts_dict else 0},{cov_per_pos[i]}\n'.encode())
         for i in range(len(cov_per_pos))]

    print(datetime.now(), f"Sample {sample_name} -- Done chromosome {ch}")
    
def process_chromosome_wrapper(param):
    process_chromosome(*param)

if __name__ == '__main__':
    main()

