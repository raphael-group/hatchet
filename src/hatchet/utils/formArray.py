from multiprocessing import Pool
import os
from datetime import datetime
import numpy as np
import pandas as pd
import gzip
import subprocess
import traceback
import tracemalloc

from .ArgParsing import parse_array_arguments
from . import Supporting as sp


def main(args=None):    
    sp.log(msg="# Parsing and checking input arguments\n", level="STEP")
    args = parse_array_arguments(args)
    sp.logArgs(args, 80)
    
    stem = args['stem']
    threads = args['processes']
    chromosomes = args['chromosomes']
    outstem = args['outstem']

    all_names = args['sample_names']
    use_chr = args['use_chr']
    compressed = args['compressed']
    centromeres = args['centromeres']

    n_workers = min(len(chromosomes), threads)

    # Read in centromere locations table
    centromeres = pd.read_table(centromeres, header = None, names = ['CHR', 'START', 'END', 'NAME', 'gieStain'])
    chr2centro = {}
    for ch in centromeres.CHR.unique():
        my_df = centromeres[centromeres.CHR == ch]
        assert (my_df.gieStain == 'acen').all()
        # Each centromere should consist of 2 adjacent segments
        assert len(my_df == 2)
        assert my_df.START.max() == my_df.END.min()
        if use_chr:
            if ch.startswith('chr'):
                chr2centro[ch] = my_df.START.min(), my_df.END.max()
            else:
                chr2centro['chr' + ch] = my_df.START.min(), my_df.END.max()
        else:
            if ch.startswith('chr'):
                chr2centro[ch[3:]] = my_df.START.min(), my_df.END.max()
            else:
                chr2centro[ch] = my_df.START.min(), my_df.END.max()
    
    for ch in chromosomes:
        if ch not in chr2centro:
            raise ValueError(sp.error(f"Chromosome {ch} not found in centromeres file. Inspect file provided as -C argument."))

    # form parameters for each worker
    params = [(stem, all_names, ch, outstem + f'.{ch}', 
               chr2centro[ch][0], chr2centro[ch][1], compressed)
              for ch in chromosomes]
    
    # dispatch workers
    with Pool(n_workers) as p:
        p.map(run_chromosome_wrapper, params)
        
    np.savetxt(outstem + '.samples', all_names, fmt = '%s')


def read_snps(baf_file, ch, all_names):
    """
    Read and validate SNP data for this patient (TSV table output from HATCHet deBAF.py).
    """
    all_names = all_names[1:] # remove normal sample -- not looking for SNP counts from normal

    # Read in HATCHet BAF table
    all_snps = pd.read_table(baf_file, names = ['CHR', 'POS', 'SAMPLE', 'ALT', 'REF'], 
                             dtype = {'CHR':object, 'POS':np.uint32, 'SAMPLE':object, 
                                      'ALT':np.uint32, 'REF':np.uint32})
    
    # Keep only SNPs on this chromosome
    snps = all_snps[all_snps.CHR == ch].sort_values(by=['POS', 'SAMPLE'])
    snps = snps.reset_index(drop = True)

    if len(snps) == 0:
        raise ValueError(sp.error(f"Chromosome {ch} not found in SNPs file (chromosomes in file: {all_snps.CHR.unique()})"))
    
    n_samples = len(all_names)
    if n_samples != len(snps.SAMPLE.unique()):
        raise ValueError(sp.error(
            f"Expected {n_samples} samples, found {len(snps.SAMPLE.unique())} samples in SNPs file."
        ))

    if set(all_names) != set(snps.SAMPLE.unique()):
        raise ValueError(sp.error(
            f"Expected sample names did not match sample names in SNPs file.\n\
                Expected: {sorted(all_names)}\n  Found:{sorted(snps.SAMPLE.unique())}"
        ))

    # Add total counts column
    snpsv = snps.copy()
    snpsv['TOTAL'] = snpsv.ALT + snpsv.REF
    
    # Create counts array and find SNPs that are not present in all samples
    snp_counts = snpsv.pivot(index = 'POS', columns = 'SAMPLE', values = 'TOTAL')
    missing_pos = snp_counts.isna().any(axis = 1)

    # Remove SNPs that are absent in any sample
    snp_counts = snp_counts.dropna(axis = 0)
    snpsv = snpsv[~snpsv.POS.isin(missing_pos[missing_pos].index)]
        
    # Pivot table for dataframe should match counts array and have no missing entries
    check_pivot = snpsv.pivot(index = 'POS', columns = 'SAMPLE', values = 'TOTAL')
    assert np.array_equal(check_pivot, snp_counts), "SNP file reading failed"
    assert not np.any(check_pivot.isna()), "SNP file reading failed"
    assert np.array_equal(all_names, list(snp_counts.columns)) # make sure that sample order is the same
    
    return np.array(snp_counts.index), np.array(snp_counts), snpsv

def form_counts_array(starts_files, perpos_files, thresholds, chromosome, chunksize = 1e5, tabix = 'tabix'):
    """
        NOTE: Assumes that starts_files[i] corresponds to the same sample as perpos_files[i]
        Parameters:
            starts_files: list of <sample>.<chromosome>.starts.gz files each containing a list of start positions
            perpos_files: list of <sample>.per-base.bed.gz files containing per-position coverage from mosdepth
            starts: list of potential bin start positions (thresholds between SNPs)
            chromosome: chromosome to extract read counts for
        
        Returns: <n> x <2d> np.ndarray 
            entry [i, 2j] contains the number of reads starting in (starts[i], starts[i + 1]) in sample j
            entry [i, 2j + 1] contains the number of reads covering position starts[i] in sample j
    """
    arr = np.zeros((thresholds.shape[0] + 1, len(starts_files) * 2)) # add one for the end of the chromosome

    for i in range(len(starts_files)):
        # populate starts in even entries
        fname = starts_files[i]
        next_idx = 1
        if fname.endswith('.gz'):
            f = gzip.open(fname)
        else:
            f = open(fname)
        
        for line in f:
            pos = int(line)
            while next_idx < len(thresholds) and pos > thresholds[next_idx]:
                next_idx += 1
                
            if next_idx == len(thresholds):
                arr[next_idx - 1, i * 2] += 1
            elif not (pos == thresholds[next_idx - 1] or pos == thresholds[next_idx]):
                assert pos > thresholds[next_idx - 1] and pos < thresholds[next_idx], (next_idx, pos)
                arr[next_idx - 1, i * 2] += 1
                
        f.close()
            
    for i in range(len(perpos_files)):
        # populate threshold coverage in odd entries
        fname = perpos_files[i]
        #print(datetime.now(), "Reading {}".format(fname))
        
        chr_sample_file = os.path.join(fname[:-3] + '.' + chromosome)
        
        if not os.path.exists(chr_sample_file):
            with open(chr_sample_file, 'w') as f:
                subprocess.run([tabix, fname, chromosome], stdout = f)
                
        with open(chr_sample_file, 'r') as records:        
            idx = 0
            last_record = None
            for line in records:
                tokens = line.split()
                if len(tokens) == 4:
                    start = int(tokens[1])
                    end = int(tokens[2])
                    nreads = int(tokens[3])

                    while idx < len(thresholds) and thresholds[idx] - 1 < end:
                        assert thresholds[idx] - 1 >= start
                        arr[idx, (2 * i) + 1] = nreads
                        idx += 1
                    last_record = line
                        
            if i == 0:
                # identify the (effective) chromosome end as the last well-formed record
                assert idx == len(thresholds)
                _, _, chr_end, end_reads = last_record.split()
                chr_end = int(chr_end)
                end_reads = int(end_reads)

                assert chr_end > thresholds[-1] 
                assert len(thresholds) == len(arr) - 1

                # add the chromosome end to thresholds
                thresholds = np.concatenate([thresholds, [chr_end]])

                # count the number of reads covering the chromosome end
                arr[idx, 0] = end_reads
        
        if os.path.exists(chr_sample_file):
            os.remove(chr_sample_file)

    return arr, thresholds

def run_chromosome(stem, all_names, chromosome, outstem, centromere_start, centromere_end, compressed):
    """
    Perform adaptive binning and infer BAFs to produce a HATCHet BB file for a single chromosome.
    """

    try:
        tracemalloc.start()

        totals_out = outstem + '.total'
        thresholds_out = outstem + '.thresholds'

        if os.path.exists(totals_out) and os.path.exists(thresholds_out):
            sp.log(msg=f"Output files already exist, skipping chromosome {chromosome}\n", level = "INFO")
            return

        sp.log(msg=f"Loading chromosome {chromosome}\n", level = "INFO")
        # Per-position coverage bed files for each sample
        perpos_files = [os.path.join(stem, 'counts', name + '.per-base.bed.gz') for name in all_names]
        
        # Identify the start-positions files for this chromosome
        starts_files = []
        for name in all_names:
            if compressed:
                starts_files.append(os.path.join(stem, 'counts', name + '.' + chromosome + '.starts.gz'))    
            else:
                starts_files.append(os.path.join(stem, 'counts', name + '.' + chromosome + '.starts'))    

        sp.log(msg=f"Reading SNPs file for chromosome {chromosome}\n", level = "INFO")
        # Load SNP positions and counts for this chromosome
        positions, snp_counts, snpsv = read_snps(os.path.join(stem, 'baf', 'bulk.1bed'), chromosome, all_names)
        
        thresholds = np.trunc(np.vstack([positions[:-1], positions[1:]]).mean(axis = 0)).astype(np.uint32)
        last_idx_p = np.argwhere(thresholds > centromere_start)[0][0]
        first_idx_q = np.argwhere(thresholds > centromere_end)[0][0]
        all_thresholds = np.concatenate([[1], thresholds[:last_idx_p], [centromere_start], 
                                        [centromere_end], thresholds[first_idx_q:]])
        
        
        sp.log(msg=f"Loading counts for chromosome {chromosome}\n", level = "INFO")
        # Load read count arrays from file (also adds end of chromosome as a threshold)
        total_counts, complete_thresholds = form_counts_array(starts_files, perpos_files, all_thresholds, chromosome) 

        np.savetxt(totals_out, total_counts, fmt = '%d')
        np.savetxt(thresholds_out, complete_thresholds, fmt = '%d')
        
        sp.log(msg=f"Done chromosome {chromosome}\n", level ="INFO")
        
        current, peak = tracemalloc.get_traced_memory()
        sp.log(msg=f"Chr {chromosome} -- Current memory usage is {int(current / 10**6)}MB; Peak was {int(peak / 10**6)}MB\n",
            level = "INFO")
        tracemalloc.stop()
    except Exception as e: 
        print(f"Error in chromosome {chromosome}:")
        print(e)
        traceback.print_exc()
        raise e
    
def run_chromosome_wrapper(param):
    run_chromosome(*param)

if __name__ == '__main__':
    main()