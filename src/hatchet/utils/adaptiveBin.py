from multiprocessing import Pool
import os
from datetime import datetime
import numpy as np
import pandas as pd
from scipy.stats import binom
import gzip
import subprocess
import traceback

import tracemalloc
from .ArgParsing import parse_abin_arguments
from . import Supporting as sp


def main(args=None):
    sp.log(msg="# Parsing and checking input arguments\n", level="STEP")
    args = parse_abin_arguments(args)
    sp.logArgs(args, 80)
    
    stem = args['stem']
    threads = args['processes']
    chromosomes = args['chromosomes']
    outfile = args['outfile']
    all_names = args['sample_names']
    msr = args['min_snp_reads']
    mtr = args['min_total_reads']
    use_chr = args['use_chr']
    centromeres = args['centromeres']
    compressed = args['compressed']

    n_workers = min(len(chromosomes), threads)
    threads_per_worker = int(threads / n_workers)

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
    params = [(stem, all_names, ch, outfile + f'.{ch}', threads_per_worker, 
               chr2centro[ch][0], chr2centro[ch][1], msr, mtr, compressed)
              for ch in chromosomes]
    
    # dispatch workers
    """
    for p in params:
        run_chromosome_wrapper(p)
    """
    with Pool(n_workers) as p:
        p.map(run_chromosome_wrapper, params)
    
    # merge all BB files together to get the one remaining BB file
    outfiles = [a[3] for a in params]
    bbs = [pd.read_table(bb) for bb in outfiles]
    big_bb = pd.concat(bbs)
    big_bb.to_csv(outfile, index = False, sep = '\t')

    # Uncommend to remove all intermediate bb files (once I know the previous steps work)
    # [os.remove(f) for f in outfiles]

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

    return arr, thresholds

def adaptive_bins_arm(snp_thresholds, total_counts, snp_positions, snp_counts,
                            min_snp_reads = 2000, min_total_reads = 5000):
    """
    Compute adaptive bins for a single chromosome arm.
    Parameters:
        snp_thresholds: length <n> array of 1-based genomic positions of candidate bin thresholds
            
        total_counts: <n> x <2d> np.ndarray 
            entry [i, 2j] contains the number of reads starting in (snp_thresholds[i], snp_thresholds[i + 1]) in sample j (only the first n-1 positions are populated)
            entry [i, 2j + 1] contains the number of reads covering position snp_thresholds[i] in sample j


        snp_positions: length <m> list of 1-based genomic positions of SNPs
        snp_counts: <m> x <d> np.ndarray containing the number of overlapping reads at each of the <n - 1> snp positions in <d> samples
                
        min_snp_reads: the minimum number of SNP-covering reads required in each bin and each sample
        min_total_reads: the minimum number of total reads required in each bin and each sample
        
    """
    assert len(snp_thresholds) == len(total_counts)
    assert len(snp_positions) == len(snp_counts)
    assert len(snp_positions) == len(snp_thresholds) - 1, (len(snp_positions), len(snp_thresholds))
    assert np.all(snp_positions > snp_thresholds[0])
    assert len(snp_positions) > 0
    assert len(snp_thresholds) >= 2
    
    n_samples = int(total_counts.shape[1] / 2)
            
    
    # number of reads that start between snp_thresholds[i] and snp_thresholds[i + 1]
    even_index = np.array([i * 2 for i in range(int(n_samples))], dtype = np.int8)
    # number of reads overlapping position snp_thresholds[i]
    odd_index = np.array([i * 2 + 1 for i in range(int(n_samples))], dtype = np.int8) 
    
    bin_total = np.zeros(n_samples) 
    bin_snp = np.zeros(n_samples - 1)
        
    starts = []
    ends = []

    my_start = snp_thresholds[0]
    prev_threshold = snp_thresholds[0]
    
    rdrs = []
    totals = []
    i = 1
    while i < len(snp_thresholds - 1):
        # Extend the current bin to the next threshold
        next_threshold = snp_thresholds[i]        
        
        # add the intervening reads to the current bin
        # adding SNP reads
        assert snp_positions[i - 1] >= snp_thresholds[i - 1]
        assert snp_positions[i - 1] <= snp_thresholds[i], (i, snp_positions[i - 1], snp_thresholds[i])

        bin_snp += snp_counts[i - 1]
        
        # adding total reads
        bin_total += total_counts[i - 1, even_index]
            
        if np.all(bin_snp >= min_snp_reads) and np.all(bin_total - total_counts[i, odd_index]  >= min_total_reads):
            
            # end this bin
            starts.append(my_start)
            ends.append(next_threshold)
            
            # to get the total reads, subtract the number of reads covering the threshold position
            bin_total -= total_counts[i, odd_index]
            totals.append(bin_total)
            
            # compute RDR
            rdrs.append(bin_total[1:] / bin_total[0])
            
            # and start a new one
            bin_total = np.zeros(n_samples)
            bin_snp = np.zeros(n_samples - 1)    
            my_start = ends[-1] + 1
            
        prev_threshold = next_threshold
        i += 1
    
    # add whatever excess at the end to the last bin
    if ends[-1] < snp_thresholds[-1]:
        # combine the last complete bin with the remainder
        last_start_idx = np.where(snp_thresholds == starts[-1] - 1)[0][0]        
        bin_total = np.sum(total_counts[last_start_idx:, even_index], 
                           axis = 0) - total_counts[-1, odd_index] 
        totals[-1] = bin_total
        rdrs[-1] = bin_total[1:] / bin_total[0]

    return starts, ends, totals, rdrs

def EM(totals_in, alts_in, start, tol=1e-6):
    """
    Adapted from chisel/Combiner.py
    """
    totals = np.array(totals_in)
    alts = np.array(alts_in)
    assert totals.size == alts.size and 0 < start < 1 and np.all(totals >= alts)
    if np.all(np.logical_or(totals == alts, alts == 0)):
        return 0.0, np.ones(len(totals_in)) * 0.5, 0.0
    else:
        baf = start
        prev = None
        while prev is None or abs(prev - baf) >= tol:
            prev = baf
            
            # E-step
            assert 0 + tol < baf < 1 - tol, (baf, totals, alts, start)
            M = (totals - 2*alts) * np.log(baf) + (2*alts - totals) * np.log(1.0 - baf)
            M = np.exp(np.clip(a = M, a_min = -100, a_max = 100))
            phases = np.reciprocal(1 + M)
            
            # M-step
            baf = float(np.sum(totals * (1 - phases) + alts * (2*phases - 1))) / float(np.sum(totals))

        assert 0 + tol < baf < 1 - tol, (baf, totals, alts, start)
        lpmf = binom.logpmf
        log_likelihood = float(np.sum(phases * lpmf(k = alts, n = totals, p = baf) + 
                             (1 - phases) * lpmf(k = alts, n = totals, p = 1-baf)))
        return baf, phases, log_likelihood

def apply_EM(totals_in, alts_in, n_samples, significance, cutoff):
    baf, phases, logl = max((EM(totals_in, alts_in, start = st) 
                    for st in np.linspace(0.01, 0.49, 50)), key=(lambda x : x[2]))
    refs = totals_in - alts_in
    phases = phases.round().astype(np.int8)
    return baf, np.sum(np.choose(phases, [refs, alts_in])), np.sum(np.choose(phases, [alts_in, refs]))
   
def compute_baf_task(bin_snps):
    samples = sorted(bin_snps.SAMPLE.unique())
    result = {}
    for sample in samples:
        # Compute BAF
        my_snps = bin_snps[bin_snps.SAMPLE == sample]

        baf, alpha, beta = apply_EM(totals_in = my_snps.TOTAL, 
                                                 alts_in = my_snps.ALT, 
                                                 n_samples = 100, 
                                                 significance = 0.05, 
                                                 cutoff = 0.3)   
        n_snps = len(my_snps)
        cov = np.sum(alpha + beta) / len(my_snps)
        
        # doing this check in bin_chromosome() instead
        #assert np.isclose(np.sum(alpha + beta), my_snps.TOTAL.sum()), (alpha + beta, my_snps.TOTAL.sum, my_snps.iloc[0].POS)
        
        result[sample] = n_snps, cov, baf, alpha, beta
    return result

def compute_bafs_parallel(snp_dfs, threads):
    #print(datetime.now(), "Starting BAF computation")
    
    #with Pool(threads) as p:
    #    result = p.map(compute_baf_task, snp_dfs)
    result = [compute_baf_task(d) for d in snp_dfs]
        
    #print(datetime.now(), "Done BAF computation")

    return result

def merge_data(bins, dfs, bafs, sample_names, chromosome):
    """
    Merge bins data (starts, ends, total counts, RDRs) with SNP data and BAF data for each bin.
    Parameters:
    bins: output from call to adaptive_bins_arm
    dfs: (only for troubleshooting) pandas DataFrame, each containing the SNP information for the corresponding bin
    bafs: the ith element is the output from compute_baf_task(dfs[i])
    
    Produces a BB file with a few additional columns.
    """

    rows = []
    sample_names = sample_names[1:] # ignore the normal sample (first in the list)
    for i in range(len(bins[0])):
        start = bins[0][i]
        end = bins[1][i]
        totals = bins[2][i]
        rdrs = bins[3][i]
        
        assert dfs[i].POS.min() >= start and dfs[i].POS.max() <= end, (start, end, i, 
                                                                       dfs[i].POS.min(), dfs[i].POS.max())
        
        snpcounts_from_df = dfs[i].pivot(index = 'POS', columns = 'SAMPLE', values = 'TOTAL').sum(axis = 0)

        for j in range(len(sample_names)):
            sample = sample_names[j]
            total = totals[j + 1]
            normal_reads = totals[0]
            rdr = rdrs[j]
            nsnps, cov, baf, alpha, beta = bafs[i][sample]
            
            assert snpcounts_from_df[sample] == alpha + beta, (i, sample)
            
            rows.append([chromosome, start, end, sample, rdr,  nsnps, cov, alpha, beta, baf, total, normal_reads])
                    
    return pd.DataFrame(rows, columns = ['CHR', 'START', 'END', 'SAMPLE', 'RD', '#SNPS', 'COV', 'ALPHA', 'BETA', 'BAF', 'TOTAL_READS', 'NORMAL_READS'])

def run_chromosome(stem, all_names, chromosome, outfile, nthreads, centromere_start, centromere_end, 
         min_snp_reads, min_total_reads, compressed):
    """
    Perform adaptive binning and infer BAFs to produce a HATCHet BB file for a single chromosome.
    """

    try:
        tracemalloc.start()

        if os.path.exists(outfile):
            sp.log(msg=f"Output file already exists, skipping chromosome {chromosome}\n", level = "INFO")
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

        #sp.log(msg=f"Reading SNPs file for chromosome {chromosome}\n", level = "INFO")
        # Load SNP positions and counts for this chromosome
        positions, snp_counts, snpsv = read_snps(os.path.join(stem, 'baf', 'bulk.1bed'), chromosome, all_names)
        
        thresholds = np.trunc(np.vstack([positions[:-1], positions[1:]]).mean(axis = 0)).astype(np.uint32)
        last_idx_p = np.argwhere(thresholds > centromere_start)[0][0]
        first_idx_q = np.argwhere(thresholds > centromere_end)[0][0]
        all_thresholds = np.concatenate([[1], thresholds[:last_idx_p], [centromere_start], 
                                        [centromere_end], thresholds[first_idx_q:]])
        
        
        #sp.log(msg=f"Loading counts for chromosome {chromosome}\n", level = "INFO")
        # Load read count arrays from file (also adds end of chromosome as a threshold)
        total_counts, complete_thresholds = form_counts_array(starts_files, perpos_files, all_thresholds, chromosome) 


        sp.log(msg=f"Binning p arm of chromosome {chromosome}\n", level = "INFO")
        if len(np.where(positions <= centromere_start)[0]) > 0:
            # There may not be a SNP between the centromere end and the next SNP threshold
            # Goal for p arm is to END at the FIRST threshold that is AFTER the LAST SNP BEFORE the centromere
            last_snp_before_centromere = positions[np.where(positions <= centromere_start)[0][-1]]
            last_threshold_before_centromere = complete_thresholds[np.where(complete_thresholds >= last_snp_before_centromere)[0][0]]  
            
            p_idx = np.where(complete_thresholds <= last_threshold_before_centromere)[0]
            p_thresholds = complete_thresholds[p_idx]
            p_counts = total_counts[p_idx]

            p_snp_idx = np.where(positions <= last_threshold_before_centromere)[0]
            p_positions = positions[p_snp_idx]
            p_snpcounts = snp_counts[p_snp_idx]

            # Identify bins 
            #print(datetime.now(), "### Binning p arm ###")
            bins_p = adaptive_bins_arm(snp_thresholds = p_thresholds, 
                                    total_counts = p_counts, 
                                    snp_positions = p_positions, 
                                    snp_counts = p_snpcounts,
                                    min_snp_reads = min_snp_reads, min_total_reads = min_total_reads)
        
            starts_p = bins_p[0]
            ends_p = bins_p[1]
            # Partition SNPs for BAF inference
            dfs_p = [snpsv[(snpsv.POS >= starts_p[i]) & (snpsv.POS <= ends_p[i])] for i in range(len(starts_p))]
            
            for i in range(len(dfs_p)):
                assert np.all(dfs_p[i].pivot(index = 'POS', columns = 'SAMPLE', values = 'TOTAL').sum(axis = 0) >= min_snp_reads), i
                    
            # Infer BAF
            bafs_p = compute_bafs_parallel(dfs_p, threads = nthreads) 
            bb_p = merge_data(bins_p, dfs_p, bafs_p, all_names, chromosome)
        else:
            sp.log(msg=f"No SNPs found in p arm for {chromosome}\n", level = "INFO")
            bb_p = None
            
            
        sp.log(msg=f"Binning q arm of chromosome {chromosome}\n", level = "INFO")
        
        if len(np.where(positions >= centromere_end)[0]) > 0:
            # There may not be a SNP between the centromere end and the next SNP threshold
            # Goal for q arm is to start at the latest threshold that is before the first SNP after the centromere
            first_snp_after_centromere = positions[np.where(positions >= centromere_end)[0][0]]
            first_threshold_after_centromere = complete_thresholds[np.where(complete_thresholds <= first_snp_after_centromere)[0][-1]]  
            
            q_idx = np.where(complete_thresholds >= first_threshold_after_centromere)[0]
            q_thresholds = complete_thresholds[q_idx]
            q_counts = total_counts[q_idx]

            q_snp_idx = np.where(positions >= first_threshold_after_centromere)[0]
            q_positions = positions[q_snp_idx]
            q_snpcounts = snp_counts[q_snp_idx]
            
            #print(datetime.now(), "### Binning q arm ###")
            bins_q = adaptive_bins_arm(snp_thresholds = q_thresholds, 
                                    total_counts = q_counts, 
                                    snp_positions = q_positions, 
                                    snp_counts = q_snpcounts,
                                    min_snp_reads = min_snp_reads, min_total_reads = min_total_reads)

            starts_q = bins_q[0]
            ends_q = bins_q[1]
                
            # Partition SNPs for BAF inference
            dfs_q = [snpsv[(snpsv.POS >= starts_q[i]) & (snpsv.POS <= ends_q[i])] for i in range(len(starts_q))]
            
            for i in range(len(dfs_q)):
                assert np.all(dfs_q[i].pivot(index = 'POS', columns = 'SAMPLE', values = 'TOTAL').sum(axis = 0) >= min_snp_reads), i
                    
            # Infer BAF
            bafs_q = compute_bafs_parallel(dfs_q, threads = nthreads)
            bb_q = merge_data(bins_q, dfs_q, bafs_q, all_names, chromosome)
        else:
            sp.log(msg=f"No SNPs found in q arm for {chromosome}\n", level = "INFO")
            bb_q = None
        
        if bb_p is None and bb_q is None:
            raise ValueError(sp.error(f"No SNPs found on either arm of chromosome {chromosome}!"))
        
        bb = pd.concat([bb_p, bb_q])
        bb.to_csv(outfile, index = False, sep = '\t')
        np.savetxt(outfile + '.totalcounts', total_counts)
        np.savetxt(outfile + '.thresholds', complete_thresholds)
        
        #sp.log(msg=f"Loading counts for chromosome {chromosome}\n", level = "INFO")
        # Load read count arrays from file (also adds end of chromosome as a threshold)
        total_counts, complete_thresholds

        sp.log(msg=f"Done chromosome {chromosome}\n", level ="INFO")
        current, peak = tracemalloc.get_traced_memory()
        sp.log(msg=f"Chr {chromosome} -- Current memory usage is {int(current / 10**6)}MB; Peak was {int(peak / 10**6)}MB",
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