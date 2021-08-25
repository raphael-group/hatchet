from multiprocessing import Pool
import os
from datetime import datetime
import numpy as np
import pandas as pd
from scipy.stats import binom, norm
import gzip
import subprocess
import traceback
from importlib.resources import path
import hatchet.resources

from .ArgParsing import parse_combine_counts_args
from . import Supporting as sp


def main(args=None):
    sp.log(msg="# Parsing and checking input arguments\n", level="STEP")
    args = parse_combine_counts_args(args)
    sp.logArgs(args, 80)
    
    baffile = args['baffile']
    threads = args['processes']
    chromosomes = args['chromosomes']
    outfile = args['outfile']
    all_names = args['sample_names']
    msr = args['min_snp_reads']
    mtr = args['min_total_reads']
    use_chr = args['use_chr']
    phase = args['phase']
    blocksize = args['blocksize']
    max_snps_per_block = args['max_snps_per_block']
    test_alpha = args['test_alpha']
    use_em = args['use_em']

    n_workers = min(len(chromosomes), threads)

    # Read in centromere locations table
    with path(hatchet.resources, f'{args["ref_version"]}.centromeres.txt') as centromeres:
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


    isX = {ch:ch.endswith('X') for ch in chromosomes}
    isY = {ch:ch.endswith('Y') for ch in chromosomes}

    # form parameters for each worker
    params = [(baffile, all_names, ch, outfile + f'.{ch}', 
               chr2centro[ch][0], chr2centro[ch][1], msr, mtr, args['array'],
               isX[ch] or isY[ch], use_em, phase, blocksize, max_snps_per_block, test_alpha)
              for ch in chromosomes]
    # dispatch workers
    """
    for p in params:
        run_chromosome_wrapper(p)
    """
    with Pool(n_workers) as p:
        p.map(run_chromosome_wrapper, params)
    
    sp.log(msg="# Merging per-chromosome bb files and correcting read counts\n", level="STEP")
    # merge all BB files together to get the one remaining BB file
    outfiles = [a[3] for a in params]
    bbs = [pd.read_table(bb, dtype = {'#CHR':str}) for bb in outfiles]
    big_bb = pd.concat(bbs)
    big_bb = big_bb.sort_values(by = ['#CHR', 'START', 'SAMPLE'])
    
    big_bb['CORRECTED_READS'] = np.NAN
    
    
    # For each sample, correct read counts to account for differences in coverage (as in HATCHet)
    # (i.e., multiply read counts by total-reads-normal/total-reads-sample)
    rc = pd.read_table(args['totalcounts'], header = None, names = ['SAMPLE', '#READS'])
    nreads_normal = rc[rc.SAMPLE == 'normal'].iloc[0]['#READS']
    for sample in rc.SAMPLE.unique():
        if sample == 'normal':
            continue
        nreads_sample = rc[rc.SAMPLE == sample].iloc[0]['#READS']
        correction = nreads_normal / nreads_sample
        my_bb = big_bb[big_bb.SAMPLE == sample]
        
        # Correct the tumor reads propotionally to the total reads in corresponding samples
        big_bb.loc[big_bb.SAMPLE == sample, 'CORRECTED_READS'] = (my_bb.TOTAL_READS * correction).astype(np.int64)
        
        # Recompute RDR according to the corrected tumor reads
        big_bb.loc[big_bb.SAMPLE == sample, 'RD'] = big_bb.loc[big_bb.SAMPLE == sample, 'CORRECTED_READS'] / big_bb.loc[big_bb.SAMPLE == sample, 'NORMAL_READS']
    
    if not "NORMAL_READS" in big_bb:
        sp.log("# NOTE: adding NORMAL_READS column to bb file", level = "INFO")
        big_bb['NORMAL_READS'] = (big_bb.CORRECTED_READS / big_bb.RD).astype(np.uint32)


    autosomes = set([ch for ch in big_bb['#CHR'] if not (ch.endswith('X') or ch.endswith('Y'))])
    big_bb[big_bb['#CHR'].isin(autosomes)].to_csv(outfile, index = False, sep = '\t')
    
    # Convert intervals from closed to half-open to match .1bed/HATCHet standard format
    big_bb.END = big_bb.END + 1 
    big_bb.to_csv(outfile + '.withXY', index = False, sep = '\t')

    # Remove intermediate BB files
    [os.remove(f) for f in outfiles]
    
    sp.log(msg="# Done\n", level="STEP")

def read_snps(baf_file, ch, all_names, phasefile = None):
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
    
    if phasefile is not None:
        # Read in phasing output
        phases = pd.read_table(phasefile, compression  = 'gzip', comment = '#', names = 'CHR\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPHASE'.split(), usecols = ['CHR', 'POS', 'PHASE'], quoting = 3, low_memory = False, dtype = {'CHR':object, 'POS':np.uint32})
        phases['FLIP'] = phases.PHASE.str.contains('1\|0').astype(np.int8)
        phases['NOFLIP'] = phases.PHASE.str.contains('0\|1').astype(np.int8)

        # Drop entries without phasing output
        phases = phases[phases.FLIP + phases.NOFLIP > 0]

        # For exact duplicate entries, drop one
        phases = phases.drop_duplicates()
        
        # For duplicate entries with the same (CHR, POS) but different phase, drop all
        phases = phases.drop_duplicates(subset = ['CHR', 'POS'], keep = False)

        # Merge tables: keep only those SNPs for which we have phasing output
        snpsv = pd.merge(snpsv, phases, on = ['CHR', 'POS'], how = 'left')
    
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
        last_start_idx = np.where((snp_thresholds == starts[-1] - 1) | (snp_thresholds == starts[-1]))[0][0]        
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
            M = (totals - 2.0*alts) * np.log(baf) + (2.0*alts - totals) * np.log(1.0 - baf)
            M = np.exp(np.clip(a = M, a_min = -100, a_max = 100))
            phases = np.reciprocal(1 + M)
            
            # M-step
            baf = float(np.sum(totals * (1 - phases) + alts * (2.0*phases - 1))) / float(np.sum(totals))

        assert 0 + tol < baf < 1 - tol, (baf, totals, alts, start)
        lpmf = binom.logpmf
        log_likelihood = float(np.sum(phases * lpmf(k = alts, n = totals, p = baf) + 
                             (1 - phases) * lpmf(k = alts, n = totals, p = 1-baf)))
        return baf, phases, log_likelihood

def apply_EM(totals_in, alts_in):
    baf, phases, logl = max((EM(totals_in, alts_in, start = st) 
                    for st in np.linspace(0.01, 0.49, 50)), key=(lambda x : x[2]))
    refs = totals_in - alts_in
    phases = phases.round().astype(np.int8)
    return baf, np.sum(np.choose(phases, [refs, alts_in])), np.sum(np.choose(phases, [alts_in, refs]))
   
def max_likelihood_hardphasing(alts, totals, n_candidates = 1000):
    """
    Essentially a "maximum-maximum" approach: computes the maximum-likelihood h for many values of p,
    then uses these h's ("phasing" for each SNP/block) for each p to evaluate its likelihood 
    (then takes p which maximizes likelihood).
    """
    candidates = np.linspace(0, 0.5, n_candidates)
    likelihoods_noflip = binom.logpmf(k = np.repeat(alts, len(candidates)), 
                                  n = np.repeat(totals, len(candidates)), 
                                  p = np.tile(candidates, len(alts))).reshape(len(alts), len(candidates))
    likelihoods_flip = binom.logpmf(k = np.repeat(alts, len(candidates)), 
                                    n = np.repeat(totals, len(candidates)), 
                                    p = 1 - np.tile(candidates, len(alts))).reshape(len(alts), len(candidates))
    
    best_phasings = np.argmax(np.stack([likelihoods_noflip, likelihoods_flip]), axis = 0)
    likelihoods_best = np.max(np.stack([likelihoods_noflip, likelihoods_flip]), axis = 0)
    best_idx = np.argmax(likelihoods_best.sum(axis = 0))
    best_p = candidates[best_idx]
    
    return best_p, 1 - best_phasings[:, best_idx]
   
def apply_MM(totals_in, alts_in):
    baf, phaess = max_likelihood_hardphasing(alts = alts_in, totals = totals_in)
    
    refs = totals_in - alts_in 
    alpha = np.sum(np.choose(phaess, [refs, alts_in]))
    beta = np.sum(np.choose(phaess, [alts_in, refs]))   
    if alpha < beta:
        return baf, alpha, beta
    else:
        return baf, beta, alpha

def compute_baf_task(bin_snps, blocksize, max_snps_per_block, test_alpha, use_em):
    """
    Estimates the BAF for the bin containing exactly <bin_snps> SNPs.
    <bin_snps> is a dataframe with at least ALT and REF columns containing read counts.
    <blocksize>, <max_snps_per_block>, and <test_alpha> are used only for constructing phase blocks.
    """
    
    samples = sorted(bin_snps.SAMPLE.unique())
    result = {}
    
    phasing = "PHASE" in bin_snps.columns
    if phasing:
        # TODO: select the highest coverage sample to use for constructing phase blocks?
        # or maybe use all samples for BAF test?
        all_phase_data = [phase_blocks_sequential(d, blocksize = blocksize, 
                                             max_snps_per_block = max_snps_per_block,
                                             alpha = test_alpha) for _, d in bin_snps.groupby('SAMPLE')]
        phase_data = merge_phasing(bin_snps, all_phase_data)
        
    for sample in samples:
        # Compute BAF
        my_snps = bin_snps[bin_snps.SAMPLE == sample]
        n_snps = len(my_snps)

        if phasing:
            my_snps = collapse_blocks(my_snps, *phase_data)
            
        if use_em:
            baf, alpha, beta = apply_EM(totals_in = my_snps.TOTAL, 
                                        alts_in = my_snps.ALT)   
        else:
            baf, alpha, beta = apply_MM(totals_in = my_snps.TOTAL, 
                                        alts_in = my_snps.ALT)   
            
        cov = np.sum(alpha + beta) / n_snps

        result[sample] = n_snps, cov, baf, alpha, beta
    return result

def merge_phasing(df, all_phase_data):
    N = len(df)
    
    if len(all_phase_data) == 1:
        return all_phase_data[0]
    
    singletons = set()
    orphans = set()
    breaks = set()
    for i in range(len(all_phase_data)):
        # Merge singletons and orphans
        singletons = singletons.union(all_phase_data[i][1])
        orphans = orphans.union(all_phase_data[i][2])
    
        # Identify breaks
        blocking = all_phase_data[i][0]
        j = 0
        for j in range(len(blocking) - 1):
            if blocking[j][-1] == blocking[j + 1][0] - 1:
                breaks.add(blocking[j][-1])
    
    # Split blocking in sample 1 to match other breaks
    blocks = []
    for b in all_phase_data[0][0]:
        br = np.argwhere([a in breaks or a in orphans for a in b[:-1]])
        blocks.extend([list(a) for a in np.split(b, br.flatten() + 1)])
    #print(blocks)
        
    #print([a[0] for a in blocks if len(a) == 1])
    orphans = orphans.union([a[0] for a in blocks if len(a) == 1])
    blocks = [b for b in blocks if len(b) > 1]
        
    return blocks, singletons, orphans, all_phase_data[0][-2]

def binom_prop_test(alt1, ref1, flip1, alt2, ref2, flip2, alpha = 0.1):
    """
    Returns True if there is sufficient evidence that SNPs 1 and 2 should not be merged, False otherwise.
    """
    if np.isnan(flip1) or np.isnan(flip2):
        return True
    
    x1 = ref1 if flip1 > 0 else alt1
    x2 = ref2 if flip2 > 0 else alt2
    n1 = max(1, alt1 + ref1)
    n2 = max(1, alt2 + ref2)
    
    p1 = x1 / n1
    p2 = x2 / n2
    
    p = (x1 + x2) / (n1 + n2)
    if p == 0:
        p = 0.001
    elif p == 1:
        p = 0.999
    denom = np.sqrt(p * (1 - p) * (1/n1 + 1/n2))
    z = (p1 - p2) / denom

    if z > 0:
        return 2 * (1 - norm.cdf(z)) < alpha
    else:
        return 2 * norm.cdf(z) < alpha
    
def block_segment(df, blocksize, max_snps_per_block):
    """
    Given a pandas dataframe containing read counts for a contiguous segment of SNPs,
    collects SNPs into phase blocks of size at most <blocksize> containing at most
    <max_snps_per_block> SNPs each.
    """
    splits = []
    df = df.copy().reset_index(drop = True)
    
    # walk along segment until either blocksize or max_snps_per_block is reached
    my_start = df.iloc[0].POS
    n = 1
    for i, r in df.iterrows():
        if i == 0:
            continue
        
        if n >= max_snps_per_block or r.POS - my_start > blocksize:
            splits.append(i)
            my_start = r.POS
            n = 1
        else:
            n += 1
    return np.split(np.arange(len(df)), splits)
        
def consecutive(data, stepsize=1):
    return np.split(data, np.where(np.diff(data) != stepsize)[0]+1)

def phase_blocks_sequential(df, blocksize = 50e3, max_snps_per_block = 10, alpha = 0):
    if len(df) == 0:
        return []
    
    """NOTE: this currently only supports single-sample data"""
    assert len(df.SAMPLE.unique()) == 1    
    sample = df.SAMPLE.unique()[0]
    
    assert len(df.CHR.unique()) == 1
    ch = df.CHR.unique()[0]
    
    df = df.copy()
    df = df.reset_index(drop = True)
    
    # Identify contiguous segments of SNPs with phasing information
    segments = consecutive(np.where(~df.PHASE.isna())[0])
    
    if alpha < 1:
        ### Use binomial test to split adjacent SNPs with very different phased BAFs
        # (further dividing the segments identified above)
        # identify indices of adj SNPs that are significantly different
        df_merged = pd.concat([df, df.shift(-1).add_prefix('next_')], axis=1)
        # no_merge is the left index of each such pair
        no_merge = np.where(df_merged[:-1].apply(lambda x: binom_prop_test(x.ALT, x.REF, x.FLIP, x.next_ALT, x.next_REF, x.next_FLIP, alpha = alpha), axis=1))[0]
        # stack the adjacent indices (top is left, bottom is right = left + 1)
        nm = np.vstack([no_merge, no_merge + 1])

        new_segments = []
        orphans = set()
        for seg in segments:
            # identify adjacent pairs that violate test (p<alpha) and are both in this segment
            violations = np.all(np.isin(nm, seg), axis = 0)

            # map violations to indices in the segment
            keys = nm[1, np.nonzero(violations)][0]
            to_split = [np.where(seg == k)[0][0] for k in keys]
            splitted = np.split(seg, to_split)
            orphans = orphans.union([a[0] for a in splitted if len(a) == 1])
            new_segments.extend([a for a in splitted if len(a) > 1])

        segments = new_segments
    
    blocks = []
    
    for segment in segments:
        if len(segment) == 0:
            pass
        elif len(segment) == 1:
            blocks.append(list(segment))
        else:
            my_df = df.iloc[segment]
            ## Run blocking on my_df
            my_blocks = block_segment(my_df, blocksize, max_snps_per_block)
            mapping = {i:my_df.iloc[i].name for i in range(len(my_df))}
            blocks.extend([[mapping[a] for a in b] for b in my_blocks])
            
    blocks = sorted(blocks, key = lambda x: x[0])
    singletons = set(np.where(df.PHASE.isna())[0])
    
    return blocks, singletons, orphans, ch, sample

def collapse_blocks(df, blocks, singletons, orphans, ch):
    ### Construct blocked 1-bed table
    # First, merge blocks for those SNPs that are in the same block
    assert len(df.SAMPLE.unique()) == 1
    sample = list(df.SAMPLE.unique())[0]
    
    rows = []
    i = 0
    j = 0
    while i < len(df):
        if i in singletons or i in orphans:
            r = df.iloc[i]
            rows.append([ch, r.POS, r.POS, sample, r.ALT, r.REF, r.TOTAL, 1])
            i += 1
        else:
            block = blocks[j]
            assert i in block, (block, i)
            i += len(block)
            j += 1
    
            my_snps = df.iloc[block]
            start = my_snps.POS.min()
            end = my_snps.POS.max()
            alt = np.sum(my_snps.FLIP * my_snps.REF + (1 - my_snps.FLIP) * my_snps.ALT).astype(np.uint64)
            ref = np.sum(my_snps.FLIP * my_snps.ALT + (1 - my_snps.FLIP) * my_snps.REF).astype(np.uint64)
            total = np.sum(my_snps.TOTAL)
            n_snps = len(my_snps)
            rows.append([ch, start, end, sample, alt, ref, total, n_snps])
        
    return pd.DataFrame(rows, columns = ['CHR', 'START', 'END', 'SAMPLE', 'ALT', 'REF', 'TOTAL', '#SNPS'])

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
        
        if dfs is not None:
            assert dfs[i].POS.min() >= start and dfs[i].POS.max() <= end, (start, end, i, 
                                                                        dfs[i].POS.min(), dfs[i].POS.max())
            snpcounts_from_df = dfs[i].pivot(index = 'POS', columns = 'SAMPLE', values = 'TOTAL').sum(axis = 0)

        for j in range(len(sample_names)):
            sample = sample_names[j]
            total = totals[j + 1]
            normal_reads = totals[0]
            rdr = rdrs[j]

            if dfs is not None:
                nsnps, cov, baf, alpha, beta = bafs[i][sample]
                assert snpcounts_from_df[sample] == alpha + beta, (i, sample)
            else:
                nsnps, cov, baf, alpha, beta = 0, 0, 0, 0, 0
        
            rows.append([chromosome, start, end, sample, rdr,  nsnps, cov, alpha, beta, baf, total, normal_reads])
                    
    return pd.DataFrame(rows, columns = ['#CHR', 'START', 'END', 'SAMPLE', 'RD', '#SNPS', 'COV', 'ALPHA', 'BETA', 'BAF', 'TOTAL_READS', 'NORMAL_READS'])

def get_chr_end(stem, chromosome):
    fname = os.path.join(stem, chromosome + '.thresholds.gz')
    zcat = subprocess.Popen(('zcat', fname), stdout= subprocess.PIPE)
    tail = subprocess.Popen(('tail', '-1'), stdin=zcat.stdout, stdout = subprocess.PIPE)
    last_start = int(tail.stdout.read().decode('utf-8').strip())      
    
    return last_start 

def run_chromosome(baffile, all_names, chromosome, outfile, centromere_start, centromere_end, 
         min_snp_reads, min_total_reads, arraystem, xy, use_em, phasefile, blocksize, 
         max_snps_per_block, test_alpha):
    """
    Perform adaptive binning and infer BAFs to produce a HATCHet BB file for a single chromosome.
    """

    try:                
        if os.path.exists(outfile):
            sp.log(msg=f"Output file already exists, skipping chromosome {chromosome}\n", level = "INFO")
            return
        
        sp.log(msg=f"Loading intermediate files for chromosome {chromosome}\n", level = "INFO")
        total_counts = np.loadtxt(os.path.join(arraystem, f'{chromosome}.total'), dtype = np.uint32)
        complete_thresholds = np.loadtxt(os.path.join(arraystem, f'{chromosome}.thresholds'), dtype = np.uint32)
        
        # TODO: identify whether XX or XY, and only avoid SNPs/BAFs for XY
        if xy:
            sp.log(msg='Running on sex chromosome -- ignoring SNPs \n', level = "INFO")
            min_snp_reads = 0
            
            # TODO: do this procedure only for XY individuals
            ### construct dummy SNP positions and all-0 snpcounts array for binning
            before_centromere = complete_thresholds[complete_thresholds <= centromere_start]
            after_centromere = complete_thresholds[complete_thresholds >= centromere_end]
            positions_p = np.mean(np.vstack([before_centromere[:-1], before_centromere[1:]]), axis = 0).astype(np.uint64)
            positions_q = np.mean(np.vstack([after_centromere[:-1], after_centromere[1:]]), axis = 0).astype(np.uint64)
            positions = np.concatenate([positions_p, positions_q])
            snp_counts = np.zeros((len(positions), len(all_names) - 1), dtype = np.int8)
            snpsv = None

        else:
            #sp.log(msg=f"Reading SNPs file for chromosome {chromosome}\n", level = "INFO")
            # Load SNP positions and counts for this chromosome
            positions, snp_counts, snpsv = read_snps(baffile, chromosome, all_names, phasefile = phasefile)
            
            


        sp.log(msg=f"Binning p arm of chromosome {chromosome}\n", level = "INFO")
        if len(np.where(positions <= centromere_start)[0]) > 0:
            # There may not be a SNP between the centromere end and the next SNP threshold
            # Goal for p arm is to END at the FIRST threshold that is AFTER the LAST SNP BEFORE the centromere
            last_snp_before_centromere = positions[np.where(positions < centromere_start)[0][-1]]
            last_threshold_before_centromere = complete_thresholds[np.where(complete_thresholds > last_snp_before_centromere)[0][0]]  
            
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
               
            # Infer BAF
            if xy:
                # TODO: compute BAFs for XX
                dfs_p = None
                bafs_p = None
            else:
                dfs_p = [snpsv[(snpsv.POS >= starts_p[i]) & (snpsv.POS <= ends_p[i])] for i in range(len(starts_p))]         
                for i in range(len(dfs_p)):
                    
                    assert np.all(dfs_p[i].pivot(index = 'POS', columns = 'SAMPLE', values = 'TOTAL').sum(axis = 0) >= min_snp_reads), i
                bafs_p = [compute_baf_task(d, blocksize, max_snps_per_block, test_alpha, use_em) for d in dfs_p] 
                
                """
                except ZeroDivisionError:    
                    import pickle
                    with open(f"dfs_p_chr{chromosome}.pickle", 'wb') as f:
                        pickle.dump((dfs_p, blocksize, max_snps_per_block, test_alpha, use_em), f)
                    """
                
            bb_p = merge_data(bins_p, dfs_p, bafs_p, all_names, chromosome)
        else:
            sp.log(msg=f"No SNPs found in p arm for {chromosome}\n", level = "INFO")
            bb_p = None
            
            
        sp.log(msg=f"Binning q arm of chromosome {chromosome}\n", level = "INFO")
        
        if len(np.where(positions >= centromere_end)[0]) > 0:
            # There may not be a SNP between the centromere end and the next SNP threshold
            # Goal for q arm is to start at the latest threshold that is before the first SNP after the centromere
            first_snp_after_centromere = positions[np.where(positions > centromere_end)[0][0]]
            first_threshold_after_centromere = complete_thresholds[np.where(complete_thresholds < first_snp_after_centromere)[0][-1]]  
            
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
                
            if xy:
                dfs_q = None
                bafs_q = None
            else:
                # Partition SNPs for BAF inference
                dfs_q = [snpsv[(snpsv.POS >= starts_q[i]) & (snpsv.POS <= ends_q[i])] for i in range(len(starts_q))]
                
                for i in range(len(dfs_q)):
                    assert np.all(dfs_q[i].pivot(index = 'POS', columns = 'SAMPLE', values = 'TOTAL').sum(axis = 0) >= min_snp_reads), i
                        
                # Infer BAF
                bafs_q = [compute_baf_task(d, blocksize, max_snps_per_block, test_alpha, use_em) for d in dfs_q]
            
            bb_q = merge_data(bins_q, dfs_q, bafs_q, all_names, chromosome)
        else:
            sp.log(msg=f"No SNPs found in q arm for {chromosome}\n", level = "INFO")
            bb_q = None
        
        if bb_p is None and bb_q is None:
            raise ValueError(sp.error(f"No SNPs found on either arm of chromosome {chromosome}!"))
        
        bb = pd.concat([bb_p, bb_q])
        bb.to_csv(outfile, index = False, sep = '\t')
        #np.savetxt(outfile + '.totalcounts', total_counts)
        #np.savetxt(outfile + '.thresholds', complete_thresholds)
        
        sp.log(msg=f"Done chromosome {chromosome}\n", level ="INFO")
    except Exception as e: 
        print(f"Error in chromosome {chromosome}:")
        print(e)
        traceback.print_exc()
        raise e
    
def run_chromosome_wrapper(param):
    run_chromosome(*param)

if __name__ == '__main__':
    main()
