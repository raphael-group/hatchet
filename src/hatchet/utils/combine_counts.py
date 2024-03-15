from multiprocessing import Pool
import os
import subprocess
import traceback
from importlib.resources import path
import numpy as np
import pandas as pd
from scipy.stats import binom, norm
from scipy.special import softmax

import hatchet.data
from hatchet.utils.ArgParsing import parse_combine_counts_args
import hatchet.utils.Supporting as sp
from hatchet.utils.pon_normalization import correct_baf, pon_normalize_rdr
from hatchet.utils.rd_gccorrect import rd_gccorrect


def main(args=None):
    sp.log(msg='# Parsing and checking input arguments\n', level='STEP')
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
    multisample = args['multisample']
    nonormalFlag = args['nonormalFlag']
    ponfile = args['ponfile']
    referencefasta = args['referencefasta']

    n_workers = min(len(chromosomes), threads)

    # Read in centromere locations table
    with path(hatchet.data, f'{args["ref_version"]}.centromeres.txt') as centromeres:
        centromeres = pd.read_table(
            centromeres,
            header=None,
            names=['CHR', 'START', 'END', 'NAME', 'gieStain'],
        )
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
            raise ValueError(
                sp.error(f'Chromosome {ch} not found in centromeres file. Inspect file provided as -C argument.')
            )

    isX = {ch: ch.endswith('X') for ch in chromosomes}
    isY = {ch: ch.endswith('Y') for ch in chromosomes}

    # form parameters for each worker
    params = [
        (
            baffile,
            all_names,
            ch,
            outfile + f'.{ch}',
            chr2centro[ch][0],
            chr2centro[ch][1],
            msr,
            mtr,
            args['array'],
            isX[ch] or isY[ch],
            multisample,
            phase,
            blocksize,
            max_snps_per_block,
            test_alpha,
            nonormalFlag,
            args['segfile'],
        )
        for ch in chromosomes
    ]
    # dispatch workers
    """
    for p in params:
        run_chromosome_wrapper(p)
    """
    with Pool(n_workers) as p:
        p.map(run_chromosome_wrapper, params)

    sp.log(
        msg='# Merging per-chromosome bb files and correcting read counts\n',
        level='STEP',
    )
    # merge all BB files together to get the one remaining BB file
    if args['segfile']:
        outfiles = [a[3] + '.segfile' for a in params]
    else:
        outfiles = [a[3] for a in params]
    bbs = [pd.read_table(bb, dtype={'CHR': str}) for bb in outfiles]
    big_bb = pd.concat(bbs)
    big_bb = big_bb.sort_values(by=['CHR', 'START', 'SAMPLE'])

    big_bb['CORRECTED_READS'] = np.NAN

    # For each sample, correct read counts to account for differences in coverage (as in HATCHet)
    # (i.e., multiply read counts by total-reads-normal/total-reads-sample)
    rc = pd.read_table(args['totalcounts'], header=None, names=['SAMPLE', '#READS'])
    normal_name = all_names[0]
    nreads_normal = rc[rc.SAMPLE == normal_name].iloc[0]['#READS']
    if nonormalFlag:
        for sample, df in big_bb.groupby('SAMPLE'):
            big_bb.loc[big_bb.SAMPLE == sample, 'RD'] = (
                df['TOTAL_READS'] / (df['END'] - df['START']) / np.median(df['TOTAL_READS'] / (df['END'] - df['START']))
            )
    else:
        for sample in rc.SAMPLE.unique():
            if sample == normal_name:
                continue
            nreads_sample = rc[rc.SAMPLE == sample].iloc[0]['#READS']
            correction = nreads_normal / nreads_sample
            my_bb = big_bb[big_bb.SAMPLE == sample]

            # Correct the tumor reads propotionally to the total reads in corresponding samples
            big_bb.loc[big_bb.SAMPLE == sample, 'CORRECTED_READS'] = (my_bb.TOTAL_READS * correction).astype(np.int64)

            # Recompute RDR according to the corrected tumor reads
            big_bb.loc[big_bb.SAMPLE == sample, 'RD'] = (
                big_bb.loc[big_bb.SAMPLE == sample, 'CORRECTED_READS']
                / big_bb.loc[big_bb.SAMPLE == sample, 'NORMAL_READS']
            )

        if 'NORMAL_READS' not in big_bb:
            sp.log('# NOTE: adding NORMAL_READS column to bb file', level='INFO')
            big_bb['NORMAL_READS'] = (big_bb.CORRECTED_READS / big_bb.RD).astype(np.uint32)

    # Correct BAF when there is no normal sample. This correction useful only in LOH regions in high-purity samples
    if nonormalFlag:
        big_bb = correct_baf(big_bb)
        big_bb['BAF'] = big_bb['BAF'].round(5)

    if ponfile is not None:
        sp.log(
            msg='# Performing panel of normal read depth correction for cell-free data\n',
            level='STEP',
        )
        big_bb = pon_normalize_rdr(big_bb, args['ponfile'])

    sp.log(
        msg='# Performing GC bias correction on read depth signal\n',
        level='STEP',
    )
    big_bb = rd_gccorrect(big_bb, referencefasta)

    # Convert intervals from closed to half-open to match .1bed/HATCHet standard format
    big_bb.END = big_bb.END + 1
    autosomes = set([ch for ch in big_bb['CHR'] if not (ch.endswith('X') or ch.endswith('Y'))])
    big_bb[big_bb['CHR'].isin(autosomes)].to_csv(outfile, index=False, sep='\t')

    big_bb.to_csv(outfile + '.withXY', index=False, sep='\t')

    # Remove intermediate BB files
    [os.remove(f) for f in outfiles]

    sp.log(msg='# Done\n', level='STEP')


def read_snps(baf_file, ch, all_names, phasefile=None):
    """
    Read and validate SNP data for this patient (TSV table output from HATCHet deBAF.py).
    """
    all_names = [
        name for name in all_names if name != 'normal'
    ]   # remove normal sample -- not looking for SNP counts from normal

    # Read in HATCHet BAF table
    all_snps = pd.read_table(
        baf_file,
        names=['CHR', 'POS', 'SAMPLE', 'REF', 'ALT', 'REFC', 'ALTC'],
        dtype={
            'CHR': object,
            'POS': np.uint32,
            'SAMPLE': object,
            'ALT': np.uint32,
            'REF': np.uint32,
            'REFC': object,
            'ALTC': object,
        },
    )

    # Keep only SNPs on this chromosome
    snps = all_snps[all_snps.CHR == ch].sort_values(by=['POS', 'SAMPLE'])
    snps = snps.reset_index(drop=True)

    if len(snps) == 0:
        raise ValueError(
            sp.error(f'Chromosome {ch} not found in SNPs file (chromosomes in file: {all_snps.CHR.unique()})')
        )

    n_samples = len(all_names)
    if n_samples != len(snps.SAMPLE.unique()):
        raise ValueError(
            sp.error(f'Expected {n_samples} samples, found {len(snps.SAMPLE.unique())} samples in SNPs file.')
        )

    if set(all_names) != set(snps.SAMPLE.unique()):
        raise ValueError(
            sp.error(
                f'Expected sample names did not match sample names in SNPs file.\n\
                Expected: {sorted(all_names)}\n  Found:{sorted(snps.SAMPLE.unique())}'
            )
        )

    # Add total counts column
    snpsv = snps.copy()
    snpsv['TOTAL'] = snpsv.ALT + snpsv.REF

    if phasefile is not None:
        # Read in phasing output
        phases = pd.read_table(
            phasefile,
            compression='gzip',
            comment='#',
            names='CHR\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPHASE'.split(),
            usecols=['CHR', 'POS', 'PHASE'],
            quoting=3,
            low_memory=False,
            dtype={'CHR': object, 'POS': np.uint32},
        )
        phases['FLIP'] = phases.PHASE.str.contains('1|0', regex=False).astype(np.int8)  # noqa: W605
        phases['NOFLIP'] = phases.PHASE.str.contains('0|1', regex=False).astype(np.int8)  # noqa: W605

        # Drop entries without phasing output
        phases = phases[phases.FLIP + phases.NOFLIP > 0]

        # For exact duplicate entries, drop one
        phases = phases.drop_duplicates()

        # For duplicate entries with the same (CHR, POS) but different phase, drop all
        phases = phases.drop_duplicates(subset=['CHR', 'POS'], keep=False)

        # Merge tables: keep only those SNPs for which we have phasing output
        snpsv = pd.merge(snpsv, phases, on=['CHR', 'POS'], how='left')

    # Create counts array and find SNPs that are not present in all samples
    snp_counts = snpsv.pivot(index='POS', columns='SAMPLE', values='TOTAL')
    missing_pos = snp_counts.isna().any(axis=1)

    # Remove SNPs that are absent in any sample
    snp_counts = snp_counts.dropna(axis=0)
    snpsv = snpsv[~snpsv.POS.isin(missing_pos[missing_pos].index)]

    # Pivot table for dataframe should match counts array and have no missing entries
    check_pivot = snpsv.pivot(index='POS', columns='SAMPLE', values='TOTAL')
    assert np.array_equal(check_pivot, snp_counts), 'SNP file reading failed'
    assert not np.any(check_pivot.isna()), 'SNP file reading failed'
    assert np.array_equal(all_names, list(snp_counts.columns))   # make sure that sample order is the same

    return np.array(snp_counts.index), np.array(snp_counts), snpsv


def adaptive_bins_arm(
    snp_thresholds,
    total_counts,
    snp_positions,
    snp_counts,
    min_snp_reads=2000,
    min_total_reads=5000,
    nonormalFlag=False,
):
    """
    Compute adaptive bins for a single chromosome arm.
    Parameters:
        snp_thresholds: length <n> array of 1-based genomic positions of candidate bin thresholds

        total_counts: <n> x <2d> np.ndarray
            entry [i, 2j] contains the number of reads starting in [snp_thresholds[i], snp_thresholds[i + 1]) in sample
            j (only the first n-1 positions are populated) entry [i, 2j + 1] contains the number of reads covering
            position snp_thresholds[i] in sample j

        snp_positions: length <m> list of 1-based genomic positions of SNPs
            NOTE: this function requires that m = n-1 for convenience of programming (could be relaxed in a different
            implementation)
        snp_counts: <m> x <d> np.ndarray containing the number of overlapping reads at each of the <n - 1> snp
            positions in <d> samples

        min_snp_reads: the minimum number of SNP-covering reads required in each bin and each sample
        min_total_reads: the minimum number of total reads required in each bin and each sample

    """
    assert len(snp_thresholds) == len(total_counts)
    assert len(snp_positions) == len(snp_counts)
    assert len(snp_positions) == len(snp_thresholds) - 1, (
        len(snp_positions),
        len(snp_thresholds),
    )
    assert np.all(snp_positions > snp_thresholds[0])
    assert len(snp_positions) > 0
    assert len(snp_thresholds) >= 2

    n_samples = int(total_counts.shape[1] / 2)

    # number of reads that start between snp_thresholds[i] and snp_thresholds[i + 1]
    even_index = np.array([i * 2 for i in range(int(n_samples))], dtype=np.int8)
    # number of reads overlapping position snp_thresholds[i]
    odd_index = np.array([i * 2 + 1 for i in range(int(n_samples))], dtype=np.int8)

    bin_total = np.zeros(n_samples)
    if nonormalFlag:
        bin_snp = np.zeros(n_samples)
    else:
        bin_snp = np.zeros(n_samples - 1)

    starts = []
    ends = []

    my_start = snp_thresholds[0]

    rdrs = []
    totals = []
    bss = []
    i = 1
    while i < len(snp_thresholds - 1):
        # Extend the current bin to the next threshold
        next_threshold = snp_thresholds[i]

        # add the intervening reads to the current bin
        # adding SNP reads
        assert snp_positions[i - 1] >= snp_thresholds[i - 1]
        assert snp_positions[i - 1] <= snp_thresholds[i], (
            i,
            snp_positions[i - 1],
            snp_thresholds[i],
        )

        bin_snp += snp_counts[i - 1]

        # adding total reads
        bin_total += total_counts[i - 1, even_index]

        if np.all(bin_snp >= min_snp_reads) and np.all(bin_total - total_counts[i, odd_index] >= min_total_reads):

            # end this bin
            starts.append(my_start)
            ends.append(next_threshold)

            bss.append(bin_snp)

            # to get the total reads, subtract the number of reads covering the threshold position
            bin_total -= total_counts[i, odd_index]
            totals.append(bin_total)

            # compute RDR
            if nonormalFlag:
                rdrs.append(bin_total[0:] / bin_total[0:])
            else:
                rdrs.append(bin_total[1:] / bin_total[0])

            # and start a new one
            bin_total = np.zeros(n_samples)
            if nonormalFlag:
                bin_snp = np.zeros(n_samples)
            else:
                bin_snp = np.zeros(n_samples - 1)
            my_start = ends[-1] + 1

        i += 1

    # handle the case of 1 bin
    if len(ends) == 0:
        sp.log(
            msg='WARNING: found only 1 bin in chromosome arm, may not meet MSR and MTR\t',
            level='WARN',
        )
        assert len(starts) == 0
        starts.append(snp_thresholds[0])
        ends.append(snp_thresholds[-1])

        bin_total = np.sum(total_counts[:, even_index], axis=0) - total_counts[-1, odd_index]
        totals.append(bin_total)
        if nonormalFlag:
            rdrs.append(bin_total[0:] / bin_total[0:])
        else:
            rdrs.append(bin_total[1:] / bin_total[0])

    # add whatever excess at the end to the last bin
    if ends[-1] < snp_thresholds[-1]:
        # combine the last complete bin with the remainder
        last_start_idx = np.where((snp_thresholds == starts[-1] - 1) | (snp_thresholds == starts[-1]))[0][0]
        bin_total = np.sum(total_counts[last_start_idx:, even_index], axis=0) - total_counts[-1, odd_index]
        totals[-1] = bin_total
        if nonormalFlag:
            rdrs[-1] = bin_total[0:] / bin_total[0:]
        else:
            rdrs[-1] = bin_total[1:] / bin_total[0]
        ends[-1] = snp_thresholds[-1]

    return starts, ends, totals, rdrs


def EM(totals_in, alts_in, start, tol=1e-6):
    """
    Adapted from chisel/Combiner.py
    """
    totals = np.array(totals_in).reshape(-1)
    alts = np.array(alts_in).reshape(-1)
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
            M = (totals - 2.0 * alts) * np.log(baf) + (2.0 * alts - totals) * np.log(1.0 - baf)
            M = M.astype(float)
            M = np.exp(np.clip(a=M, a_min=-100, a_max=100))
            phases = np.reciprocal(1 + M)

            # M-step
            baf = float(np.sum(totals * (1 - phases) + alts * (2.0 * phases - 1))) / float(np.sum(totals))

        assert 0 + tol < baf < 1 - tol, (baf, totals, alts, start)
        lpmf = binom.logpmf
        log_likelihood = float(
            np.sum(
                phases * lpmf(k=list(alts), n=list(totals), p=baf)
                + (1 - phases) * lpmf(k=list(alts), n=list(totals), p=1 - baf)
            )
        )
        return [baf], phases, log_likelihood


def multisample_EM(totals, alts, start, tol=10e-6):

    refs = totals - alts
    assert refs.shape == alts.shape, 'Alternate and reference count arrays must have the same shape'
    assert 0 < start <= 0.5, 'Initial estimate must be in (0, 0.5]'

    n_samples, n_snps = alts.shape
    totals_sum = np.sum(totals, axis=1)
    assert np.all(totals_sum > 0), 'Every bin must have >0 SNP-covering reads in each sample!'

    e_arr = np.zeros((2, n_snps))

    if np.all(np.logical_or(refs == 0, alts == 0)):
        return np.array([0.0] * n_samples), np.ones(n_snps) * 0.5, 0.0
    else:
        theta = np.repeat(start, n_samples)
        prev = None
        while prev is None or np.all(np.abs(prev - theta) >= tol):
            prev = theta

            # Ensure that theta is not exactly 0 or 1 to avoid log(0)
            theta = np.clip(theta, tol, 1 - tol)

            # E-step in log space
            e_arr[0] = np.log(theta) @ refs + np.log(1 - theta) @ alts
            e_arr[1] = np.log(1 - theta) @ refs + np.log(theta) @ alts
            phases = softmax(e_arr, axis=0)[0, :]
            assert not np.any(np.isnan(phases)), (phases, e_arr)

            # M-step
            t1 = refs @ phases + alts @ (1 - phases)
            theta = t1 / totals_sum
            assert not np.any(np.isnan(theta)), (theta, t1)

    # If mean(BAF) > 0.5, flip phases accordingly
    if np.mean(theta) > 0.5:
        theta = np.clip(theta, tol, 1 - tol)

        theta = 1 - theta
        t1 = np.sum(np.log(theta) @ refs + np.log(1 - theta) @ alts, axis=0)
        t2 = np.sum(np.log(1 - theta) @ refs + np.log(theta) @ alts, axis=0)
        phases = softmax(np.vstack([t1, t2]), axis=0)[0, :]

    t1 = np.sum(np.log(theta) * (phases * refs).T, axis=0)
    t2 = np.sum(np.log(1 - theta) * (phases * alts).T, axis=0)
    t3 = np.sum(np.log(theta) * ((1 - phases) * alts).T, axis=0)
    t4 = np.sum(np.log(1 - theta) * ((1 - phases) * refs).T, axis=0)
    log_likelihood = np.sum(t1 + t2 + t3 + t4)

    return theta, phases, log_likelihood


def apply_EM(totals, alts, reference_based_haplo_blocks, multisample_em):
    if multisample_em:
        runs = (multisample_EM(totals, alts, b) for b in np.arange(0.05, 0.5, 0.05))
    else:
        runs = (EM(totals, alts, start=st) for st in np.linspace(0.01, 0.49, 50))
    baf, phases, _ = max(runs, key=lambda x: x[-1])
    refs = totals - alts
    phases = phases.round().astype(np.int8).reshape(-1)
    inverse_reference_haplo = pd.Series([[1 - ph for ph in hap] for hap in reference_based_haplo_blocks])
    haplo = np.choose(phases, [reference_based_haplo_blocks, inverse_reference_haplo])
    haploflat = [int(item) for sublist in haplo for item in sublist]
    haplostring = ','.join(list(map(str, haploflat)))
    # if not np.all(phases == (alts < refs)):
    #     sp.log(totals)
    return (
        baf,
        [int(np.sum(np.choose(phases, [alts[i], refs[i]]))) for i in range(len(alts))],
        [int(np.sum(np.choose(phases, [refs[i], alts[i]]))) for i in range(len(alts))],
        haplo,
        haploflat,
        haplostring,
    )


def compute_baf_task(bin_snps, blocksize, max_snps_per_block, test_alpha, multisample_em):
    """
    Estimates the BAF for the bin containing exactly <bin_snps> SNPs.
    <bin_snps> is a dataframe with at least ALT and REF columns containing read counts.
    <blocksize>, <max_snps_per_block>, and <test_alpha> are used only for constructing phase blocks.
    """

    samples = sorted(bin_snps.SAMPLE.unique())
    bafs, alphas, betas = [0] * len(samples), [0] * len(samples), [0] * len(samples)
    result = {}

    phasing = 'PHASE' in bin_snps.columns
    if phasing:
        # TODO: select the highest coverage sample to use for constructing phase blocks?
        # or maybe use all samples for BAF test?
        all_phase_data = [
            phase_blocks_sequential(
                d,
                blocksize=blocksize,
                max_snps_per_block=max_snps_per_block,
                alpha=test_alpha,
            )
            for _, d in bin_snps.groupby('SAMPLE')
        ]
        phase_data = merge_phasing(bin_snps, all_phase_data)
        grouped_snps = collapse_blocks(bin_snps, *phase_data, bin_snps.iloc[0].CHR)
    else:
        # each snp is it's on haplo block. [1] represents the minor count.
        # After EM algorithm, the haplostring tracks the haplotype with the minor
        # allele count (beta)
        grouped_snps = bin_snps.copy()
        grouped_snps['START'] = grouped_snps['POS']
        grouped_snps['HAPLO'] = len(grouped_snps) * [[1]]
    reference_based_haplo_blocks = grouped_snps[grouped_snps.SAMPLE == samples[0]].HAPLO

    totals = grouped_snps.pivot(index='SAMPLE', columns='START', values='TOTAL').to_numpy().astype(np.uint64)
    alts = grouped_snps.pivot(index='SAMPLE', columns='START', values='ALT').to_numpy().astype(np.uint64)

    # compute phasing once using all samples
    if multisample_em:
        bafs, alphas, betas, _, haploflat, haplostring = apply_EM(
            totals, alts, reference_based_haplo_blocks, multisample_em
        )

    for i in range(len(samples)):
        sample = samples[i]
        sample_snps = bin_snps[bin_snps.SAMPLE == sample]

        # compute phases if not multisample
        if not multisample_em:
            bafs[i], alphas[i], betas[i], _, haploflat, haplostring = apply_EM(
                totals[i].reshape(1, -1), alts[i].reshape(1, -1), reference_based_haplo_blocks, multisample_em
            )
            bafs[i] = bafs[i][0]
            alphas[i] = alphas[i][0]
            betas[i] = betas[i][0]

        assert np.sum(np.choose(haploflat, [sample_snps.ALT, sample_snps.REF])) == betas[i]

        snp_pos = ','.join(map(str, sample_snps.POS))
        snp_ref_counts = ','.join(map(str, sample_snps.REF))
        snp_alt_counts = ','.join(map(str, sample_snps.ALT))
        n_snps = len(sample_snps)

        baf = bafs[i]
        cov = (alphas[i] + betas[i]) / n_snps
        tot = alphas[i] + betas[i]

        result[sample] = (
            n_snps,
            cov,
            baf,
            alphas[i],
            betas[i],
            tot,
            snp_pos,
            snp_ref_counts,
            snp_alt_counts,
            haplostring,
        )
    return result


def merge_phasing(_, all_phase_data):
    """
    Merge phasing results across all samples:
    if a pair of SNPs is split in any sample, they won't be split.
    """

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

    orphans = orphans.union([a[0] for a in blocks if len(a) == 1])
    blocks = [b for b in blocks if len(b) > 1]

    return blocks, singletons, orphans


def binom_prop_test(alt1, ref1, flip1, alt2, ref2, flip2, alpha=0.1):
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
    denom = np.sqrt(p * (1 - p) * (1 / n1 + 1 / n2))
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
    df = df.copy().reset_index(drop=True)

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
    return np.split(data, np.where(np.diff(data) != stepsize)[0] + 1)


def phase_blocks_sequential(df, blocksize=50e3, max_snps_per_block=10, alpha=0):
    if len(df) == 0:
        return []

    """NOTE: this currently only supports single-sample data"""
    assert len(df.SAMPLE.unique()) == 1
    assert len(df.CHR.unique()) == 1

    df = df.copy()
    df = df.reset_index(drop=True)

    # Identify contiguous segments of SNPs with phasing information
    segments = consecutive(np.where(~df.PHASE.isna())[0])

    orphans = set()
    if alpha < 1:
        ### Use binomial test to split adjacent SNPs with very different phased BAFs
        # (further dividing the segments identified above)
        # identify indices of adj SNPs that are significantly different
        df_merged = pd.concat([df, df.shift(-1).add_prefix('next_')], axis=1)
        # no_merge is the left index of each such pair
        no_merge = np.where(
            df_merged[:-1].apply(
                lambda x: binom_prop_test(
                    x.ALT,
                    x.REF,
                    x.FLIP,
                    x.next_ALT,
                    x.next_REF,
                    x.next_FLIP,
                    alpha=alpha,
                ),
                axis=1,
            )
        )[0]
        # stack the adjacent indices (top is left, bottom is right = left + 1)
        nm = np.vstack([no_merge, no_merge + 1])

        new_segments = []
        for seg in segments:
            # identify adjacent pairs that violate test (p<alpha) and are both in this segment
            violations = np.all(np.isin(nm, seg), axis=0)

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
            mapping = {i: my_df.iloc[i].name for i in range(len(my_df))}
            blocks.extend([[mapping[a] for a in b] for b in my_blocks])

    blocks = sorted(blocks, key=lambda x: x[0])
    singletons = set(np.where(df.PHASE.isna())[0])

    return blocks, singletons, orphans


def collapse_blocks(df, blocks, singletons, orphans, ch):
    ### Construct blocked 1-bed table
    # First, merge blocks for those SNPs that are in the same block

    rows = []
    for sample, df0 in df.groupby('SAMPLE'):
        i = 0
        j = 0
        while i < len(df0):
            if i in singletons or i in orphans:
                r = df0.iloc[i]
                rows.append([ch, r.POS, r.POS, sample, r.ALT, r.REF, r.TOTAL, 1, [1]])
                i += 1
            else:
                block = blocks[j]
                assert i in block, (block, i)
                i += len(block)
                j += 1

                my_snps = df0.iloc[block]
                start = my_snps.POS.min()
                end = my_snps.POS.max()
                alt = np.sum(my_snps.FLIP * my_snps.REF + (1 - my_snps.FLIP) * my_snps.ALT).astype(np.uint64)
                ref = np.sum(my_snps.FLIP * my_snps.ALT + (1 - my_snps.FLIP) * my_snps.REF).astype(np.uint64)
                total = np.sum(my_snps.TOTAL)
                n_snps = len(my_snps)

                rows.append(
                    [
                        ch,
                        start,
                        end,
                        sample,
                        alt,
                        ref,
                        total,
                        n_snps,
                        list(my_snps.NOFLIP),
                    ]
                )

    return pd.DataFrame(
        rows,
        columns=[
            'CHR',
            'START',
            'END',
            'SAMPLE',
            'ALT',
            'REF',
            'TOTAL',
            '#SNPS',
            'HAPLO',
        ],
    ).sort_values(by=['CHR', 'START', 'SAMPLE'])


def merge_data(bins, dfs, bafs, sample_names, chromosome, unit):
    """
    Merge bins data (starts, ends, total counts, RDRs) with SNP data and BAF data for each bin.
    Parameters:
    bins: output from call to adaptive_bins_arm
    dfs: (only for troubleshooting) pandas DataFrame, each containing the SNP information for the corresponding bin
    bafs: the ith element is the output from compute_baf_task(dfs[i])

    Produces a BB file with a few additional columns.
    """

    rows = []
    if 'normal' in sample_names:
        sample_names = sample_names[1:]   # ignore the normal sample (first in the list)
        nonormalFlag = False
    else:
        nonormalFlag = True

    for i in range(len(bins[0])):
        start = bins[0][i]
        end = bins[1][i]
        totals = bins[2][i]
        rdrs = bins[3][i]

        if dfs is not None:
            assert dfs[i].POS.min() >= start and dfs[i].POS.max() <= end, (
                start,
                end,
                i,
                dfs[i].POS.min(),
                dfs[i].POS.max(),
            )
            snpcounts_from_df = dfs[i].pivot(index='POS', columns='SAMPLE', values='TOTAL').sum(axis=0)

        for j in range(len(sample_names)):
            sample = sample_names[j]
            if not nonormalFlag:
                total = int(totals[j + 1])
                normal_reads = int(totals[0])
            else:
                total = int(totals[j])
                normal_reads = 0
            rdr = rdrs[j]

            if dfs is not None:
                # the second variable below is "cov" from the old code. Should be removed.
                (nsnps, _, baf, alpha, beta, tot_snp, snp_pos, snp_ref_counts, snp_alt_counts, haplo,) = bafs[
                    i
                ][sample]
                assert snpcounts_from_df[sample] == alpha + beta, (i, sample)
            else:
                (
                    nsnps,
                    _,
                    baf,
                    alpha,
                    beta,
                    tot_snp,
                    snp_pos,
                    snp_ref_counts,
                    snp_alt_counts,
                    haplo,
                ) = (0, 0, 0, 0, 0, 0, '', '', '', '')

            rows.append(
                [
                    chromosome,
                    unit,
                    start,
                    end,
                    sample,
                    rdr,
                    total,
                    normal_reads,
                    nsnps,
                    beta,
                    tot_snp,
                    haplo,
                    snp_pos,
                    snp_ref_counts,
                    snp_alt_counts,
                    baf,
                ]
            )

    return pd.DataFrame(
        rows,
        columns=[
            'CHR',
            'UNIT',
            'START',
            'END',
            'SAMPLE',
            'RD',
            'TOTAL_READS',
            'NORMAL_READS',
            'SNPS',
            'BCOUNT',
            'TOTAL_SNP_READS',
            'HAPLO',
            'SNP_POS',
            'SNP_REF_COUNTS',
            'SNP_ALT_COUNTS',
            'BAF',
        ],
    )


def get_chr_end(stem, chromosome):
    fname = os.path.join(stem, chromosome + '.thresholds.gz')
    zcat = subprocess.Popen(('zcat', fname), stdout=subprocess.PIPE)
    tail = subprocess.Popen(('tail', '-1'), stdin=zcat.stdout, stdout=subprocess.PIPE)
    last_start = int(tail.stdout.read().decode('utf-8').strip())

    return last_start


def backtrack(bp):
    n, p = bp.shape

    starts = [bp[-1, -1]]
    for i in range(p - 1, 0, -1):
        starts.append(bp[starts[-1] - 1, i - 1])

    thresholds = starts[::-1] + [n]
    return thresholds


def segmented_piecewise(X, pieces=2):
    n, s = X.shape
    segcost_memo = {}

    def segment_cost(i, j):
        if (i, j) in segcost_memo:
            return segcost_memo[i, j]
        else:
            my_mean = np.mean(X[i:j], axis=0) if j > i else 0
            result = np.sum(np.square(X[i:j] - my_mean))
            segcost_memo[i, j] = result
            return result

    # penalty A[i, p] is the minimum error possible for fitting X[0] with p+1 pieces
    A = np.zeros((n, pieces))
    backpoint = np.zeros((n, pieces), dtype=int)

    A[:, 0] = [segment_cost(0, i) for i in range(n)]

    for p in range(1, pieces):
        for t in range(1, n):
            # search over t' < t
            best_cost = np.inf
            best_tprime = None
            for tprime in range(t + 1):
                # compute cost for making a new segment from tprime through t
                cost = A[tprime, p - 1] + segment_cost(tprime, t)
                if cost < best_cost:
                    best_cost = cost
                    best_tprime = tprime

            A[t, p] = best_cost
            backpoint[t, p] = best_tprime   # if this throws a TypeError for int(None), then X may have NaNs
    return A, backpoint


def correct_haplotypes(
    orig_bafs,
    min_prop_switch=0.01,
    n_segments=10,
    min_switch_density=0.1,
    min_mean_baf=0.45,
    minmax_al_imb=0.02,
):
    # Count switches using only samples with mean allelic imbalance above <minmax_al_imb>
    imb_samples = np.where(np.mean(np.abs(orig_bafs - 0.5), axis=0) > minmax_al_imb)[0]

    if len(imb_samples) == 0:
        sp.log(
            msg=f'No sample with avg. allelic imbalance above [{minmax_al_imb}], skipping correction.\n',
            level='INFO',
        )
        return orig_bafs, None

    bafs = orig_bafs[:, imb_samples]

    # Look for haplotype switches
    above_mid = bafs > 0.5
    is_alternating = np.concatenate(
        [
            np.zeros((1, bafs.shape[1]), dtype=bool),
            np.logical_xor(above_mid[1:], above_mid[:-1]),
        ]
    )
    haplotype_switches = np.where(np.all(is_alternating, axis=1))[0]
    prop_switched = len(haplotype_switches) / len(is_alternating)
    sp.log(
        msg=f'Found haplotype switches in {prop_switched *100 : .2f}% of bins\n',
        level='INFO',
    )

    if prop_switched > min_prop_switch:
        # If sufficient switches are found, run segmentation to identify segments w/ many switches

        sp.log(
            msg=f'Checking haplotype switching using [{n_segments}] segments\n',
            level='INFO',
        )
        # Segment using the mean BAFs only - faster and more reliable than using all samples
        A, bp = segmented_piecewise(np.mean(bafs, axis=1).reshape(-1, 1), pieces=n_segments)

        ts = backtrack(bp[:, :n_segments])

        for idx in np.where(np.diff(ts) == 0)[0]:
            del ts[idx]

        segments = [orig_bafs[ts[i] : ts[i + 1]] for i in range(len(ts) - 1)]

        # Identify problematic segments as those with many switches and mean near 0.5
        # (note that mean(BAF_i) across samples i is always <= 0.5 by def. from EM function)
        switch_densities = np.array(
            [
                len(haplotype_switches[(haplotype_switches >= ts[i]) & (haplotype_switches < ts[i + 1])])
                / (ts[i + 1] - ts[i])
                for i in range(len(ts) - 1)
            ]
        )
        segment_means = np.array([(np.mean(s) if len(s) > 0 else 0.5) for s in segments])

        # ALSO only correct segments with allelic imbalance at least <min_al_imb> in at least 1 sample
        segment_imbalances = np.array(
            [(np.max(np.abs(0.5 - np.mean(np.minimum(s, 1 - s), axis=0))) if len(s) > 0 else 0) for s in segments]
        )

        """
        segment_lengths = [len(s) for s in segments]
        [sp.log(msg=f'Segment {i}: length {a},\tmean {b:.3f},\timbalance {c:.3f},\tswitch prop. {d:.3f}\n',
                level = 'INFO')
               for i, (a,b,c, d)
               in enumerate(zip(segment_lengths, segment_means, segment_imbalances, switch_densities))]
        """

        bad_segments = np.where(
            np.logical_and(
                np.logical_and(
                    switch_densities >= min_switch_density,
                    segment_means >= min_mean_baf,
                ),
                segment_imbalances >= minmax_al_imb,
            )
        )[0]

        sp.log(
            msg=f'Identified {len(bad_segments)} segments with haplotype switching\n',
            level='INFO',
        )
        final_segments = []

        for s_idx in range(len(segments)):
            if s_idx in bad_segments:
                seg = segments[s_idx].copy()
                assert len(seg) > 0, 'Haplotype switch correction flagged a length-0 segment'

                # Identify the sample with the most extreme allelic imbalance for this segment
                minseg = np.minimum(seg, 1 - seg)
                mmeans = np.mean(minseg, axis=0)
                extreme_sample = np.argmin(mmeans)

                # Flip all of those bins with mhBAF > 0.5 in the corresponding sample
                my_flips = np.where(seg[:, extreme_sample] > 0.5)

                seg[my_flips] = 1 - seg[my_flips]

                my_mean = np.mean(seg)

                # If this ends up choosing the haplotype that is more abundant on average,
                if my_mean > 0.5:
                    # flip all bins to select the minor haplotype instead
                    seg = 1 - seg

                final_segments.append(seg)
            else:
                final_segments.append(segments[s_idx])

        return np.concatenate(final_segments), ts
    else:
        sp.log(
            msg=f'Insufficient haplotype switching detected (<[{min_prop_switch}]), skipping correction.\n',
            level='INFO',
        )
        return orig_bafs, None


def run_chromosome(
    baffile,
    all_names,
    chromosome,
    outfile,
    centromere_start,
    centromere_end,
    min_snp_reads,
    min_total_reads,
    arraystem,
    xy,
    multisample,
    phasefile,
    blocksize,
    max_snps_per_block,
    test_alpha,
    nonormalFlag,
    segfile,
):
    """
    Perform adaptive binning and infer BAFs to produce a HATCHet `bblock` and `btrack` file for a single chromosome.
    """

    try:
        if os.path.exists(outfile):
            sp.log(
                msg=f'Output file already exists, skipping chromosome {chromosome}\n',
                level='INFO',
            )
            return

        sp.log(
            msg=f'Loading intermediate files for chromosome {chromosome}\n',
            level='INFO',
        )
        if segfile:
            total_name, thresholds_name = (
                f'{chromosome}.segfile_total.gz',
                f'{chromosome}.segfile_thresholds.gz',
            )
        else:
            total_name, thresholds_name = (
                f'{chromosome}.total.gz',
                f'{chromosome}.thresholds.gz',
            )

        total_file = os.path.join(arraystem, total_name)
        thresholds_file = os.path.join(arraystem, thresholds_name)

        if not os.path.exists(total_file) or not os.path.exists(thresholds_file):
            raise ValueError(
                sp.error(
                    f'input files {total_file} or {thresholds_file} for custom segmentation not found!'
                    ' Make sure you have run both count-reads and combine-counts with or without the segmentation file'
                )
            )

        total_counts = np.loadtxt(total_file, dtype=np.uint32)
        complete_thresholds = np.loadtxt(thresholds_file, dtype=np.uint32)

        if segfile:
            sp.log(
                msg='# Collecting read depth anf BAF info for pre-specified segments!\n',
                level='STEP',
            )
            sp.log(
                msg=f'Reading SNPs file for chromosome {chromosome}\n',
                level='INFO',
            )
            # Load SNP positions and counts for this chromosome
            # snp_counts contains TOTAL depth for each sample at each site
            # snpsv is DataFrame with ALT/REF counts, phase info
            positions, snp_counts, snpsv = read_snps(baffile, chromosome, all_names, phasefile=phasefile)
            bins = bins_from_segfile(total_counts, complete_thresholds, positions)

            starts = bins[0]
            ends = bins[1]

            if xy:
                dfs = None
                bafs = None
            else:
                # Partition SNPs for BAF inference
                dfs = [snpsv[(snpsv.POS >= starts[i]) & (snpsv.POS <= ends[i])] for i in range(len(starts))]

                # if len(dfs) > 1:
                #    for i in range(len(dfs)):
                #        assert np.all(
                #            dfs[i].pivot(index='POS', columns='SAMPLE', values='TOTAL').sum(axis=0) >= min_snp_reads
                #        ), i

                # Infer BAF
                bafs = [
                    compute_baf_task(
                        d,
                        blocksize,
                        max_snps_per_block,
                        test_alpha,
                        multisample,
                    )
                    for d in dfs
                ]

            bb = merge_data(bins, dfs, bafs, all_names, chromosome)
            bb.to_csv(f'{outfile}.segfile', index=False, sep='\t')
            # np.savetxt(outfile + '.totalcounts', total_counts)
            # np.savetxt(outfile + '.thresholds', complete_thresholds)

            sp.log(
                msg=f'Done with custom segmentation on chromosome {chromosome}\n',
                level='INFO',
            )
        else:
            # adaptive binning
            # TODO: identify whether XX or XY, and only avoid SNPs/BAFs for XY
            if xy:
                sp.log(
                    msg='Running on sex chromosome -- ignoring SNPs \n',
                    level='INFO',
                )
                min_snp_reads = 0

                # TODO: do this procedure only for XY individuals
                ### construct dummy SNP positions and all-0 snpcounts array for binning
                before_centromere = complete_thresholds[complete_thresholds <= centromere_start]
                after_centromere = complete_thresholds[complete_thresholds >= centromere_end]
                positions_p = np.mean(
                    np.vstack([before_centromere[:-1], before_centromere[1:]]),
                    axis=0,
                ).astype(np.uint64)
                positions_q = np.mean(
                    np.vstack([after_centromere[:-1], after_centromere[1:]]),
                    axis=0,
                ).astype(np.uint64)
                positions = np.concatenate([positions_p, positions_q])
                snp_counts = np.zeros((len(positions), len(all_names) - 1), dtype=np.int8)
                snpsv = None

            else:
                # sp.log(msg=f"Reading SNPs file for chromosome {chromosome}\n", level = "INFO")
                # Load SNP positions and counts for this chromosome
                positions, snp_counts, snpsv = read_snps(baffile, chromosome, all_names, phasefile=phasefile)

            def arm_indices(arm):
                if arm == 'p':
                    # There may not be a SNP between the centromere end and the next SNP threshold
                    # Goal for p arm is to END at the FIRST threshold that is AFTER the LAST SNP BEFORE the centromere
                    if len(np.where(positions < centromere_start)[0]) > 0:
                        last_snp_before_centromere = positions[np.where(positions < centromere_start)[0][-1]]
                        last_threshold_before_centromere = complete_thresholds[
                            np.where(complete_thresholds > last_snp_before_centromere)[0][0]
                        ]
                        idx = np.where(complete_thresholds <= last_threshold_before_centromere)[0]
                        snp_idx = np.where(positions <= last_threshold_before_centromere)[0]
                    else:
                        idx, snp_idx = [], []
                else:   # arm == "q"
                    if len(np.where(positions > centromere_end)[0]) > 0:
                        first_snp_after_centromere = positions[np.where(positions > centromere_end)[0][0]]
                        first_threshold_after_centromere = complete_thresholds[
                            np.where(complete_thresholds < first_snp_after_centromere)[0][-1]
                        ]
                        idx = np.where(complete_thresholds >= first_threshold_after_centromere)[0]
                        snp_idx = np.where(positions >= first_threshold_after_centromere)[0]
                    else:
                        idx, snp_idx = [], []
                return idx, snp_idx

            armbbs = []
            for arm in ['p', 'q']:
                sp.log(
                    msg=f'Binning {arm} arm of chromosome {chromosome}\n',
                    level='INFO',
                )
                idx, snp_idx = arm_indices(arm)
                if len(snp_idx) > 0:
                    thresholds = complete_thresholds[idx]
                    counts = total_counts[idx]

                    arm_positions = positions[snp_idx]
                    arm_snp_counts = snp_counts[snp_idx]

                    # Identify bins
                    bins = adaptive_bins_arm(
                        snp_thresholds=thresholds,
                        total_counts=counts,
                        snp_positions=arm_positions,
                        snp_counts=arm_snp_counts,
                        min_snp_reads=min_snp_reads,
                        min_total_reads=min_total_reads,
                        nonormalFlag=nonormalFlag,
                    )

                    starts = bins[0]
                    ends = bins[1]
                    # Partition SNPs for BAF inference

                    # Infer BAF
                    if xy:
                        # TODO: compute BAFs for XX
                        dfs = None
                        bafs = None
                    else:
                        dfs = [snpsv[(snpsv.POS >= starts[i]) & (snpsv.POS <= ends[i])] for i in range(len(starts))]

                        if len(dfs) > 1:
                            for i in range(len(dfs)):
                                assert np.all(
                                    dfs[i]
                                    .pivot(
                                        index='POS',
                                        columns='SAMPLE',
                                        values='TOTAL',
                                    )
                                    .sum(axis=0)
                                    >= min_snp_reads
                                ), i

                        bafs = [
                            compute_baf_task(
                                d,
                                blocksize,
                                max_snps_per_block,
                                test_alpha,
                                multisample,
                            )
                            for d in dfs
                        ]
                    bb = merge_data(bins, dfs, bafs, all_names, chromosome, arm)

                    bafs = bb.pivot(index=['CHR', 'START'], columns='SAMPLE', values='BAF').to_numpy()
                    if bafs.shape[0] > 2:
                        sp.log(
                            msg=f'Correcting haplotype switches on {arm} arm...\n',
                            level='STEP',
                        )
                        # TODO: pass through other parameters to correct_haplotypes
                        corrected_bafs, _ = correct_haplotypes(bafs)

                        # flatten these results out and put them back into the BAF array
                        bb['ORIGINAL_BAF'] = bb.BAF
                        bb['BAF'] = corrected_bafs.flatten()
                else:
                    sp.log(
                        msg=f'No SNPs found in {arm} arm for {chromosome}\n',
                        level='INFO',
                    )
                    bb = None
                armbbs.append(bb)

            if all(element is None for element in armbbs):
                raise ValueError(sp.error(f'No SNPs found on either arm of chromosome {chromosome}!'))

            bb = pd.concat(armbbs)

            bb.to_csv(outfile, index=False, sep='\t', float_format='%.5f')
            # np.savetxt(outfile + '.totalcounts', total_counts)
            # np.savetxt(outfile + '.thresholds', complete_thresholds)

            sp.log(msg=f'Done chromosome {chromosome}\n', level='INFO')
    except Exception as e:
        sp.log(msg=f'Error in chromosome {chromosome}:', level='ERROR')
        sp.log(msg=str(e), level='ERROR')
        traceback.print_exc()
        raise e


def run_chromosome_wrapper(param):
    run_chromosome(*param)


def bins_from_segfile(total_counts, thresholds, snp_positions):
    n_samples = int(total_counts.shape[1] / 2)
    # number of reads that start between snp_thresholds[i] and snp_thresholds[i + 1]
    even_index = np.array([i * 2 for i in range(int(n_samples))], dtype=np.int8)
    # number of reads overlapping position snp_thresholds[i]
    odd_index = np.array([i * 2 + 1 for i in range(int(n_samples))], dtype=np.int8)
    bin_total = np.zeros(n_samples)
    # bin_snp = np.zeros(n_samples - 1)
    starts = []
    ends = []

    my_start = thresholds[0]

    rdrs = []
    totals = []
    i = 1
    j = 0
    while i < len(thresholds) - 1:
        # Extend the current bin to the next threshold
        next_threshold = thresholds[i]

        # add the intervening reads to the current bin
        # adding SNP reads
        # assert snp_positions[i - 1] <= thresholds[i], (
        #    i,
        #    snp_positions[i - 1],
        #    thresholds[i],
        # )
        # adding total reads
        bin_total += total_counts[i - 1, even_index]

        # get all SNPs
        snps_overlapping_interval = 0
        while (snp_positions[j] <= thresholds[i]) and (j < len(snp_positions) - 1):
            snps_overlapping_interval += 1
            j += 1

        # avoid creating bins without any BAF info
        if snps_overlapping_interval >= 1:
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
            my_start = ends[-1] + 1

        i += 1

    # handle the case of 1 bin
    if len(ends) == 0:
        sp.log(
            msg='WARNING: found only 1 bin in chromosome arm, may not meet MSR and MTR\t',
            level='WARN',
        )
        assert len(starts) == 0
        starts.append(thresholds[0])
        ends.append(thresholds[-1])

        bin_total = np.sum(total_counts[:, even_index], axis=0) - total_counts[-1, odd_index]
        totals.append(bin_total)
        rdrs.append(bin_total[1:] / bin_total[0])

    # add whatever excess at the end to the last bin
    if ends[-1] < thresholds[-1]:
        # combine the last complete bin with the remainder
        last_start_idx = np.where((thresholds == starts[-1] - 1) | (thresholds == starts[-1]))[0][0]
        bin_total = np.sum(total_counts[last_start_idx:, even_index], axis=0) - total_counts[-1, odd_index]
        totals[-1] = bin_total
        rdrs[-1] = bin_total[1:] / bin_total[0]
        ends[-1] = thresholds[-1]

    return starts, ends, totals, rdrs


if __name__ == '__main__':
    main()
