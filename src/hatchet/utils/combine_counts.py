from multiprocessing import Pool
import os
import subprocess
import traceback
from importlib.resources import path
import numpy as np
import pandas as pd

import hatchet.data
from hatchet.utils.ArgParsing import parse_combine_counts_args
import hatchet.utils.Supporting as sp
from hatchet.utils.pon_normalization import correct_baf, pon_normalize_rdr
from hatchet.utils.rd_gccorrect import rd_gccorrect
from hatchet.utils.compute_baf import compute_baf_task
from hatchet.utils.correct_haplotypes import correct_haplotypes


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
    XX = args['XX']

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

    # check if isX true for any of the chromosomes
    if any(isX.values()):
        sp.log(
            msg='Determining the allosomes for handling X chromosome\n',
            level='INFO',
        )
        # X chromosome is present. We must determine sex of samples if not provided already.
        if XX == 'auto':
            # if read counts for Y chromosome is not provided, we cannot determine sex.
            # We assume the sample is from a female
            if not any(isY.values()):
                xy = False
            else:
                # find total read counts for X and Y chromosomes.
                # if the ratio of X to Y is larger than 50, we assume the sample is from a female

                # find the first chromosome such that isX is true
                x_ch = next((ch for ch in chromosomes if isX[ch]), None)
                y_ch = next((ch for ch in chromosomes if isY[ch]), None)

                x_tc, _ = read_total_and_thresholds(x_ch, args['array'], args['segfile'])
                y_tc, _ = read_total_and_thresholds(y_ch, args['array'], args['segfile'])
                total_x_reads = x_tc[:, 0].sum()
                total_y_reads = y_tc[:, 0].sum()
                if total_y_reads == 0:
                    xy = False
                else:
                    xy = total_x_reads / total_y_reads <= 50

        elif XX == 'True':
            xy = False
        elif XX == 'False':
            xy = True

        sp.log(
            msg='Allosomes are determined as ' + ('XY' if xy else 'XX') + '\n',
            level='INFO',
        )
    else:
        # there is no X chromosome, no need to determine sex
        # it does not matter what sex is given
        xy = False
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
            xy,
            multisample,
            phase,
            blocksize,
            max_snps_per_block,
            test_alpha,
            nonormalFlag,
            args['segfile'],
        )
        for ch in chromosomes if not ch.endswith('Y') and not (ch.endswith('X') and xy)
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
                df['TOTAL_READS'] / (df['END'] - df['START']) / np.mean(df['TOTAL_READS'] / (df['END'] - df['START']))
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

    # autosomes = set([ch for ch in big_bb['CHR'] if not (ch.endswith('X') or ch.endswith('Y'))])
    # autosomal_bb = big_bb[big_bb['CHR'].isin(autosomes)].copy()
    # autosomal_bb = rd_gccorrect(autosomal_bb, referencefasta)
    big_bb = rd_gccorrect(big_bb, referencefasta)

    if xy and any(isX.values()):
        x_ch = next((ch for ch in chromosomes if isX[ch]), None)
        # set BAF to 0.5 for X chromosome
        # also double the RD
        big_bb.loc[big_bb.CHR == x_ch, 'BAF'] = 0.5

    # Convert intervals from closed to half-open to match .1bed/HATCHet standard format
    # autosomal_bb.END = autosomal_bb.END + 1
    # autosomal_bb.to_csv(outfile, index=False, sep='\t')

    big_bb.END = big_bb.END + 1
    big_bb.to_csv(outfile, index=False, sep='\t')

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
    chromosome,
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
    assert len(snp_positions) == len(snp_thresholds) - 1 or chromosome.endswith('X'), (
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
    j = 0
    while i < len(snp_thresholds - 1) and j < len(snp_positions):
        # Extend the current bin to the next threshold
        next_threshold = snp_thresholds[i]
        if next_threshold < snp_positions[j]:
            # adding total reads
            bin_total += total_counts[i - 1, even_index]
            i += 1
            continue
        # add the intervening reads to the current bin
        # adding SNP reads
        while j < len(snp_positions) and snp_positions[j] <= next_threshold:
            bin_snp += snp_counts[j]
            j += 1
        assert snp_positions[j - 1] >= snp_thresholds[i - 1]
        assert snp_positions[j - 1] <= snp_thresholds[i], (
            i,
            snp_positions[j - 1],
            snp_thresholds[i],
        )

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
        total_counts, complete_thresholds = read_total_and_thresholds(chromosome, arraystem, segfile)

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

            if xy and chromosome.endswith('X'):
                # the sample is from a male. Don't compute BAF on X chromosome
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
            if xy and chromosome.endswith('X'):
                sp.log(
                    msg='Running on sex chromosome X in a male -- ignoring SNPs \n',
                    level='INFO',
                )
                min_snp_reads = 0

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
                        chromosome=chromosome,
                        min_snp_reads=min_snp_reads,
                        min_total_reads=min_total_reads,
                        nonormalFlag=nonormalFlag,
                    )

                    starts = bins[0]
                    ends = bins[1]
                    # Partition SNPs for BAF inference

                    # Infer BAF
                    if xy and chromosome.endswith('X'):
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


def read_total_and_thresholds(chromosome, arraystem, segfile):
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
    return total_counts, complete_thresholds


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
