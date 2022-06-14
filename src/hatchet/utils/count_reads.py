import sys
import os
from multiprocessing import Pool
import subprocess as sp
from shutil import which

from datetime import datetime
import numpy as np
import pandas as pd
import gzip
import subprocess
import traceback
from importlib.resources import path
import hatchet.data

from hatchet.utils.ArgParsing import parse_count_reads_args
from hatchet.utils.Supporting import log, logArgs, error
import hatchet.utils.TotalCounting as tc


def main(args=None):
    log(msg='# Parsing and checking input arguments\n', level='STEP')
    args = parse_count_reads_args(args)
    logArgs(args, 80)

    bams = args['bams']
    names = args['names']

    tbams, tnames = zip(
        *sorted(zip(*(bams[1:], names[1:])), key=lambda x: x[1])
    )
    bams = [bams[0]] + list(tbams)
    names = [names[0]] + list(tnames)

    chromosomes = args['chromosomes']
    samtools = args['samtools']
    processes = args['j']
    outdir = args['outdir']
    mosdepth = args['mosdepth']
    tabix = args['tabix']
    readquality = args['readquality']

    if len(check_array_files(outdir, chromosomes)) == 0:
        log(
            msg='# Found all array files, skipping to total read counting. \n',
            level='STEP',
        )
    else:

        if len(check_counts_files(outdir, chromosomes, names)) == 0:
            log(
                msg='# Found all count files, skipping to forming arrays\n',
                level='STEP',
            )

        else:
            params = zip(
                np.repeat(chromosomes, len(bams)),
                [outdir] * len(bams) * len(chromosomes),
                [samtools] * len(bams) * len(chromosomes),
                bams * len(chromosomes),
                names * len(chromosomes),
                [readquality] * len(bams) * len(chromosomes),
            )

            n_workers_samtools = min(processes, len(bams) * len(chromosomes))
            with Pool(
                n_workers_samtools
            ) as p:   # divide by 2 because each worker starts 2 processes
                p.map(count_chromosome_wrapper, params)

            n_workers_mosdepth = min(processes, len(bams))

            # compute number of decompression threads to use for each call to mosdepth
            if processes > len(bams):
                base = min(int(processes / n_workers_mosdepth), 4)
                threads_per_task = [base] * len(bams)

                if base < 4:
                    remainder = processes % n_workers_mosdepth
                    i = 0
                    while remainder > 0:
                        threads_per_task[i] += 1
                        i += 1
                        remainder -= 1
            else:
                threads_per_task = [1] * len(bams)

            # Note: These function calls are the only section that uses mosdepth
            mosdepth_params = [
                (
                    outdir,
                    names[i],
                    bams[i],
                    threads_per_task[i],
                    mosdepth,
                    readquality,
                )
                for i in range(len(bams))
            ]
            with Pool(n_workers_mosdepth) as p:
                p.map(mosdepth_wrapper, mosdepth_params)

            if len(check_counts_files(outdir, chromosomes, names)) > 0:
                raise ValueError(error('Missing some counts files!'))

        ### Aggregate count files into count arrays for adaptive binning ###
        # (formerly formArray)
        use_chr = args['use_chr']

        n_workers = min(len(chromosomes), processes)

        # Read in centromere locations table
        with path(
            hatchet.data, f'{args["refversion"]}.centromeres.txt'
        ) as centromeres:
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
                    error(
                        f'Chromosome {ch} not found in centromeres file. Inspect file provided as -C argument.'
                    )
                )

        # Use Tabix to index per-position coverage bed files for each sample
        for name in names:
            perpos_file = os.path.join(outdir, name + '.per-base.bed.gz')
            subprocess.run([tabix, '-f', perpos_file])

        # form parameters for each worker
        params = [
            (
                outdir,
                names,
                ch,
                chr2centro[ch][0],
                chr2centro[ch][1],
                args['baf_file'],
                tabix,
            )
            for ch in chromosomes
        ]

        # dispatch workers
        with Pool(n_workers) as p:
            p.map(run_chromosome_wrapper, params)

        np.savetxt(os.path.join(outdir, 'samples.txt'), names, fmt='%s')

        if len(check_array_files(outdir, chromosomes)) > 0:
            raise ValueError(error('Missing some output arrays!'))
        else:
            log(
                msg='# Array forming completed successfully, removing intermediate count files. \n',
                level='STEP',
            )
            [
                os.remove(f)
                for f in expected_counts_files(outdir, chromosomes, names)
            ]

    totals_file = os.path.join(outdir, 'total.tsv')
    if os.path.exists(totals_file):
        log(msg='# Found total reads file, exiting. \n', level='STEP')
        return
    else:
        # TODO: take -q option and pass in here
        log(
            msg='# Counting total number of reads for normal and tumor samples\n',
            level='STEP',
        )
        total_counts = tc.tcount(
            samtools=samtools,
            samples=[(bams[i], names[i]) for i in range(len(names))],
            chromosomes=chromosomes,
            num_workers=processes,
            q=readquality,
        )

        try:
            total = {
                name: sum(
                    total_counts[name, chromosome]
                    for chromosome in chromosomes
                )
                for name in names
            }
        except KeyError:
            raise KeyError(
                error(
                    'Either a chromosome or a sample has not been considered in the total counting!'
                )
            )

        log(
            msg='# Writing the total read counts for all samples in {}\n'.format(
                totals_file
            ),
            level='STEP',
        )
        with open(totals_file, 'w') as f:
            for name in names:
                f.write('{}\t{}\n'.format(name, total[name]))


def mosdepth_wrapper(params):
    run_mosdepth(*params)


def run_mosdepth(outdir, sample_name, bam, threads, mosdepth, readquality):
    try:
        last_file = os.path.join(outdir, sample_name + '.per-base.bed.gz.csi')
        if os.path.exists(last_file):
            log(
                f'Skipping mosdepth on sample {sample_name} (output file {last_file} exists)\n',
                level='STEP',
            )
            return

        log(
            f'Starting mosdepth on sample {sample_name} with {threads} threads\n',
            level='STEP',
        )
        sys.stderr.flush()

        md = sp.run(
            [
                mosdepth,
                '-t',
                str(threads),
                '-Q',
                str(readquality),
                os.path.join(outdir, sample_name),
                bam,
            ]
        )
        md.check_returncode()
        log(f'Done mosdepth on sample {sample_name}\n', level='STEP')

    except Exception as e:
        log('Exception in countPos: {}\n'.format(e), level='ERROR')
        raise e


def count_chromosome(
    ch, outdir, samtools, bam, sample_name, readquality, compression_level=6
):
    try:
        outfile = os.path.join(outdir, f'{sample_name}.{ch}.starts')
        if os.path.exists(outfile):
            log(
                f'Skipping sample {sample_name} chromosome {ch} (output file exists)\n',
                level='STEP',
            )
            return

        log(
            f'Sample {sample_name} -- Starting chromosome {ch}\n', level='STEP'
        )
        sys.stderr.flush()

        # Get start positions
        st = sp.Popen(
            (samtools, 'view', '-q', str(readquality), bam, ch), stdout=sp.PIPE
        )
        cut = sp.Popen(('cut', '-f', '4'), stdin=st.stdout, stdout=sp.PIPE)
        gzip = sp.Popen(
            ('gzip', '-{}'.format(compression_level)),
            stdin=cut.stdout,
            stdout=open(outfile + '.gz', 'w'),
        )
        st.wait()
        cut.wait()
        gzip.wait()
        if st.returncode != 0:
            raise ValueError(
                'samtools subprocess returned nonzero value: {}'.format(
                    st.returncode
                )
            )
        if cut.returncode != 0:
            raise ValueError(
                'cut subprocess returned nonzero value: {}'.format(
                    cut.returncode
                )
            )

        log(f'Sample {sample_name} -- Done chromosome {ch}\n', level='STEP')

    except Exception as e:
        log('Exception in countPos: {}\n'.format(e), level='ERROR')
        raise e


def count_chromosome_wrapper(param):
    count_chromosome(*param)


def read_snps(baf_file, ch, all_names):
    """
    Read and validate SNP data for this patient (TSV table output from HATCHet deBAF.py).
    """
    all_names = all_names[
        1:
    ]   # remove normal sample -- not looking for SNP counts from normal

    # Read in HATCHet BAF table
    all_snps = pd.read_table(
        baf_file,
        names=['CHR', 'POS', 'SAMPLE', 'REF', 'ALT'],
        dtype={
            'CHR': object,
            'POS': np.uint32,
            'SAMPLE': object,
            'ALT': np.uint32,
            'REF': np.uint32,
        },
    )

    # Keep only SNPs on this chromosome
    snps = all_snps[all_snps.CHR == ch].sort_values(by=['POS', 'SAMPLE'])
    snps = snps.reset_index(drop=True)

    if len(snps) == 0:
        raise ValueError(
            error(
                f'Chromosome {ch} not found in SNPs file (chromosomes in file: {all_snps.CHR.unique()})'
            )
        )

    n_samples = len(all_names)
    if n_samples != len(snps.SAMPLE.unique()):
        raise ValueError(
            error(
                f'Expected {n_samples} samples, found {len(snps.SAMPLE.unique())} samples in SNPs file.'
            )
        )

    if set(all_names) != set(snps.SAMPLE.unique()):
        raise ValueError(
            error(
                f'Expected sample names did not match sample names in SNPs file.\n\
                Expected: {sorted(all_names)}\n  Found:{sorted(snps.SAMPLE.unique())}'
            )
        )

    # Add total counts column
    snpsv = snps.copy()
    snpsv['TOTAL'] = snpsv.ALT + snpsv.REF

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
    assert np.array_equal(all_names, list(snp_counts.columns)), (
        all_names,
        list(snp_counts.columns),
    )   # make sure that sample order is the same

    return np.array(snp_counts.index), np.array(snp_counts), snpsv


def form_counts_array(
    starts_files, perpos_files, thresholds, chromosome, tabix, chunksize=1e5
):
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
    arr = np.zeros(
        (thresholds.shape[0] + 1, len(starts_files) * 2)
    )   # add one for the end of the chromosome

    for i in range(len(starts_files)):
        # populate starts in even entries
        fname = starts_files[i]
        if fname.endswith('.gz'):
            f = gzip.open(fname)
        else:
            f = open(fname)

        read_starts = [int(a) for a in gzip.open(fname, 'r')]
        end = thresholds[0]
        j = 0
        for idx in range(1, len(thresholds)):
            start = end
            end = thresholds[idx]

            my_count = 0
            while j < len(read_starts) and read_starts[j] < end:
                assert read_starts[j] >= start
                my_count += 1
                j += 1

            arr[idx - 1, i * 2] = my_count

        arr[idx - 1, i * 2] += len(read_starts) - j
        f.close()

    for i in range(len(perpos_files)):
        # populate threshold coverage in odd entries
        fname = perpos_files[i]
        # print(datetime.now(), "Reading {}".format(fname))

        chr_sample_file = os.path.join(fname[:-3] + '.' + chromosome)

        if not os.path.exists(chr_sample_file):
            with open(chr_sample_file, 'w') as f:
                subprocess.run([tabix, fname, chromosome], stdout=f)

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


def get_chr_end(stem, all_names, chromosome):
    starts_files = []
    for name in all_names:
        starts_files.append(
            os.path.join(stem, name + '.' + chromosome + '.starts.gz')
        )

    last_start = 0
    for sfname in starts_files:
        zcat = subprocess.Popen(('zcat', sfname), stdout=subprocess.PIPE)
        tail = subprocess.Popen(
            ('tail', '-1'), stdin=zcat.stdout, stdout=subprocess.PIPE
        )
        my_last = int(tail.stdout.read().decode('utf-8').strip())

        if my_last > last_start:
            last_start = my_last

    return last_start


def run_chromosome(
    outdir,
    all_names,
    chromosome,
    centromere_start,
    centromere_end,
    baf_file,
    tabix,
):
    """
    Construct arrays that contain all counts needed to perform adaptive binning for a single chromosome (across all samples).
    """

    try:
        totals_out = os.path.join(outdir, f'{chromosome}.total.gz')
        thresholds_out = os.path.join(outdir, f'{chromosome}.thresholds.gz')

        if os.path.exists(totals_out) and os.path.exists(thresholds_out):
            log(
                msg=f'Output files already exist, skipping chromosome {chromosome}\n',
                level='INFO',
            )
            return

        log(msg=f'Loading chromosome {chromosome}\n', level='INFO')
        # Per-position coverage bed files for each sample
        perpos_files = [
            os.path.join(outdir, name + '.per-base.bed.gz')
            for name in all_names
        ]

        # Identify the start-positions files for this chromosome
        starts_files = []
        for name in all_names:
            starts_files.append(
                os.path.join(outdir, name + '.' + chromosome + '.starts.gz')
            )

        log(
            msg=f'Reading SNPs file for chromosome {chromosome}\n',
            level='INFO',
        )
        # Load SNP positions and counts for this chromosome

        if chromosome.endswith('X') or chromosome.endswith('Y'):
            log(
                msg='Running on sex chromosome -- ignoring SNPs and min SNP reads\n',
                level='INFO',
            )

            # TODO: do this procedure only for XY
            last_start = get_chr_end(outdir, all_names, chromosome)
            positions = np.arange(5000, last_start, 5000)

        else:
            positions, _, _ = read_snps(baf_file, chromosome, all_names)

        thresholds = np.trunc(
            np.vstack([positions[:-1], positions[1:]]).mean(axis=0)
        ).astype(np.uint32)
        last_idx_p = np.argwhere(thresholds > centromere_start)[0][0]
        first_idx_q = np.argwhere(thresholds > centromere_end)[0][0]
        all_thresholds = np.concatenate(
            [
                [1],
                thresholds[:last_idx_p],
                [centromere_start],
                [centromere_end],
                thresholds[first_idx_q:],
            ]
        )

        log(msg=f'Loading counts for chromosome {chromosome}\n', level='INFO')
        # Load read count arrays from file (also adds end of chromosome as a threshold)
        total_counts, complete_thresholds = form_counts_array(
            starts_files, perpos_files, all_thresholds, chromosome, tabix=tabix
        )

        np.savetxt(totals_out, total_counts, fmt='%d')
        np.savetxt(thresholds_out, complete_thresholds, fmt='%d')

        log(msg=f'Done chromosome {chromosome}\n', level='INFO')
    except Exception as e:
        print(f'Error in chromosome {chromosome}:')
        print(e)
        traceback.print_exc()
        raise e


def run_chromosome_wrapper(param):
    run_chromosome(*param)


def expected_counts_files(dcounts, chrs, all_names):
    expected = []

    for chr in chrs:
        for name in all_names:
            fname = os.path.join(dcounts, '.'.join([name, chr, 'starts.gz']))
            expected.append(fname)

    # mosdepth output
    for name in all_names:
        fname = os.path.join(dcounts, name + '.mosdepth.global.dist.txt')
        expected.append(fname)

        fname = os.path.join(dcounts, name + '.mosdepth.summary.txt')
        expected.append(fname)

        fname = os.path.join(dcounts, name + '.per-base.bed.gz')
        expected.append(fname)

        fname = os.path.join(dcounts, name + '.per-base.bed.gz.csi')
        expected.append(fname)

    return expected


def check_counts_files(dcounts, chrs, all_names):
    return [
        a
        for a in expected_counts_files(dcounts, chrs, all_names)
        if not os.path.exists(a)
    ]


def expected_arrays(darray, chrs):
    expected = []
    # formArray (abin/<chr>.<total/thresholds> files))

    fname = os.path.join(darray, 'samples.txt')
    expected.append(fname)

    for ch in chrs:
        fname = os.path.join(darray, f'{ch}.total.gz')
        expected.append(fname)

        fname = os.path.join(darray, f'{ch}.thresholds.gz')
        expected.append(fname)

    return expected


def check_array_files(darray, chrs):
    return [a for a in expected_arrays(darray, chrs) if not os.path.exists(a)]


if __name__ == '__main__':
    main()
