import sys
import os.path
from os.path import isfile, isdir
import argparse
import subprocess

from hatchet.utils.Supporting import (
    ensure,
    log,
    error,
    which,
    url_exists,
    bcolors,
    numericOrder,
)
from hatchet import config, __version__


def parse_plot_bins_1d2d_args(args=None):
    """
    Parse command line arguments for auxiliary cluster plotting command (1D and 2D plot with matching colors and
    optional labeled centers)
    """
    description = ''
    parser = argparse.ArgumentParser(
        prog='hatchet plot-cn',
        description=description,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        '-b',
        '--bbc',
        help='Filename for BBC table (e.g., bbc/bulk.bbc)',
        required=True,
        type=str,
    )
    parser.add_argument(
        '-s',
        '--seg',
        help='Filename for SEG table (e.g., bbc/bulk.seg) optional but required to show cluster centers',
        required=False,
        type=str,
        default=config.plot_bins_1d2d.seg,
    )
    parser.add_argument(
        '-O',
        '--outdir',
        required=True,
        type=str,
        help='Directory for output files',
    )
    parser.add_argument(
        '--baflim',
        required=False,
        type=str,
        help=(
            "Axis limits for mirrored BAF as comma-separated values, e.g., '0,0.51' (default: None -- show full "
            'range of data)'
        ),
        default=config.plot_bins_1d2d.baflim,
    )
    parser.add_argument(
        '--rdrlim',
        required=False,
        type=str,
        help=(
            "Axis limits for read-depth ratio as comma-separated values, e.g., '0,3' (default: None -- show full "
            'range of data)'
        ),
        default=config.plot_bins_1d2d.fcnlim,
    )
    parser.add_argument(
        '--centers',
        required=False,
        action='store_true',
        help='Show cluster centers (requires SEG file provided via -s/--seg)',
        default=config.plot_bins_1d2d.centers,
    )
    parser.add_argument(
        '--centromeres',
        required=False,
        action='store_true',
        help='Mark centromere locations with grey rectangles',
        default=config.plot_bins_1d2d.centromeres,
    )
    parser.add_argument(
        '-a',
        '--alpha',
        required=False,
        type=float,
        help='Opacity alpha (default 1)',
        default=config.plot_bins_1d2d.alpha,
    )
    args = parser.parse_args(args)

    ensure(isfile(args.bbc), f'Input BBC file [{args.bbc}] not found.')
    ensure(
        (args.seg is None) or (isfile(args.seg)),
        f'Input SEG file [{args.seg}] not found.',
    )

    if args.baflim is not None:
        tkns = args.baflim.split(',')
        try:
            minbaf = float(tkns[0].strip())
            maxbaf = float(tkns[1].strip())
        except ValueError:
            error(
                '--baflim must be comma-separated float or integer values.',
                raise_exception=True,
            )
        ensure(minbaf < maxbaf, 'Minimum BAF must be < maximum BAF.')
    else:
        minbaf = maxbaf = None

    if args.rdrlim is not None:
        tkns = args.rdrlim.split(',')
        try:
            minrdr = float(tkns[0].strip())
            maxrdr = float(tkns[1].strip())
        except ValueError:
            error(
                '--rdrlim must be comma-separated float or integer values.',
                raise_exception=True,
            )
        ensure(minrdr < maxrdr, 'Minimum RDR must be < maximum RDR.')
    else:
        minrdr = maxrdr = None

    return {
        'bbc': args.bbc,
        'seg': args.seg,
        'outdir': args.outdir,
        'minbaf': minbaf,
        'maxbaf': maxbaf,
        'minrdr': minrdr,
        'maxrdr': maxrdr,
        'centers': args.centers,
        'centromeres': args.centromeres,
        'alpha': args.alpha,
    }


def parse_plot_cn_1d2d_args(args=None):
    """
    Parse command line arguments for auxiliary plotting command (1D and 2D plot with labeled copy states)
    """
    description = ''
    parser = argparse.ArgumentParser(
        prog='hatchet plot-cn',
        description=description,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument('INPUT', help='Filename for BBC table (e.g., results/best.bbc.ucn)')
    parser.add_argument(
        '-O',
        '--outdir',
        required=True,
        type=str,
        help='Directory for output files',
    )
    parser.add_argument(
        '--baflim',
        required=False,
        type=str,
        help=(
            "Axis limits for mirrored BAF values to show as comma-separated values, e.g., '0,0.51' (default: None "
            '-- show full range of data)'
        ),
        default=config.plot_cn_1d2d.baflim,
    )
    parser.add_argument(
        '--fcnlim',
        required=False,
        type=str,
        help=(
            "Axis limits for fractional copy number values to show as comma-separated values, e.g., '0,3' (default: "
            'None -- show full range of data)'
        ),
        default=config.plot_cn_1d2d.fcnlim,
    )
    parser.add_argument(
        '--centromeres',
        required=False,
        action='store_true',
        help='Mark centromere locations with grey rectangles',
        default=config.plot_cn_1d2d.centromeres,
    )
    parser.add_argument(
        '--bysample',
        required=False,
        action='store_true',
        help='Write each sample to a separate file rather than combining all into 2 file',
        default=config.plot_cn_1d2d.bysample,
    )

    args = parser.parse_args(args)

    ensure(isfile(args.INPUT), f'Input file [{args.INPUT}] not found.')

    if args.baflim is not None:
        tkns = args.baflim.split(',')
        try:
            minbaf = float(tkns[0].strip())
            maxbaf = float(tkns[1].strip())
        except ValueError:
            error(
                '--baflim must be comma-separated float or integer values.',
                raise_exception=True,
            )
    else:
        minbaf = maxbaf = None

    if args.fcnlim is not None:
        tkns = args.fcnlim.split(',')
        try:
            minfcn = float(tkns[0].strip())
            maxfcn = float(tkns[1].strip())
        except ValueError:
            error(
                '--fcnlim must be comma-separated float or integer values.',
                raise_exception=True,
            )
    else:
        minfcn = maxfcn = None

    return {
        'input': args.INPUT,
        'outdir': args.outdir,
        'minbaf': minbaf,
        'maxbaf': maxbaf,
        'minfcn': minfcn,
        'maxfcn': maxfcn,
        'bysample': args.bysample,
        'centromeres': args.centromeres,
    }


def parse_cluster_bins_args(args=None):
    """
    Parse command line arguments
    Returns:
    """
    description = (
        'Combine tumor bin counts, normal bin counts, and tumor allele counts to obtain the read-depth '
        'ratio and the mean B-allel frequency (BAF) of each bin. Optionally, the normal allele counts can '
        'be provided to add the BAF of each bin scaled by the normal BAF. The output is written on stdout.'
    )
    parser = argparse.ArgumentParser(
        prog='hatchet cluster-bins',
        description=description,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        'BBFILE',
        help=(
            'A BB file containing a line for each bin in each sample and the corresponding values of read-depth '
            'ratio and B-allele frequency (BAF)'
        ),
    )
    parser.add_argument(
        '-o',
        '--outsegments',
        required=False,
        default=config.cluster_bins.outsegments,
        type=str,
        help='Output filename for the segments computed by clustering bins (default: stdout)',
    )
    parser.add_argument(
        '-O',
        '--outbins',
        required=False,
        default=config.cluster_bins.outbins,
        type=str,
        help='Output filename for a BB file adding the clusters (default: stdout)',
    )
    parser.add_argument(
        '-d',
        '--diploidbaf',
        type=float,
        required=False,
        default=config.cluster_bins.diploidbaf,
        help=(
            'Maximum diploid-BAF shift used to determine the largest copy-neutral cluster and to rescale all the '
            'cluster inside this threshold accordingly (default: 0.1)'
        ),
    )
    parser.add_argument(
        '--minK',
        type=int,
        required=False,
        default=config.cluster_bins.mink,
        help=f'Minimum number of clusters to infer (default = {config.cluster_bins.mink})',
    )
    parser.add_argument(
        '--maxK',
        type=int,
        required=False,
        default=config.cluster_bins.maxk,
        help=f'Maximum number of clusters to infer (default = {config.cluster_bins.maxk})',
    )
    parser.add_argument(
        '--exactK',
        type=int,
        required=False,
        default=config.cluster_bins.exactk,
        help='Skip model selection and infer exactly this many clusters (default: None)',
    )
    parser.add_argument(
        '-t',
        '--transmat',
        type=str,
        required=False,
        default=config.cluster_bins.transmat,
        help='Form of transition matrix to infer: fixed, diag (1-parameter), or full (default: diag)',
    )
    parser.add_argument(
        '--tau',
        type=float,
        required=False,
        default=config.cluster_bins.tau,
        help=f'Off-diagonal value for initializing transition matrix (default: {config.cluster_bins.tau})',
    )
    parser.add_argument(
        '-c',
        '--covar',
        type=str,
        required=False,
        default=config.cluster_bins.covar,
        help='Form of covariance matrix: spherical, diag, full, or tied (default: diag)',
    )
    parser.add_argument(
        '-x',
        '--decoding',
        type=str,
        required=False,
        default=config.cluster_bins.decoding,
        help='Decoding algorithm to use: map or viterbi (default: map)',
    )
    parser.add_argument(
        '-R',
        '--restarts',
        type=int,
        required=False,
        default=config.cluster_bins.restarts,
        help='Number of restarts performed by the clustering to choose the best (default: 10)',
    )
    parser.add_argument(
        '-S',
        '--subset',
        required=False,
        default=config.cluster_bins.subset,
        type=str,
        nargs='+',
        help=(
            'List of sample names to use as a subset of those included in binning' '(default: none, run on all samples)'
        ),
    )
    parser.add_argument('-V', '--version', action='version', version=f'%(prog)s {__version__}')
    args = parser.parse_args(args)

    ensure(isfile(args.BBFILE), 'The specified BB file does not exist!')
    ensure(
        (args.diploidbaf is None) or (0.0 <= args.diploidbaf <= 0.5),
        'The specified maximum for diploid-BAF shift must be a value in [0.0, 0.5]',
    )

    if args.exactK > 0:
        log(
            msg=f'WARNING: using fixed K (K={args.exactK}) will bypass model selection.\n',
            level='WARN',
        )

    ensure(args.minK >= 0, 'Min. number of clusters must be non-negative.')
    ensure(args.maxK >= 0, 'Max. number of clusters must be non-negative.')
    ensure(
        args.minK <= args.maxK,
        'Min number of clusters must be smaller than max number of clusters.',
    )

    if args.maxK == args.minK:
        log(
            msg=f'WARNING: minK = maxK = {args.minK}, so no model selection will be performed.\n',
            level='WARN',
        )

    valid_covars = ('spherical', 'diag', 'full', 'tied')
    ensure(
        args.covar in valid_covars,
        f'Invalid -c/--covar argument: {args.covar}. Valid values are {valid_covars}.',
    )

    valid_decodings = ('map', 'viterbi')
    ensure(
        args.decoding in valid_decodings,
        f'Invalid -x/--decoding argument: {args.decoding}. Valid values are {valid_decodings}.',
    )

    valid_transmats = ('fixed', 'diag', 'tied')
    ensure(
        args.transmat in valid_transmats,
        f'Invalid -t/--transmat argument: {args.transmat}. Valid values are {valid_transmats}.',
    )
    ensure(args.restarts >= 1, 'Number of restarts must be positive.')
    ensure(args.tau >= 0, 'Transition parameter --tau must be non-negative.')
    ensure(args.tau >= 6e-17, 'Transition parameter --tau must be at least 6e-17 to ensure numerical precision.')

    if args.subset is not None:
        import pandas as pd

        args.subset = args.subset.split()
        bb = pd.read_table(args.BBFILE)
        samples = bb.SAMPLE.unique()
        ensure(
            all([a in samples for a in args.subset]),
            'Samples indicated in "subset" must be present in input BB file. BB file:'
            + f'{samples}, argument: {args.subset}',
        )

    return {
        'bbfile': args.BBFILE,
        'diploidbaf': args.diploidbaf,
        'decoding': args.decoding,
        'minK': args.minK,
        'maxK': args.maxK,
        'exactK': args.exactK,
        'covar': args.covar,
        'transmat': args.transmat,
        'tau': args.tau,
        'outsegments': args.outsegments,
        'outbins': args.outbins,
        'restarts': args.restarts,
        'subset': args.subset,
    }


def parse_count_reads_args(args=None):
    description = (
        'Count the mapped sequencing reads in bins of fixed and given length, uniformly for a BAM file of '
        'a normal sample and one or more BAM files of tumor samples. This program supports both data from '
        'whole-genome sequencing (WGS) and whole-exome sequencing (WES), but the a BED file with targeted '
        'regions is required when considering WES.'
    )
    parser = argparse.ArgumentParser(prog='hatchet count-reads', description=description)
    parser.add_argument(
        '-N',
        '--normal',
        required=True,
        type=str,
        help='BAM file corresponding to matched normal sample',
    )
    parser.add_argument(
        '-T',
        '--tumor',
        required=True,
        type=str,
        nargs='+',
        help='BAM files corresponding to samples from the same tumor',
    )
    parser.add_argument(
        '-b',
        '--baffile',
        required=True,
        type=str,
        help='1bed file containing SNP information from tumor samples (i.e., baf/bulk.1bed)',
    )
    parser.add_argument(
        '-V',
        '--refversion',
        required=True,
        type=str,
        help='Version of reference genome used in BAM files',
    )
    parser.add_argument(
        '-O',
        '--outdir',
        required=True,
        type=str,
        help='Directory for output files',
    )
    parser.add_argument(
        '-S',
        '--samples',
        required=False,
        default=config.count_reads.samples,
        type=str,
        nargs='+',
        help=(
            'Sample names for each BAM, given in the same order where the normal name is first (default: inferred '
            'from file names)'
        ),
    )
    parser.add_argument(
        '-st',
        '--samtools',
        required=False,
        type=str,
        default=config.paths.samtools,
        help='Path to samtools executable',
    )
    parser.add_argument(
        '-md',
        '--mosdepth',
        required=False,
        type=str,
        default=config.paths.mosdepth,
        help='Path to mosdepth executable',
    )
    parser.add_argument(
        '-tx',
        '--tabix',
        required=False,
        type=str,
        default=config.paths.tabix,
        help='Path to tabix executable',
    )
    parser.add_argument(
        '-j',
        '--processes',
        required=False,
        default=config.count_reads.processes,
        type=int,
        help='Number of available parallel processes (default: 2)',
    )
    parser.add_argument(
        '-q',
        '--readquality',
        required=False,
        default=config.count_reads.readquality,
        type=int,
        help=(
            'Minimum mapping quality for an aligned read to be considered (default: ',
            f'{config.count_reads.readquality})',
        ),
    )
    parser.add_argument(
        '-i',
        '--intermediates',
        required=False,
        action='store_true',
        default=config.count_reads.intermediates,
        help='Produce intermediate counts files only and do not proceed to forming arrays (default: False)',
    )
    parser.add_argument(
        '--chromosomes',
        required=False,
        nargs='*',
        help='One or more chromosomes to process (default: blank to process all chromosomes)',
    )

    args = parser.parse_args(args)

    # Parse BAM files, check their existence, and infer or parse the corresponding sample names
    bams = [args.normal] + args.tumor
    for bamfile in bams:
        ensure(isfile(bamfile), 'The specified tumor BAM file does not exist')

    names = args.samples
    ensure(
        (names is None) or (len(bams) == len(names)),
        'A sample name must be provided for each corresponding BAM: both for each normal sample and each tumor sample',
    )

    if names is None:
        names = ['normal']
        for bamfile in bams[1:]:
            names.append(os.path.splitext(os.path.basename(bamfile))[0])

    # In default mode, check the existence and compatibility of samtools and bcftools
    samtools = os.path.join(args.samtools, 'samtools')
    ensure(
        which(samtools) is not None,
        'samtools has not been found or is not executable.',
    )

    # Extract the names of the chromosomes and check their consistency across the given BAM files and the reference
    chromosomes = extractChromosomes(samtools, [bams[0], 'normal'], [(x, '') for x in bams[1:]])
    if args.chromosomes:
        chromosomes = [c for c in chromosomes if c in args.chromosomes]

    # Check that chr notation is consistent across chromosomes
    using_chr = [a.startswith('chr') for a in chromosomes]
    if any(using_chr):
        ensure(
            all(using_chr),
            'Some chromosomes use chr notation while others do not.',
        )
        use_chr = True
    else:
        use_chr = False

    ver = args.refversion
    ensure(
        ver in ('hg19', 'hg38'),
        'Invalid reference genome version. Supported versions are hg38 and hg19.',
    )
    ensure(os.path.exists(args.baffile), f'BAF file not found: {args.baffile}')
    ensure(
        args.processes > 0,
        'The number of parallel processes must be greater than 0',
    )
    ensure(
        os.path.exists(args.outdir),
        f'Output directory <{args.outdir}> not found!',
    )

    mosdepth = os.path.join(args.mosdepth, 'mosdepth')
    tabix = os.path.join(args.tabix, 'tabix')
    ensure(
        which(mosdepth) is not None,
        (
            'The mosdepth executable was not found or is not executable. Please install mosdepth (e.g., conda install '
            '-c bioconda mosdepth) and/or supply the path to the executable.'
        ),
    )
    ensure(
        which(tabix) is not None,
        (
            'The tabix executable was not found or is not executable. Please install tabix (e.g., conda install -c '
            'bioconda tabix) and/or supply the path to the executable.'
        ),
    )

    return {
        'bams': bams,
        'names': names,
        'chromosomes': chromosomes,
        'samtools': samtools,
        'mosdepth': mosdepth,
        'tabix': tabix,
        'j': args.processes,
        'outdir': args.outdir,
        'use_chr': use_chr,
        'refversion': ver,
        'baf_file': args.baffile,
        'readquality': args.readquality,
    }


def parse_combine_counts_args(args=None):
    parser = argparse.ArgumentParser(
        description='Perform adaptive binning, compute RDR and BAF for each bin, and produce a BB file.'
    )
    parser.add_argument(
        '-A',
        '--array',
        type=str,
        required=True,
        help='Directory containing array files (output from "count_reads" command)',
    )
    parser.add_argument(
        '-t',
        '--totalcounts',
        required=True,
        type=str,
        help=(
            'Total read counts in the format "SAMPLE\tCOUNT" used to normalize by the different number of reads '
            'extracted from each sample'
        ),
    )
    parser.add_argument(
        '-b',
        '--baffile',
        required=True,
        type=str,
        help='1bed file containing SNP information from tumor samples (i.e., baf/bulk.1bed)',
    )
    parser.add_argument('-o', '--outfile', required=True, type=str, help='Filename for output')
    parser.add_argument(
        '--msr',
        type=int,
        help=f'Minimum SNP reads per bin (default {config.combine_counts.msr})',
        default=config.combine_counts.msr,
    )
    parser.add_argument(
        '--mtr',
        type=int,
        help=f'Minimum total reads per bin (default {config.combine_counts.mtr})',
        default=config.combine_counts.mtr,
    )
    parser.add_argument(
        '-j',
        '--processes',
        type=int,
        help=f'Number of parallel processes to use (default {config.combine_counts.processes})',
        default=config.combine_counts.processes,
    )
    parser.add_argument(
        '-p',
        '--phase',
        required=False,
        default=config.combine_counts.phase,
        type=str,
        help='VCF file containing phasing for heterozygous germline SNPs',
    )
    parser.add_argument(
        '-s',
        '--max_blocksize',
        required=False,
        default=config.combine_counts.blocksize,
        type=int,
        help=f'Maximum size of phasing block (default {config.combine_counts.blocksize})',
    )
    parser.add_argument(
        '-m',
        '--max_spb',
        required=False,
        default=config.combine_counts.max_spb,
        type=int,
        help=f'Maximum number of SNPs per phasing block (default {config.combine_counts.max_spb})',
    )
    parser.add_argument(
        '-a',
        '--alpha',
        required=False,
        default=config.combine_counts.alpha,
        type=float,
        help=(
            'Significance level for phase blocking adjacent SNPs. Higher means less trust in phasing. '
            f'(default {config.combine_counts.alpha})'
        ),
    )
    parser.add_argument(
        '--ss_em',
        action='store_true',
        default=config.combine_counts.ss_em,
        required=False,
        help='Use single-sample EM BAF inference (instead of multi-sample)',
    )
    parser.add_argument(
        '-z',
        '--not_compressed',
        action='store_true',
        default=config.combine_counts.not_compressed,
        required=False,
        help='Non-compressed intermediate files',
    )
    parser.add_argument(
        '-V',
        '--refversion',
        required=True,
        type=str,
        help='Version of reference genome used in BAM files',
    )
    args = parser.parse_args(args)

    ensure(os.path.exists(args.baffile), f'BAF file not found: {args.baffile}')

    if args.totalcounts is not None and not isfile(args.totalcounts):
        raise ValueError(error('The specified file for total read counts does not exist!'))
    ensure(
        (args.phase is None) or isfile(args.phase),
        'The specified phasing file does not exist!',
    )
    ensure(args.max_blocksize > 0, 'The max_blocksize argument must be positive.')
    ensure(args.max_spb > 0, 'The max_snps_per_bin argument must be positive.')
    ensure(
        0 <= args.alpha <= 1,
        'The alpha argument must be between 0 and 1, inclusive.',
    )
    ensure(
        os.path.exists(args.array),
        f'The provided array directory does not exist: {args.array}',
    )

    namesfile = os.path.join(args.array, 'samples.txt')
    ensure(
        os.path.exists(namesfile),
        'Missing file containing sample names (1 per line): {namesfile}',
    )
    names = open(namesfile).read().split()

    chromosomes = set()
    for f in os.listdir(args.array):
        tkns = f.split('.')
        if len(tkns) > 1 and (tkns[1] == 'thresholds' or tkns[1] == 'total'):
            chromosomes.add(tkns[0])
    chromosomes = sorted(chromosomes)

    for ch in chromosomes:
        if args.not_compressed:
            totals_arr = os.path.join(args.array, f'{ch}.total')
            thresholds_arr = os.path.join(args.array, f'{ch}.thresholds')
        else:
            totals_arr = os.path.join(args.array, f'{ch}.total.gz')
            thresholds_arr = os.path.join(args.array, f'{ch}.thresholds.gz')
        if not os.path.exists(totals_arr):
            raise ValueError(error('Missing array file: {}'.format(totals_arr)))
        if not os.path.exists(thresholds_arr):
            raise ValueError(error('Missing array file: {}'.format(thresholds_arr)))

    log(msg=f'Identified {len(chromosomes)} chromosomes.\n', level='INFO')

    using_chr = [a.startswith('chr') for a in chromosomes]
    if any(using_chr):
        if not all(using_chr):
            raise ValueError(error("Some starts files use 'chr' notation while others do not."))
        use_chr = True
    else:
        use_chr = False

    ensure(args.processes > 0, 'The number of jobs must be positive.')
    ensure(
        args.msr > 0,
        'The minimum number of SNP-covering reads must be positive.',
    )
    ensure(args.mtr > 0, 'The minimum number of total reads must be positive.')

    ver = args.refversion
    ensure(
        ver in ('hg19', 'hg38'),
        'Invalid reference genome version. Supported versions are hg38 and hg19.',
    )

    outdir = os.sep.join(args.outfile.split(os.sep)[:-1])
    ensure(
        os.path.exists(outdir),
        f'Directory for output file does not exist: <{outdir}>',
    )

    return {
        'baffile': args.baffile,
        'outfile': args.outfile,
        'sample_names': names,
        'min_snp_reads': args.msr,
        'min_total_reads': args.mtr,
        'use_chr': use_chr,
        'processes': args.processes,
        'chromosomes': chromosomes,
        'array': args.array,
        'totalcounts': args.totalcounts,
        'phase': args.phase,
        'blocksize': args.max_blocksize,
        'max_snps_per_block': args.max_spb,
        'test_alpha': args.alpha,
        'multisample': not args.ss_em,
        'ref_version': ver,
    }


def parse_genotype_snps_arguments(args=None):
    description = 'Genotype and call SNPs in a matched-normal sample.'
    parser = argparse.ArgumentParser(prog='hatchet genotype-snps', description=description)
    parser.add_argument(
        '-N',
        '--normal',
        required=True,
        type=str,
        default = "None",
        help='BAM file corresponding to matched normal sample',
    )
    parser.add_argument(
        '-T',
        '--tumor',
        required=False,
        nargs = "+",
        type=str,
        default = [],
        help='BAM file corresponding to matched tumor sample(s)',
    )
    parser.add_argument(
        '-r',
        '--reference',
        required=True,
        type=str,
        help='Human reference genome of BAMs',
    )
    parser.add_argument(
        '-st',
        '--samtools',
        required=False,
        default=config.paths.samtools,
        type=str,
        help=(
            'Path to the directory to "samtools" executable, required in default mode (default: samtools is '
            'directly called as it is in user $PATH)'
        ),
    )
    parser.add_argument(
        '-bt',
        '--bcftools',
        required=False,
        default=config.paths.bcftools,
        type=str,
        help=(
            'Path to the directory of "bcftools" executable, required in default mode (default: bcftools is '
            'directly called as it is in user $PATH)'
        ),
    )
    parser.add_argument(
        '-R',
        '--snps',
        required=False,
        default=config.genotype_snps.snps,
        type=str,
        help=(
            'List of SNPs to consider in the normal sample (default: heterozygous SNPs are inferred from the normal '
            'sample)'
        ),
    )
    parser.add_argument(
        '-j',
        '--processes',
        required=False,
        default=config.genotype_snps.processes,
        type=int,
        help='Number of available parallel processes (default: 2)',
    )
    parser.add_argument(
        '-q',
        '--readquality',
        required=False,
        default=config.genotype_snps.readquality,
        type=int,
        help='Minimum mapping quality for an aligned read to be considered (default: 0)',
    )
    parser.add_argument(
        '-Q',
        '--basequality',
        required=False,
        default=config.genotype_snps.basequality,
        type=int,
        help='Minimum base quality for a base to be considered (default: 11)',
    )
    parser.add_argument(
        '-c',
        '--mincov',
        required=False,
        default=config.genotype_snps.mincov,
        type=int,
        help='Minimum coverage for SNPs to be considered (default: 0)',
    )
    parser.add_argument(
        '-C',
        '--maxcov',
        required=False,
        default=config.genotype_snps.maxcov,
        type=int,
        help=(
            'Maximum coverage for SNPs to be considered (default: 1000, suggested: twice the values of expected '
            'average coverage to avoid aligning artefacts)'
        ),
    )
    parser.add_argument(
        '-E',
        '--newbaq',
        required=False,
        action='store_true',
        default=config.genotype_snps.newbaq,
        help='Recompute alignment of reads on the fly during SNP calling (default: false)',
    )
    parser.add_argument(
        '-o',
        '--outputsnps',
        required=False,
        default=config.genotype_snps.outputsnps,
        type=str,
        help='Output folder for SNPs separated by chromosome (default: ./)',
    )
    parser.add_argument(
        '--chromosomes',
        required=False,
        nargs='*',
        help='One or more chromosomes to process (default: blank to process all chromosomes)',
    )
    parser.add_argument(
        '-v',
        '--verbose',
        action='store_true',
        default=config.genotype_snps.verbose,
        required=False,
        help='Use verbose log messages',
    )
    parser.add_argument('-V', '--version', action='version', version=f'%(prog)s {__version__}')
    args = parser.parse_args(args)

    # Parse BAM files, check their existence, and infer or parse the corresponding sample names
    nonormalFlag = False
    normalbaf = args.normal
    if normalbaf == "None":
        log(
            msg='Normal BAM file is not provided. Therefore the analysis will run without a matched normal.'
                'het SNP positions will be inferred from the first tumor sample\n',
            level='WARN',
        )
        normalbaf = args.tumor[0]
        nonormalFlag = True
    if not isfile(os.path.abspath(normalbaf)):
        raise ValueError(error(f'The specified normal BAM file {normalbaf} does not exist'))
    normal = (os.path.abspath(normalbaf), 'Normal')

    # In default mode, check the existence and compatibility of samtools and bcftools
    samtools = os.path.join(args.samtools, 'samtools')
    bcftools = os.path.join(args.bcftools, 'bcftools')
    ensure(
        which(samtools) is not None,
        'samtools has not been found or is not executable!',
    )
    ensure(
        which(bcftools) is not None,
        'bcftools has not been found or is not executable!',
    )

    # Check that SNP, reference, and region files exist when given in input
    ensure(
        isfile(args.reference),
        'The provided file for human reference genome does not exist!',
    )
    ensure(
        (args.outputsnps is None) or (isdir(args.outputsnps)),
        'The folder for output SNPs does not exist!',
    )

    if args.snps is not None and len(args.snps) < 2:
        args.snps = None
    if args.snps is not None and not (isfile(args.snps) or url_exists(args.snps)):
        error('The provided list of SNPs does not exist!', raise_exception=True)

    # Extract the names of the chromosomes and check their consistency across the given BAM files and the reference
    chromosomes = extractChromosomes(samtools, normal, [], args.reference)
    if args.chromosomes:
        chromosomes = [c for c in chromosomes if c in args.chromosomes]

    ensure(
        args.processes > 0,
        'The number of parallel processes must be greater than 0',
    )
    ensure(args.readquality >= 0, 'The read mapping quality must be positive')
    ensure(args.basequality >= 0, 'The base quality quality must be positive')
    ensure(args.mincov >= 0, 'The minimum-coverage value must be positive')
    ensure(args.maxcov >= 0, 'The maximum-coverage value must be positive')

    if args.verbose:
        log(
            msg='stderr of samtools and bcftools will be collected in the following file "samtools.log"\n',
            level='WARN',
        )
        with open('samtools.log', 'w') as f:
            f.write('')

    return {
        'normal': normal,
        'nonormal': nonormalFlag,
        'chromosomes': chromosomes,
        'samtools': samtools,
        'bcftools': bcftools,
        'snps': args.snps,
        'reference': args.reference,
        'j': args.processes,
        'q': args.readquality,
        'Q': args.basequality,
        'E': args.newbaq,
        'mincov': args.mincov,
        'maxcov': args.maxcov,
        'outputsnps': os.path.abspath(args.outputsnps),
        'verbose': args.verbose,
    }


def parse_download_panel_arguments(args=None):
    description = 'Download and prepare files for phasing germline SNPs using a reference panel'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        '-D',
        '--refpaneldir',
        required=False,
        type=str,
        default=config.download_panel.refpaneldir,
        help='Path to Reference Panel',
    )
    parser.add_argument(
        '-R',
        '--refpanel',
        default=config.download_panel.refpanel,
        required=False,
        type=str,
        help=(
            'Reference Panel; specify "1000GP_Phase3" to automatically download and use the panel form the 1000 '
            'genomes project'
        ),
    )
    args = parser.parse_args(args)

    ensure(
        args.refpaneldir is not None,
        'The command "download_panel" requires a path for the variable "refpaneldir".',
    )

    return {
        'refpanel': args.refpanel,
        'refpaneldir': os.path.abspath(args.refpaneldir),
    }


def parse_phase_snps_arguments(args=None):
    description = 'Phase germline SNPs using a reference panel'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        '-D',
        '--refpaneldir',
        required=False,
        type=str,
        default=config.download_panel.refpaneldir,
        help='Path to Reference Panel',
    )
    parser.add_argument(
        '-g',
        '--refgenome',
        required=True,
        type=str,
        help='Path to Reference genome used in BAM files',
    )
    parser.add_argument(
        '-V',
        '--refversion',
        required=True,
        type=str,
        help='Version of reference genome used in BAM files',
    )
    parser.add_argument(
        '-N',
        '--chrnotation',
        required=False,
        action='store_true',
        help='Use this flag to indicate that chromosomes are named with "chr" (default: no "chr")',
        default=False,
    )
    parser.add_argument(
        '-o',
        '--outdir',
        required=False,
        type=str,
        help='Output folder for phased VCFs',
    )
    parser.add_argument(
        '-L',
        '--snps',
        required=True,
        type=str,
        nargs='+',
        help='List of SNPs in the normal sample to phase',
    )
    parser.add_argument(
        '-j',
        '--processes',
        required=False,
        default=config.genotype_snps.processes,
        type=int,
        help='Number of available parallel processes (default: 2)',
    )
    parser.add_argument(
        '-si',
        '--shapeit',
        required=False,
        type=str,
        default=config.paths.shapeit,
        help='Path to shapeit executable (default: look in $PATH)',
    )
    parser.add_argument(
        '-pc',
        '--picard',
        required=False,
        type=str,
        default=config.paths.picard,
        help='Path to picard executable or jar (default: look in $PATH)',
    )
    parser.add_argument(
        '-bt',
        '--bcftools',
        required=False,
        default=config.paths.bcftools,
        type=str,
        help='Path to the directory of "bcftools" executable (default: look in $PATH)',
    )
    parser.add_argument(
        '-bg',
        '--bgzip',
        required=False,
        default=config.paths.bgzip,
        type=str,
        help='Path to the directory of "bcftools" executable (default: look in $PATH)',
    )
    args = parser.parse_args(args)

    # add safety checks for custom ref panel
    # if args.refpanel != "1000GP_Phase3":

    bgzip = os.path.join(args.bgzip, 'bgzip')
    if which(bgzip) is None:
        raise ValueError(
            error(
                "The 'bgzip' executable was not found or is not executable. \
            Please install bgzip and/or supply the path to the executable."
            )
        )

    shapeit = os.path.join(args.shapeit, 'shapeit')
    if which(shapeit) is None:
        raise ValueError(
            error(
                "The 'shapeit' executable was not found or is not executable. \
            Please install shapeit (e.g., conda install -c dranew shapeit) and/or supply the path to the executable."
            )
        )

    # We may have access to a picard wrapper executable or picard.jar. Take care of both scenarios.
    picard_dir = args.picard
    picard_java_flags = config.phase_snps.picard_java_flags
    picard_bin_path = os.path.join(picard_dir, 'picard')
    picard_jar_path = os.path.join(picard_dir, 'picard.jar')
    if which(picard_bin_path) is not None:
        picard = f'{picard_bin_path} {picard_java_flags}'
    elif os.path.exists(picard_jar_path):
        picard = f'java {picard_java_flags} -jar {picard_jar_path}'
    else:
        error(
            (
                'Neither picard nor picard.jar were found in the given location. Please install picard'
                '(e.g., conda install -c bioconda picard) and/or supply the path to the executable'
                ' or JAR file.'
            ),
            raise_exception=True,
        )

    bcftools = os.path.join(args.bcftools, 'bcftools')
    if which(bcftools) is None:
        raise ValueError(error('bcftools has not been found or is not executable!'))

    if args.refversion not in ('hg19', 'hg38'):
        error(
            'The reference genome version of your samples is not "hg19" or "hg38".',
            raise_exception=True,
        )

    rpd = args.refpaneldir
    if not os.path.exists(os.path.join(rpd, '1000GP_Phase3', '1000GP_Phase3.sample')):
        error(
            'Please download the 1000GP reference panel before proceeding',
            raise_exception=True,
        )

    # Check that SNP files exist when given in input
    snplists = {}
    for f in args.snps:
        if not isfile(f):
            raise ValueError(error('The specified SNP file {} does not exist!'.format(f)))
        # use keys that correspond to chromosomes names used (below)
        snplists = {os.path.basename(f).split('.')[0].replace('chr', ''): f for f in args.snps}

    # Get chromosome names from vcf file names, since they're named according to chromosome in SNPCaller
    # rename chromosomes if they have chr prefix; used to locate ref panel files!
    chromosomes = [os.path.basename(snplists[i]).replace('.vcf.gz', '').replace('chr', '') for i in snplists]

    return {
        'refpaneldir': args.refpaneldir,
        'j': args.processes,
        'chromosomes': chromosomes,
        'snps': snplists,
        'refvers': args.refversion,
        'chrnot': args.chrnotation,
        'refgenome': args.refgenome,
        'outdir': os.path.abspath(args.outdir),
        'shapeit': shapeit,
        'picard': picard,
        'bcftools': bcftools,
        'bgzip': bgzip,
    }


def parse_count_alleles_arguments(args=None):
    description = (
        'Count the A/B alleles from a matched-normal BAM file and multiple tumor BAM files in specified '
        'SNP positions or estimated heterozygous SNPs in the normal genome. This tool can be applied both '
        'to whole-genome sequencing (WGS) data or whole-exome sequencing (WES) data, but coding regions '
        'must be specified as a BED file in the case of WES.'
    )
    parser = argparse.ArgumentParser(prog='hatchet count-alleles', description=description)
    parser.add_argument(
        '-N',
        '--normal',
        required=True,
        type=str,
        help='BAM file corresponding to matched normal sample',
    )
    parser.add_argument(
        '-T',
        '--tumors',
        required=True,
        type=str,
        nargs='+',
        help='BAM files corresponding to samples from the same tumor',
    )
    parser.add_argument(
        '-r',
        '--reference',
        required=True,
        type=str,
        help='Human reference genome of BAMs',
    )
    parser.add_argument(
        '-L',
        '--snps',
        required=True,
        type=str,
        nargs='+',
        help='List of SNPs to consider in the normal sample',
    )
    parser.add_argument(
        '-S',
        '--samples',
        required=False,
        default=config.count_alleles.samples,
        type=str,
        nargs='+',
        help='Sample names for each BAM (given in the same order where the normal name is first)',
    )
    parser.add_argument(
        '-st',
        '--samtools',
        required=False,
        default=config.paths.samtools,
        type=str,
        help=(
            'Path to the directory to "samtools" executable, required in default mode (default: samtools is '
            'directly called as it is in user $PATH)'
        ),
    )
    parser.add_argument(
        '-bt',
        '--bcftools',
        required=False,
        default=config.paths.bcftools,
        type=str,
        help=(
            'Path to the directory of "bcftools" executable, required in default mode (default: bcftools is '
            'directly called as it is in user $PATH)'
        ),
    )
    parser.add_argument(
        '-j',
        '--processes',
        required=False,
        default=config.count_alleles.processes,
        type=int,
        help='Number of available parallel processes (default: 2)',
    )
    parser.add_argument(
        '-q',
        '--readquality',
        required=False,
        default=config.count_alleles.readquality,
        type=int,
        help='Minimum mapping quality for an aligned read to be considered (default: 0)',
    )
    parser.add_argument(
        '-Q',
        '--basequality',
        required=False,
        default=config.count_alleles.basequality,
        type=int,
        help='Minimum base quality for a base to be considered (default: 11)',
    )
    parser.add_argument(
        '-U',
        '--snpquality',
        required=False,
        default=config.count_alleles.snpquality,
        type=int,
        help='Minimum SNP-variant quality, QUAL, for a variant to be considered (default: 11)',
    )
    parser.add_argument(
        '-g',
        '--gamma',
        required=False,
        default=config.count_alleles.gamma,
        type=float,
        help='Level of confidence to determine heterozigosity of SNPs (default: 0.05)',
    )
    parser.add_argument(
        '-b',
        '--maxshift',
        required=False,
        default=config.count_alleles.maxshift,
        type=float,
        help=(
            'Maximum allowed absolute difference of BAF from 0.5 for selected heterozygous SNPs in the normal '
            'sample (default: 0.5)'
        ),
    )
    parser.add_argument(
        '-c',
        '--mincov',
        required=False,
        default=config.count_alleles.mincov,
        type=int,
        help='Minimum coverage for SNPs to be considered (default: 0)',
    )
    parser.add_argument(
        '-C',
        '--maxcov',
        required=False,
        default=config.count_alleles.maxcov,
        type=int,
        help=(
            'Maximum coverage for SNPs to be considered (default: 1000, suggested: twice the values of expected '
            'average coverage to avoid aligning artefacts)'
        ),
    )
    parser.add_argument(
        '-E',
        '--newbaq',
        required=False,
        action='store_true',
        default=config.count_alleles.newbaq,
        help='Recompute alignment of reads on the fly during SNP calling (default: false)',
    )
    parser.add_argument(
        '-O',
        '--outputnormal',
        required=False,
        default=config.count_alleles.outputnormal,
        type=str,
        help='Filename of output for allele counts in the normal sample (default: standard output)',
    )
    parser.add_argument(
        '-o',
        '--outputtumors',
        required=False,
        default=config.count_alleles.outputtumors,
        type=str,
        help='Output filename for allele counts in tumor samples (default: standard output)',
    )
    parser.add_argument(
        '-l',
        '--outputsnps',
        required=False,
        default=config.count_alleles.outputsnps,
        type=str,
        help='Output directory for lists of selected SNPs (default: ./)',
    )
    parser.add_argument(
        '--chromosomes',
        required=False,
        nargs='*',
        help='One or more chromosomes to process (default: blank to process all chromosomes)',
    )
    parser.add_argument(
        '-v',
        '--verbose',
        action='store_true',
        default=config.count_alleles.verbose,
        required=False,
        help='Use verbose log messages',
    )
    parser.add_argument('-V', '--version', action='version', version=f'%(prog)s {__version__}')
    args = parser.parse_args(args)

    # Parse BAM files, check their existence, and infer or parse the corresponding sample names
    normalbaf = args.normal
    if not isfile(normalbaf):
        raise ValueError(error('The specified normal BAM file does not exist'))
    tumors = args.tumors
    for tumor in tumors:
        if not isfile(tumor):
            raise ValueError(error(f'The specified tumor BAM file does not exist: {tumor}'))
    names = args.samples
    ensure(
        (names is None) or ((len(tumors) + 1) == len(names)),
        (
            'A sample name must be provided for each corresponding BAM: both for each normal sample and each ',
            'tumor sample',
        ),
    )

    samples = set()
    if names is None:
        normal = normalbaf, os.path.splitext(os.path.basename(normalbaf))[0]
        for tumor in tumors:
            samples.add((tumor, os.path.splitext(os.path.basename(tumor))[0]))
    else:
        normal = normalbaf, names[0]
        for i in range(len(tumors)):
            samples.add((tumors[i], names[i + 1]))

    # In default mode, check the existence and compatibility of samtools and bcftools
    samtools = os.path.join(args.samtools, 'samtools')
    bcftools = os.path.join(args.bcftools, 'bcftools')
    if which(samtools) is None:
        raise ValueError(error('{}samtools has not been found or is not executable!{}'))
    elif which(bcftools) is None:
        raise ValueError(error('{}bcftools has not been found or is not executable!{}'))

    # Check that SNP, reference, and region files exist when given in input
    snplists = {}
    for f in args.snps:
        if not isfile(f):
            raise ValueError(error('The specified SNP file {} does not exist!'.format(f)))
        snplists = {os.path.basename(f).split('.')[0]: f for f in args.snps}
    if not isfile(args.reference):
        raise ValueError(error('The provided file for human reference genome does not exist!'))

    # Extract the names of the chromosomes and check their consistency across the given BAM files and the reference
    chromosomes = extractChromosomes(samtools, normal, samples, args.reference)
    if args.chromosomes:
        chromosomes = [c for c in chromosomes if c in args.chromosomes]
    snplists = {c: snplists.get(c, []) for c in chromosomes}

    ensure(
        args.processes > 0,
        'The number of parallel processes must be greater than 0',
    )
    ensure(args.readquality >= 0, 'The read mapping quality must be positive')
    ensure(args.basequality >= 0, 'The base quality quality must be positive')
    ensure(0 <= args.gamma <= 1, 'Gamma must be a floating value between 0 and 1')
    ensure(
        0 <= args.maxshift <= 0.5,
        'Max BAF shift must be a floating value between 0 and 0.5',
    )
    ensure(args.mincov >= 0, 'The minimum-coverage value must be positive')
    ensure(args.maxcov >= 0, 'The maximum-coverage value must be positive')
    ensure(args.snpquality >= 0, 'The QUAL value must be positive')

    if args.verbose:
        log(
            msg='stderr of samtools and bcftools will be collected in the following file "samtools.log"\n',
            level='WARN',
        )
        with open('samtools.log', 'w') as f:
            f.write('')

    return {
        'normal': normal,
        'samples': samples,
        'chromosomes': chromosomes,
        'samtools': samtools,
        'bcftools': bcftools,
        'snps': snplists,
        'reference': args.reference,
        'j': args.processes,
        'q': args.readquality,
        'Q': args.basequality,
        'qual': args.snpquality,
        'E': args.newbaq,
        'gamma': args.gamma,
        'maxshift': args.maxshift,
        'mincov': args.mincov,
        'maxcov': args.maxcov,
        'outputNormal': args.outputnormal,
        'outputTumors': args.outputtumors,
        'outputSnps': args.outputsnps,
        'verbose': args.verbose,
    }


def parse_count_reads_fw_arguments(args=None):
    description = (
        'Count the mapped sequencing reads in bins of fixed and given length, uniformly for a BAM file of '
        'a normal sample and one or more BAM files of tumor samples. This program supports both data from '
        'whole-genome sequencing (WGS) and whole-exome sequencing (WES), but the a BED file with targeted '
        'regions is required when considering WES.'
    )
    parser = argparse.ArgumentParser(prog='hatchet count-reads', description=description)
    parser.add_argument(
        '-N',
        '--normal',
        required=True,
        type=str,
        help='BAM file corresponding to matched normal sample',
    )
    parser.add_argument(
        '-T',
        '--tumors',
        required=True,
        type=str,
        nargs='+',
        help='BAM files corresponding to samples from the same tumor',
    )
    parser.add_argument(
        '-b',
        '--size',
        required=True,
        type=str,
        help='Size of the bins, specified as a full number or using the notations either "kb" or "Mb"',
    )
    parser.add_argument(
        '-S',
        '--samples',
        required=False,
        default=config.count_reads_fw.samples,
        type=str,
        nargs='+',
        help=(
            'Sample names for each BAM, given in the same order where the normal name is first (default: inferred '
            'from file names)'
        ),
    )
    parser.add_argument(
        '-st',
        '--samtools',
        required=False,
        default=config.paths.samtools,
        type=str,
        help=(
            'Path to the directory to "samtools" executable, required in default mode (default: samtools is '
            'directly called as it is in user $PATH)'
        ),
    )
    parser.add_argument(
        '-r',
        '--regions',
        required=False,
        default=config.count_reads_fw.regions,
        type=str,
        help=(
            'BED file containing the a list of genomic regions to consider in the format "CHR  START  END", '
            'REQUIRED for WES data (default: none, consider entire genome)'
        ),
    )
    parser.add_argument(
        '-g',
        '--reference',
        required=False,
        default=config.paths.reference,
        type=str,
        help=(
            'Reference genome, note that reference must be indexed and the dictionary must exist in the same '
            'directory with the same name and .dict extension'
        ),
    )
    parser.add_argument(
        '-j',
        '--processes',
        required=False,
        default=config.count_reads_fw.processes,
        type=int,
        help='Number of available parallel processes (default: 2)',
    )
    parser.add_argument(
        '-q',
        '--readquality',
        required=False,
        default=config.count_reads_fw.readquality,
        type=int,
        help=(
            'Minimum mapping quality for an aligned read to be considered ',
            f'(default: {config.count_reads_fw.readquality})',
        ),
    )
    parser.add_argument(
        '-O',
        '--outputnormal',
        required=False,
        default=config.count_reads_fw.outputnormal,
        type=str,
        help='Filename of output for allele counts in the normal sample (default: standard output)',
    )
    parser.add_argument(
        '-o',
        '--outputtumors',
        required=False,
        default=config.count_reads_fw.outputtumors,
        type=str,
        help='Output filename for allele counts in tumor samples (default: standard output)',
    )
    parser.add_argument(
        '-t',
        '--outputtotal',
        required=False,
        default=config.count_reads_fw.outputtotal,
        type=str,
        help='Output filename for total read counts in all tumor samples (default: "total_read.counts")',
    )
    parser.add_argument(
        '-v',
        '--verbose',
        action='store_true',
        default=config.count_reads_fw.verbose,
        required=False,
        help='Use verbose log messages',
    )
    parser.add_argument(
        '--chromosomes',
        required=False,
        nargs='*',
        help='One or more chromosomes to process (default: blank to process all chromosomes)',
    )
    parser.add_argument('-V', '--version', action='version', version=f'%(prog)s {__version__}')
    args = parser.parse_args(args)

    # Parse BAM files, check their existence, and infer or parse the corresponding sample names
    normalbaf = args.normal
    ensure(isfile(normalbaf), 'The specified normal BAM file does not exist')
    tumors = args.tumors
    for tumor in tumors:
        ensure(
            isfile(tumor),
            f'The specified tumor BAM file {tumor} does not exist',
        )

    names = args.samples
    ensure(
        names is None or (len(tumors) + 1) == len(names),
        'A sample name must be provided for each corresponding BAM: both for each normal sample and each tumor sample',
    )

    samples = set()
    if names is None:
        normal = (normalbaf, os.path.splitext(os.path.basename(normalbaf))[0])
        for tumor in tumors:
            samples.add((tumor, os.path.splitext(os.path.basename(tumor))[0]))
    else:
        normal = (normalbaf, names[0])
        for i in range(len(tumors)):
            samples.add((tumors[i], names[i + 1]))

    # Check the region file
    ensure(
        args.regions is None or isfile(args.regions),
        'The specified region file does not exist',
    )

    # In default mode, check the existence and compatibility of samtools and bcftools
    samtools = os.path.join(args.samtools, 'samtools')
    if which(samtools) is None:
        raise ValueError(error('samtools has not been found or is not executable!'))

    # Check and parse the given size
    size = 0
    try:
        if args.size[-2:] == 'kb':
            size = int(args.size[:-2]) * 1000
        elif args.size[-2:] == 'Mb':
            size = int(args.size[:-2]) * 1000000
        else:
            size = int(args.size)
    except (IndexError, ValueError):
        error(
            'Size must be a number, optionally ending with either "kb" or "Mb"!',
            raise_exception=True,
        )

    # Check that either region file or available reference name are available
    ensure(
        any((args.reference, args.regions)),
        (
            'Please either provide a BED file of regions or specify a name of an available references for inferring '
            'maximum-chromosome lengths',
        ),
    )
    ensure(
        args.reference is None or isfile(args.reference),
        'The specified reference genome does not exist!',
    )

    refdict = os.path.splitext(args.reference)[0] + '.dict'
    ensure(
        args.reference is None or isfile(args.reference),
        (
            'The dictionary of the reference genome has not been found! Reference genome must be indexed and its '
            'dictionary must exist in the same directory with same name but extension .dict'
        ),
    )

    # Extract the names of the chromosomes and check their consistency across the given BAM files and the reference
    chromosomes = extractChromosomes(samtools, normal, samples)
    if args.chromosomes:
        chromosomes = [c for c in chromosomes if c in args.chromosomes]

    ensure(
        args.processes > 0,
        'The number of parallel processes must be greater than 0',
    )
    ensure(args.readquality >= 0, 'The read mapping quality must be positive')

    return {
        'normal': normal,
        'samples': samples,
        'chromosomes': chromosomes,
        'samtools': samtools,
        'regions': args.regions,
        'size': size,
        'reference': args.reference,
        'refdict': refdict,
        'j': args.processes,
        'q': args.readquality,
        'outputNormal': args.outputnormal,
        'outputTumors': args.outputtumors,
        'outputTotal': args.outputtotal,
        'verbose': args.verbose,
    }


def parse_combine_counts_fw_args(args=None):
    parser = argparse.ArgumentParser(
        prog='hatchet combine-counts',
        description=(
            'Combine tumor bin counts, normal bin counts, and tumor allele counts to obtain the read-depth '
            'ratio and the mean B-allel frequency (BAF) of each bin. Optionally, the normal allele counts '
            'can be provided to add the BAF of each bin scaled by the normal BAF. The output is written on '
            'stdout.'
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        '-c',
        '--normalbins',
        required=True,
        type=str,
        help='Normal bin counts in the format "SAMPLE\tCHR\tSTART\tEND\tCOUNT"',
    )
    parser.add_argument(
        '-C',
        '--tumorbins',
        required=True,
        type=str,
        help='Tumor bin counts in the format "SAMPLE\tCHR\tSTART\tEND\tCOUNT"',
    )
    parser.add_argument(
        '-B',
        '--tumorbafs',
        required=True,
        type=str,
        help='Tumor allele counts in the format "SAMPLE\tCHR\tPOS\tREF-COUNT\tALT-COUNT"',
    )
    parser.add_argument(
        '-p',
        '--phase',
        required=False,
        default=config.combine_counts.phase,
        type=str,
        help='Phasing of heterozygous germline SNPs in the format "CHR\tPOS\t<string containing 0|1 or 1|0>"',
    )
    parser.add_argument(
        '-d',
        '--diploidbaf',
        type=float,
        required=False,
        default=0.1,
        help=(
            'Maximum diploid-BAF shift used to select the bins whose BAF should be normalized by the normal when '
            'normalbafs is given (default: 0.1)'
        ),
    )
    parser.add_argument(
        '-l',
        '--blocklength',
        required=False,
        default=config.combine_counts.blocksize,
        type=str,
        help=(
            'Size of the haplotype blocks, specified as a full number or using the notations either "kb" or "Mb" ',
            '(default: 50kb)',
        ),
    )
    parser.add_argument(
        '-t',
        '--totalcounts',
        required=False,
        default=None,
        type=str,
        help=(
            'Total read counts in the format "SAMPLE\tCOUNT" used to normalize by the different number of reads '
            'extracted from each sample (default: none)'
        ),
    )
    parser.add_argument(
        '-g',
        '--gamma',
        type=float,
        required=False,
        default=0.05,
        help=(
            'Confidence level used to determine if a bin is copy neutral with BAF of 0.5 in the BINOMIAL_TEST mode '
            '(default: 0.05)'
        ),
    )
    parser.add_argument(
        '-e',
        '--seed',
        type=int,
        required=False,
        default=0,
        help='Random seed used for the normal distributions used in the clouds (default: 0)',
    )
    parser.add_argument(
        '-v',
        '--verbose',
        action='store_true',
        required=False,
        help='Use verbose log messages',
    )
    parser.add_argument(
        '-r',
        '--disablebar',
        action='store_true',
        required=False,
        help='Disable progress bar',
    )
    parser.add_argument('-V', '--version', action='version', version=f'%(prog)s {__version__}')
    args = parser.parse_args(args)

    ensure(
        isfile(args.normalbins),
        'The specified file for normal bin counts does not exist!',
    )
    ensure(
        isfile(args.tumorbins),
        'The specified file for tumor bin counts does not exist!',
    )
    ensure(
        isfile(args.normalbins),
        'The specified file for normal bin counts does not exist!',
    )
    if args.phase == 'None':  # TODO - why do we need this?
        args.phase = None
    ensure(
        (args.phase is None) or (isfile(args.phase)),
        'The specified file for phase does not exist!',
    )
    ensure(
        0.0 <= args.diploidbaf <= 0.5,
        'The specified maximum for diploid-BAF shift must be a value in [0.0, 0.5]',
    )
    ensure(
        (args.totalcounts is None) or (isfile(args.totalcounts)),
        'The specified file for total read counts does not exist!',
    )
    ensure(
        0.0 <= args.gamma <= 0.1,
        'The specified gamma must be a value in [0.0, 0.1]',
    )
    ensure(args.seed >= 0, 'Seed parameter must be positive!')

    args.blocklength = str(args.blocklength)
    try:
        if args.blocklength[-2:] == 'kb':
            size = int(args.blocklength[:-2]) * 1000
        elif args.blocklength[-2:] == 'Mb':
            size = int(args.blocklength[:-2]) * 1000000
        else:
            size = int(args.blocklength)
    except (IndexError, ValueError):
        error(
            'Size must be a number, optionally ending with either "kb" or "Mb"!',
            raise_exception=True,
        )

    return {
        'normalbins': args.normalbins,
        'tumorbins': args.tumorbins,
        'tumorbafs': args.tumorbafs,
        'phase': args.phase,
        'block': size,
        'diploidbaf': args.diploidbaf,
        'totalcounts': args.totalcounts,
        'gamma': args.gamma,
        'seed': args.seed,
        'verbose': args.verbose,
        'disable': args.disablebar,
    }


def parse_cluster_bins_gmm_args(args=None):
    parser = argparse.ArgumentParser(
        prog='hatchet cluster-bins',
        description=(
            'Combine tumor bin counts, normal bin counts, and tumor allele counts to obtain the read-depth '
            'ratio and the mean B-allel frequency (BAF) of each bin. Optionally, the normal allele counts '
            'can be provided to add the BAF of each bin scaled by the normal BAF. The output is written on '
            'stdout.'
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        'BBFILE',
        help=(
            'A BB file containing a line for each bin in each sample and the corresponding values of read-depth '
            'ratio and B-allele frequency (BAF)'
        ),
    )
    parser.add_argument(
        '-o',
        '--outsegments',
        required=False,
        default=config.cluster_bins_gmm.outsegments,
        type=str,
        help='Output filename for the segments computed by clustering bins (default: stdout)',
    )
    parser.add_argument(
        '-O',
        '--outbins',
        required=False,
        default=config.cluster_bins_gmm.outbins,
        type=str,
        help='Output filename for a BB file adding the clusters (default: stdout)',
    )
    parser.add_argument(
        '-d',
        '--diploidbaf',
        type=float,
        required=False,
        default=config.cluster_bins_gmm.diploidbaf,
        help=(
            'Maximum diploid-BAF shift used to determine the largest copy-neutral cluster and to rescale all the '
            'cluster inside this threshold accordingly (default: None, scaling is not performed)'
        ),
    )
    parser.add_argument(
        '-tR',
        '--tolerancerdr',
        type=float,
        required=False,
        default=config.cluster_bins_gmm.tolerancerdr,
        help=(
            'Refine the clustering merging the clusters with this maximum difference in RDR values (default: None, '
            'gurobipy required)'
        ),
    )
    parser.add_argument(
        '-tB',
        '--tolerancebaf',
        type=float,
        required=False,
        default=config.cluster_bins_gmm.tolerancebaf,
        help=(
            'Refine the clustering merging the clusters with this maximum difference in BAF values (default: None, ',
            'gurobipy required)',
        ),
    )
    parser.add_argument(
        '-u',
        '--bootclustering',
        type=int,
        required=False,
        default=config.cluster_bins_gmm.bootclustering,
        help=(
            'Number of points to add for bootstraping each bin to improve the clustering. Each point is generated '
            'by drawing its values from a normal distribution centered on the values of the bin. This can help the '
            'clustering when the input number of bins is low (default: 0)'
        ),
    )
    parser.add_argument(
        '-dR',
        '--ratiodeviation',
        type=float,
        required=False,
        default=config.cluster_bins_gmm.ratiodeviation,
        help='Standard deviation of the read ratios used to generate the points in the clouds (default: 0.02)',
    )
    parser.add_argument(
        '-dB',
        '--bafdeviation',
        type=float,
        required=False,
        default=config.cluster_bins_gmm.bafdeviation,
        help='Standard deviation of the BAFs used to generate the points in the clouds (default: 0.02)',
    )
    parser.add_argument(
        '-e',
        '--seed',
        type=int,
        required=False,
        default=config.cluster_bins_gmm.seed,
        help='Random seed used for clustering AND the normal distributions used in the clouds (default: 0)',
    )
    parser.add_argument(
        '-K',
        '--initclusters',
        type=int,
        required=False,
        default=config.cluster_bins_gmm.initclusters,
        help='The maximum number of clusters to infer (default: 50)',
    )
    parser.add_argument(
        '-c',
        '--concentration',
        type=float,
        required=False,
        default=config.cluster_bins_gmm.concentration,
        help=(
            'Tuning parameter for clustering (concentration parameter for Dirichlet process prior). Higher favors '
            'more clusters, lower favors fewer clusters (default 0.02 = 1/K).'
        ),
    )
    parser.add_argument(
        '-R',
        '--restarts',
        type=int,
        required=False,
        default=config.cluster_bins_gmm.restarts,
        help='Number of restarts performed by the clustering to choose the best (default: 10)',
    )
    parser.add_argument(
        '-v',
        '--verbose',
        action='store_true',
        default=config.cluster_bins_gmm.verbose,
        required=False,
        help='Use verbose log messages',
    )
    parser.add_argument(
        '--disablebar',
        action='store_true',
        default=config.cluster_bins_gmm.disablebar,
        required=False,
        help='Disable progress bar',
    )
    parser.add_argument('-V', '--version', action='version', version=f'%(prog)s {__version__}')
    args = parser.parse_args(args)

    ensure(isfile(args.BBFILE), 'The specified BB file does not exist!')
    ensure(
        (args.diploidbaf is None) or (0.0 <= args.diploidbaf <= 0.5),
        'The specified maximum for diploid-BAF shift must be a value in [0.0, 0.5]',
    )
    ensure(args.tolerancerdr >= 0, 'Tolerance-RDR parameter must be positive!')
    ensure(args.tolerancebaf >= 0, 'Tolerance-BAF parameter must be positive!')
    ensure(args.bootclustering >= 0, 'Bootclustering parameter must be positive!')
    ensure(args.ratiodeviation >= 0, 'Ratio-deviation parameter must be positive!')
    ensure(args.bafdeviation >= 0, 'BAF-deviation parameter must be positive!')
    ensure(args.seed >= 0, 'Seed parameter must be positive!')
    ensure(args.initclusters >= 0, 'Init-cluster parameter must be positive!')
    ensure(args.concentration >= 0, 'Concentration parameter must be positive!')
    ensure(args.restarts >= 1, 'Number of restarts must be positive!')

    return {
        'bbfile': args.BBFILE,
        'cloud': args.bootclustering,
        'diploidbaf': args.diploidbaf,
        'rdtol': args.tolerancerdr,
        'baftol': args.tolerancebaf,
        'ratiodeviation': args.ratiodeviation,
        'bafdeviation': args.bafdeviation,
        'seed': args.seed,
        'initclusters': args.initclusters,
        'concentration': args.concentration,
        'restarts': args.restarts,
        'verbose': args.verbose,
        'disable': args.disablebar,
        'outsegments': args.outsegments,
        'outbins': args.outbins,
    }


def parse_plot_bins_args(args=None):
    parser = argparse.ArgumentParser(
        prog='hatchet plot-bins',
        description=(
            'Generate plots for read-depth ratio (RD), B-allele frequency (BAF), and clusters for genomic ',
            'bins in multiple samples using .bb, .cbb, .seg files.',
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument('INPUT', help='Input BBC file with RDR and BAF')

    b, e = bcolors.BOLD, bcolors.ENDC
    parser.add_argument(
        '-c',
        '--command',
        required=False,
        type=str,
        default=config.plot_bins.command,
        help=(
            f'The command determining the plots to generate (default: all)\n\n\t{b}RD{e}: Plot the read-depth ratio ',
            f'(RD) values of the genomes for each sample.\n\n\t{b}CRD{e}: Plot the read-depth ratio (CRD) values of ',
            f'the genomes for each sample colored by corresponding cluster.\n\n\t{b}BAF{e}: Plot the B-allele ',
            f'frequency (BAF) values of the genomes for each sample.\n\n\t{b}CBAF{e}: Plot BAF values for each ',
            f'sample colored by corresponding cluster.\n\n\t{b}BB{e}: Plot jointly the values of read-depth ratio ',
            f'(RD) and B-allele frequency (BAF) for each bin in all samples and their density.\n\n\t{b}CBB{e}: Plot ',
            'jointly the values of read-depth ratio (RD) and B-allele frequency (BAF) for each bin in all samples ',
            f'by coloring the bins depending on their cluster.\n\n\t{b}CLUSTER{e}: Plot jointly the values of '
            f'read-depth ratio (RD) and B-allele frequency (BAF) for each cluster in all samples where the size of ',
            'the markers is proportional to the number of bins.',
        ),
    )
    parser.add_argument(
        '-s',
        '--segfile',
        required=False,
        type=str,
        default=config.plot_bins.segfile,
        help='When the corresponding seg file is provided the clusters are also plotted (default: none)',
    )
    parser.add_argument(
        '-m',
        '--colormap',
        required=False,
        type=str,
        default=config.plot_bins.colormap,
        help=(
            'Colormap to use for the colors in the plots, the available colormaps are the following '
            '{Set1, Set2, Paired, Dark2, tab10, tab20}'
        ),
    )
    parser.add_argument(
        '-tC',
        '--chrthreshold',
        required=False,
        type=int,
        default=config.plot_bins.chrthreshold,
        help='Only covering at least this number of chromosomes are considered (default: None)',
    )
    parser.add_argument(
        '-tS',
        '--sizethreshold',
        required=False,
        type=float,
        default=config.plot_bins.sizethreshold,
        help='Only covering at least this genome proportion (default: None)',
    )
    parser.add_argument(
        '--resolution',
        required=False,
        default=config.plot_bins.resolution,
        type=int,
        help='Resolution of bins (default: bins are not merged)',
    )
    parser.add_argument(
        '--xmin',
        required=False,
        default=config.plot_bins.xmin,
        type=float,
        help='Minimum value on x-axis for supported plots (default: inferred from data)',
    )
    parser.add_argument(
        '--xmax',
        required=False,
        default=config.plot_bins.xmax,
        type=float,
        help='Maximum value on x-axis for supported plots (default: inferred from data)',
    )
    parser.add_argument(
        '--ymin',
        required=False,
        default=config.plot_bins.ymin,
        type=float,
        help='Minimum value on y-axis for supported plots (default: inferred from data)',
    )
    parser.add_argument(
        '--ymax',
        required=False,
        default=config.plot_bins.ymax,
        type=float,
        help='Maximum value on y-axis for supported plots (default: inferred from data)',
    )
    parser.add_argument(
        '--figsize',
        required=False,
        default=config.plot_bins.figsize,
        type=str,
        help='Size of the plotted figures in the form "(X-SIZE, Y-SIZE)"',
    )
    parser.add_argument(
        '--markersize',
        required=False,
        default=config.plot_bins.markersize,
        type=int,
        help='Size of the markers (default: values inferred for each plot)',
    )
    parser.add_argument(
        '--colwrap',
        required=False,
        default=config.plot_bins.colwrap,
        type=int,
        help='Wrapping the plots in this number of columnes (default: 2)',
    )
    parser.add_argument(
        '--fontscale',
        required=False,
        default=config.plot_bins.fontscale,
        type=float,
        help='Font scale (default: 1)',
    )
    parser.add_argument(
        '-x',
        '--rundir',
        required=False,
        default=config.plot_bins.rundir,
        type=str,
        help='Running dirrectory where output the results (default: current directory)',
    )
    parser.add_argument(
        '--pdf',
        action='store_true',
        default=config.plot_bins.pdf,
        required=False,
        help='Output the bb_clustered figure in PDF format (default: PNG)',
    )
    parser.add_argument(
        '--dpi',
        required=False,
        default=config.plot_bins.dpi,
        type=int,
        help='DPI of PNG images (default: 900)',
    )
    parser.add_argument('-V', '--version', action='version', version=f'%(prog)s {__version__}')
    args = parser.parse_args(args)

    ensure(isfile(args.INPUT), 'The specified BB file does not exist!')
    ensure(
        args.command is None
        or args.command
        in (
            'RD',
            'CRD',
            'BAF',
            'CBAF',
            'BB',
            'CBB',
            'CLUSTER',
        ),
        'Unrecognized COMMAND!',
    )
    ensure(
        args.segfile is None or isfile(args.segfile),
        'Specified seg-file does not exist!',
    )

    ensure(
        args.colormap
        in (
            'Set1',
            'Set2',
            'Set3',
            'Paired',
            'Accent',
            'Dark2',
            'tab10',
            'tab20',
            'husl',
            'hls',
            'muted',
            'colorblind',
            'Pastel1',
            'Pastel2',
        ),
        'Unrecognized colormap!',
    )
    ensure(
        args.resolution is None or args.resolution >= 1,
        'Resolution must be greater than 1!',
    )
    ensure(
        args.chrthreshold is None or (0 <= args.chrthreshold <= 22),
        'The chromosome threshold must be a integer in [0, 22]',
    )
    ensure(args.colwrap >= 1, 'Colwrap must be grater than 1!')
    ensure(args.fontscale >= 0, 'Font scale must be positive!')
    ensure(
        args.sizethreshold is None or (0 <= args.sizethreshold <= 1.0),
        'The size threshold must be a real value in [0, 1]!',
    )
    if args.figsize is not None:
        try:
            parsed = args.figsize.strip().split(',')
            figsize = (float(parsed[0]), float(parsed[1]))
        except IndexError:
            error('Wrong format of figsize!', raise_exception=True)
    else:
        figsize = None
    ensure(
        isdir(args.rundir),
        'Running directory either does not exists or is not a directory!',
    )

    return {
        'input': args.INPUT,
        'command': args.command,
        'segfile': args.segfile,
        'ct': args.chrthreshold,
        'st': args.sizethreshold,
        'cmap': args.colormap,
        'resolution': args.resolution,
        'xmin': args.xmin,
        'xmax': args.xmax,
        'ymin': args.ymin,
        'ymax': args.ymax,
        'x': args.rundir,
        'figsize': figsize,
        'markersize': args.markersize,
        'colwrap': args.colwrap,
        'fontscale': args.fontscale,
        'pdf': args.pdf,
        'dpi': args.dpi,
    }


def extractChromosomes(samtools, normal, tumors, reference=None):
    """
    Parameters:
        samtools: path to samtools executable
        normal: tuple of (path to normal BAM file, string name)
        tumor: list of tuples (path to BAM file, string name)
        reference: path to FASTA file
    """

    # Read the names of sequences in normal BAM file
    normal_sq = getSQNames(samtools, normal[0])

    # Extract only the names of chromosomes in standard formats
    chrm = set()
    no_chrm = set()
    for i in range(1, 23):
        if str(i) in normal_sq:
            no_chrm.add(str(i))
        elif 'chr' + str(i) in normal_sq:
            chrm.add('chr' + str(i))
        else:
            sys.stderr.write(
                (
                    f'WARNING: a chromosome named either {i} or a variant of CHR{i} cannot be found in the normal ',
                    'BAM file\n',
                )
            )

    for c in ['X', 'Y']:
        if c in normal_sq:
            no_chrm.add(c)
        elif 'chr' + c in normal_sq:
            chrm.add('chr' + c)
        else:
            sys.stderr.write(
                (
                    f'WARNING: a chromosome named either {c} or a variant of CHR{c} cannot be found in the normal ',
                    'BAM file\n',
                )
            )

    if len(chrm) == 0 and len(no_chrm) == 0:
        raise ValueError('No chromosomes found in the normal BAM')
    chromosomes = chrm if len(chrm) > len(no_chrm) else no_chrm

    # Check that chromosomes with the same names are present in each tumor BAM contain
    for tumor in tumors:
        tumor_sq = getSQNames(samtools, tumor[0])
        if not chromosomes <= tumor_sq:
            sys.stderr.write(
                'WARNING: chromosomes {} are not present in the tumor sample {}\n'.format(chromosomes - tumor_sq, tumor)
            )

    # Check consistency of chromosome names with the reference
    if reference is not None:
        stdout, stderr = subprocess.Popen(
            'grep -e "^>" {}'.format(reference),
            stdout=subprocess.PIPE,
            shell=True,
            universal_newlines=True,
        ).communicate()
        if stderr is not None:
            raise ValueError('Error in reading the reference: {}'.format(reference))
        else:
            ref = set(c[1:].strip().split()[0] for c in stdout.strip().split('\n'))
        ensure(
            chromosomes <= ref,
            (
                'The given reference cannot be used because the chromosome names are inconsistent!\nChromosomes found '
                f'in BAF files: {chromosomes}\nChromosomes with the same name found in reference genome: {ref}',
            ),
        )

    return sorted(list(chromosomes), key=numericOrder)


def getSQNames(samtools, bamfile):
    header, stderr = subprocess.Popen(
        [samtools, 'view', '-H', bamfile],
        stdout=subprocess.PIPE,
        shell=False,
        universal_newlines=True,
    ).communicate()
    if stderr is not None:
        raise ValueError('The header of the normal-sample BAM cannot be read with samtools!')
    names = set()
    for line in header.strip().split('\n'):
        line = line.split()
        if len(line) > 0 and line[0] == '@SQ':
            names.add(line[1].split(':')[1].strip())
    return names


def parseRegions(region_file, chromosomes):
    res = {}
    for chro in chromosomes:
        res[chro] = []
    nofound = set()
    with open(region_file, 'r') as f:
        for line in f:
            split = line.strip().split()
            chro = split[0]
            if chro in chromosomes:
                res[chro].append((int(split[1]), int(split[2])))
            else:
                nofound.add(chro)

    for c in nofound:
        log(
            msg=(
                f'The chromosome {c} present in the provided regions is non-autosome or is not present in the given '
                'BAM files\n'
            ),
            level='WARN',
        )

    for key in res:
        ensure(
            len(res[key]) != 0,
            f'The regions of chromosome {key} could not be determined.',
        )
        res[key].sort(key=lambda x: x[0])
        if not all(a[0] <= a[1] <= b[0] <= b[1] for a, b in zip(res[key], res[key][1:])):
            error(
                (
                    f'The regions provided for chromosome {key} are non-disjoint or a region start is greater than ',
                    'corresponding region end',
                ),
                raise_exception=True,
            )

    return res
