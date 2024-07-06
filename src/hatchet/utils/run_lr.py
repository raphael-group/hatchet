from io import StringIO
import os.path
import glob
import argparse

from hatchet import config
from hatchet.utils.count_reads import main as count_reads
from hatchet.utils.count_reads_fw import main as count_reads_fw
from hatchet.utils.lr_functions import genotype_snps, phase_snps, combine_counts
from hatchet.utils.count_alleles import main as count_alleles
from hatchet.utils.combine_counts_fw import main as combine_counts_fw
from hatchet.utils.cluster_bins_gmm import main as cluster_bins_gmm
from hatchet.utils.cluster_bins import main as cluster_bins
from hatchet.utils.plot_bins import main as plot_bins
from hatchet.utils.plot_bins_1d2d import main as plot_bins_1d2d
from hatchet.bin.HATCHet import main as hatchet_main
from hatchet.utils.plot_cn import main as plot_cn
from hatchet.utils.plot_cn_1d2d import main as plot_cn_1d2d
from hatchet.utils.Supporting import log, error


def main(args=None):

    parser = argparse.ArgumentParser(prog='hatchet run-ont', description='Run HATCHet pipeline on ONT long read data')
    parser.add_argument('inifile', help='.ini file for run configuration')
    args = parser.parse_args(args)

    config.read(args.inifile)

    output = config.run.output
    output = output.rstrip('/')
    os.makedirs(output, exist_ok=True)

    try:
        chromosomes = [c for c in config.run.chromosomes.split()]
    except (KeyError, AttributeError):  # if key is absent or is blank (None)
        chromosomes = []  # process all

    extra_args = []
    try:
        if config.run.processes is not None:
            extra_args = ['-j', str(config.run.processes)]
    except KeyError:
        pass

    # ----------------------------------------------------

    if config.run.genotype_snps:
        genotype_file = config.run.ont_genotype_file
        if not os.path.exists(genotype_file):
            raise RuntimeError(
                (
                    'Please specify a valid VCF file (genotype_file) which'
                    'contains genotyped & phased SNPs for the normal'
                )
            )
        os.makedirs(f'{output}/snps', exist_ok=True)
        genotype_snps(
            args=[
                '-N',
                config.run.normal,
                '-r',
                config.run.reference,
                '-R',
                genotype_file,
                '-o',
                f'{output}/snps/',
                '--chromosomes',
            ]
            + (chromosomes or [])  # important to keep this as a list here to allow proper argparse parsing
            + extra_args
        )

    # ----------------------------------------------------

    if config.run.phase_snps:
        genotype_file = config.run.ont_genotype_file
        if not os.path.exists(genotype_file):
            raise RuntimeError(
                (
                    'Please specify a valid VCF file (genotype_file) which'
                    'contains genotyped & phased SNPs for the normal'
                )
            )

        os.makedirs(f'{output}/phase', exist_ok=True)
        phase_snps(
            args=[
                '-N',
                config.run.normal,
                '-r',
                config.run.reference,
                '-R',
                genotype_file,
                '-o',
                f'{output}/phase/',
                '--chromosomes',
            ]
            + (chromosomes or [])  # important to keep this as a list here to allow proper argparse parsing
            + extra_args
        )

    # ----------------------------------------------------
    if config.run.count_alleles:
        os.makedirs(f'{output}/baf', exist_ok=True)
        count_alleles(
            args=['-N', config.run.normal, '-T']
            + config.run.bams.split()
            + ['-S']
            + ('normal ' + config.run.samples).split()
            + ['-r', config.run.reference, '-L']
            + glob.glob(f'{output}/snps/*.vcf.gz')
            + [
                '-O',
                f'{output}/baf/normal.1bed',
                '-o',
                f'{output}/baf/tumor.1bed',
                '-l',
                f'{output}',
                '--chromosomes',
            ]
            + (chromosomes or [])  # important to keep this as a list here to allow proper argparse parsing
            + extra_args
        )

    if config.run.fixed_width is None or not config.run.fixed_width:
        # ----------------------------------------------------
        # ----------------------------------------------------
        # Variable-width/adaptive binning

        if config.run.count_reads:

            os.makedirs(f'{output}/rdr', exist_ok=True)
            params = [
                '-N', config.run.normal, '-T',
                *config.run.bams.split(),
                '-S',
                *('normal ' + config.run.samples).split(),
                '-V', config.genotype_snps.reference_version,
                '-b', f'{output}/baf/tumor.1bed',
                '-O', f'{output}/rdr',
                '--segfile', config.count_reads.segfile,
                '--chromosomes',
                *(chromosomes or []),  # important to keep this as a list here to allow proper argparse parsing
                *extra_args
            ]
            count_reads(args=params)

        # ----------------------------------------------------

        if config.run.combine_counts:
            haplotype_file = config.run.ont_haplotype_file
            mosdepth_files = config.run.ont_mosdepth_files.split()
            if not os.path.exists(haplotype_file):
                raise RuntimeError(
                    (
                        'Please specify a valid GTF file (haplotype_file) which'
                        'contains haplotype blocks for the patient'
                    )
                )
            import sys

            _stdout = sys.stdout
            sys.stdout = StringIO()
            os.makedirs(f'{output}/bb', exist_ok=True)

            phasefile = f'{output}/phase/phased.vcf.gz'
            args = [
                '-A',
                f'{output}/rdr',
                '-b',
                f'{output}/baf/tumor.1bed',
                '-t',
                f'{output}/rdr/total.tsv',
                '-V',
                config.genotype_snps.reference_version,
                '-o',
                f'{output}/bb/bulk.bb',
                '-r',
                config.run.reference,
            ] + extra_args

            if os.path.exists(phasefile):
                log(
                    msg='Found phasing file, including phasing in binning process.\n',
                    level='INFO',
                )
                args = ['-p', f'{output}/phase/phased.vcf.gz'] + args
            else:
                log(
                    msg=f'NO PHASING FILE FOUND at {phasefile}. Not including phasing in binning process.\n',
                    level='WARNING',
                )
                raise RuntimeError(
                    (
                        'Please run phase_snps step.'
                    )
                )

            combine_counts(args, haplotype_file, mosdepth_files, config.run.bams)

    else:
        # ----------------------------------------------------
        # ----------------------------------------------------
        # Old fixed-width binning

        # throw an exception with the error message saying that fixed_width is not supported with long reads
        raise RuntimeError(
            (
                'Fixed-width binning is not supported with long reads. Please use the adaptive binning method.'
            )
        )

    if config.run.cluster_bins:
        os.makedirs(f'{output}/bbc', exist_ok=True)

        if config.run.loc_clust:
            cluster_bins(
                args=[
                    f'{output}/bb/bulk.bb',
                    '-o',
                    f'{output}/bbc/bulk.seg',
                    '-O',
                    f'{output}/bbc/bulk.bbc',
                ]
            )
        else:
            cluster_bins_gmm(
                args=[
                    f'{output}/bb/bulk.bb',
                    '-o',
                    f'{output}/bbc/bulk.seg',
                    '-O',
                    f'{output}/bbc/bulk.bbc',
                ]
            )

    # ----------------------------------------------------

    if config.run.plot_bins:
        os.makedirs(f'{output}/plots', exist_ok=True)
        plot_bins(
            args=[
                f'{output}/bbc/bulk.bbc',
                '--rundir',
                f'{output}/plots',
                '--ymin',
                '0',
                '--ymax',
                '3',
            ]
        )

        os.makedirs(f'{output}/plots/1d2d', exist_ok=True)
        plot_bins_1d2d(
            args=[
                '-b',
                f'{output}/bbc/bulk.bbc',
                '-s',
                f'{output}/bbc/bulk.seg',
                '--outdir',
                f'{output}/plots/1d2d',
                '--centers',
                '--centromeres',
            ]
        )

    # ----------------------------------------------------

    if config.run.compute_cn:
        os.makedirs(f'{output}/results', exist_ok=True)
        hatchet_main(args=['-x', f'{output}/results', '-i', f'{output}/bbc/bulk'] + extra_args)

    # ----------------------------------------------------

    if config.run.plot_cn:
        os.makedirs(f'{output}/summary', exist_ok=True)
        plot_cn(
            args=[
                f'{output}/results/best.bbc.ucn',
                '--rundir',
                f'{output}/summary',
            ]
        )

        os.makedirs(f'{output}/summary/1d2d', exist_ok=True)
        plot_cn_1d2d(
            args=[f'{output}/results/best.bbc.ucn', '--outdir', f'{output}/summary/1d2d', '--bysample', '--centromeres']
        )


if __name__ == '__main__':
    main()
