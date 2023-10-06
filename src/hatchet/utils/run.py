import sys
from io import StringIO
import os.path
import glob
import argparse

from hatchet import config
from hatchet.utils.count_reads import main as count_reads
from hatchet.utils.count_reads_fw import main as count_reads_fw
from hatchet.utils.genotype_snps import main as genotype_snps
from hatchet.utils.count_alleles import main as count_alleles
from hatchet.utils.combine_counts import main as combine_counts
from hatchet.utils.combine_counts_fw import main as combine_counts_fw
from hatchet.utils.cluster_bins_gmm import main as cluster_bins_gmm
from hatchet.utils.cluster_bins import main as cluster_bins
from hatchet.utils.plot_bins import main as plot_bins
from hatchet.utils.plot_bins_1d2d import main as plot_bins_1d2d
from hatchet.bin.HATCHet import main as hatchet_main
from hatchet.utils.plot_cn import main as plot_cn
from hatchet.utils.plot_cn_1d2d import main as plot_cn_1d2d
from hatchet.utils.download_panel import main as download_panel
from hatchet.utils.phase_snps import main as phase_snps
from hatchet.utils.Supporting import log, error


def main(args=None):

    parser = argparse.ArgumentParser(prog='hatchet run', description='Run HATCHet pipeline')
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

    if config.run.download_panel:
        if not config.download_panel.refpaneldir:
            raise ValueError(
                error(
                    (
                        'The step "download_panel" requires that the variable "refpaneldir" indicates the directory in '
                        'which to store the reference panel.'
                    )
                )
            )

        download_panel(
            args=[
                '-D',
                config.download_panel.refpaneldir,
                '-R',
                config.download_panel.refpanel,
            ]
        )

    # ----------------------------------------------------

    if config.run.genotype_snps:
        snps = ''
        if config.genotype_snps.snps:
            snps = config.genotype_snps.snps
        elif config.genotype_snps.reference_version:
            snps_mapping = {
                (
                    'hg19',
                    True,
                ): 'https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/GATK/00-All.vcf.gz',
                (
                    'hg19',
                    False,
                ): 'https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/00-All.vcf.gz',
                (
                    'hg38',
                    True,
                ): 'https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/GATK/00-All.vcf.gz',
                (
                    'hg38',
                    False,
                ): 'https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/00-All.vcf.gz',
            }

            if (
                config.genotype_snps.reference_version,
                config.genotype_snps.chr_notation,
            ) not in snps_mapping:
                raise RuntimeError(
                    (
                        'Please specify valid values of reference_version and chr_notation. '
                        f'Valid pairs include: {snps_mapping.keys()}'
                    )
                )
            else:
                snps = snps_mapping[
                    config.genotype_snps.reference_version,
                    config.genotype_snps.chr_notation,
                ]

        os.makedirs(f'{output}/snps', exist_ok=True)
        genotype_snps(
            args=[
                '-N',
                config.run.normal,
                '-T',
                config.run.bams,
                '-r',
                config.run.reference,
                '-R',
                snps,
                '-o',
                f'{output}/snps/',
                '--chromosomes',
            ]
            + (chromosomes or [])  # important to keep this as a list here to allow proper argparse parsing
            + extra_args
        )

    # ----------------------------------------------------

    if config.run.phase_snps:
        if len(glob.glob(f'{output}/snps/*.vcf.gz')) == 0:
            raise ValueError(
                error(
                    (
                        'No SNP files were found for phasing. Are there any .vcf.gz files in the snps subdirectory of '
                        'the output folder? Try running genotype_snps.'
                    )
                )
            )

        if not config.download_panel.refpaneldir:
            raise ValueError(
                error(
                    (
                        'The step "phase_snps" requires that the config variable "download_panel.refpaneldir" '
                        'indicates the directory where the reference panel is located.'
                    )
                )
            )

        os.makedirs(f'{output}/phase', exist_ok=True)
        phase_snps(
            args=[
                '-D',
                config.download_panel.refpaneldir,
                '-g',
                config.run.reference,
                '-V',
                config.genotype_snps.reference_version,
                '-o',
                f'{output}/phase/',
                '-L',
            ]
            + glob.glob(f'{output}/snps/*.vcf.gz')
            + (['-N'] if config.genotype_snps.chr_notation else [])
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
            count_reads(
                args=['-N', config.run.normal, '-T']
                + config.run.bams.split()
                + ['-S']
                + ('normal ' + config.run.samples).split()
                + [
                    '-V',
                    config.genotype_snps.reference_version,
                    '-b',
                    f'{output}/baf/tumor.1bed',
                    '-O',
                    f'{output}/rdr',
                    '--segfile',
                    config.count_reads.segfile,
                    '--chromosomes',
                ]
                + (chromosomes or [])  # important to keep this as a list here to allow proper argparse parsing
                + extra_args
            )

        # ----------------------------------------------------

        if config.run.combine_counts:
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
                    level='INFO',
                )

            combine_counts(args)

    else:
        # ----------------------------------------------------
        # ----------------------------------------------------
        # Old fixed-width binning

        if config.run.count_reads:
            os.makedirs(f'{output}/rdr', exist_ok=True)
            count_reads_fw(
                args=[
                    '-N',
                    config.run.normal,
                    '-g',
                    config.run.reference,
                    '-T',
                ]
                + config.run.bams.split()
                + ['-b', config.count_reads_fw.size, '-S']
                + ('Normal ' + config.run.samples).split()
                + [
                    '-O',
                    f'{output}/rdr/normal.1bed',
                    '-o',
                    f'{output}/rdr/tumor.1bed',
                    '-t',
                    f'{output}/rdr/total.tsv',
                    '--chromosomes',
                ]
                + (chromosomes or [])  # important to keep this as a list here to allow proper argparse parsing
                + extra_args
            )

        # ----------------------------------------------------

        if config.run.combine_counts:
            _stdout = sys.stdout
            sys.stdout = StringIO()

            combine_counts_fw(
                args=[
                    '-c',
                    f'{output}/rdr/normal.1bed',
                    '-C',
                    f'{output}/rdr/tumor.1bed',
                    '-B',
                    f'{output}/baf/tumor.1bed',
                    '-t',
                    f'{output}/rdr/total.tsv',
                ]
            )
            out = sys.stdout.getvalue()
            sys.stdout.close()
            sys.stdout = _stdout

            os.makedirs(f'{output}/bb', exist_ok=True)
            with open(f'{output}/bb/bulk.bb', 'w') as f:
                f.write(out)

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
