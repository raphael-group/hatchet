import sys
from io import StringIO
import os.path
import glob
import argparse

from hatchet import config
from hatchet.utils.count_reads import main as count_reads
from hatchet.utils.genotype_snps import main as genotype_snps
from hatchet.utils.count_alleles import main as count_alleles
from hatchet.utils.combine_counts import main as combine_counts
from hatchet.utils.cluster_bins import main as cluster_bins
from hatchet.utils.plot_bins import main as plot_bins
from hatchet.bin.HATCHet import main as hatchet_main
from hatchet.utils.plot_cn import main as plot_cn


def main(args=None):

    parser = argparse.ArgumentParser(prog='hatchet run', description='Run HATCHet pipeline')
    parser.add_argument("inifile", help=".ini file for run configuration")
    args = parser.parse_args(args)

    config.read(args.inifile)

    output = config.run.output
    output = output.rstrip('/')
    os.makedirs(output, exist_ok=True)

    extra_args = []
    try:
        if config.run.processes is not None:
            extra_args = ['-j', str(config.run.processes)]
    except KeyError:
        pass

    # ----------------------------------------------------

    if config.run.count_reads:
        os.makedirs(f'{output}/rdr', exist_ok=True)
        count_reads(
            args=[
                '-N', config.run.normal,
                '-g', config.run.reference,
                '-T'
            ] + config.run.bams.split() + [
                '-b', config.count_reads.size,
                '-S'
            ] + ('Normal ' + config.run.samples).split() + [
                '-O', f'{output}/rdr/normal.1bed',
                '-o', f'{output}/rdr/tumor.1bed',
                '-t', f'{output}/rdr/total.tsv'
            ] + extra_args
        )

    # ----------------------------------------------------

    if config.run.genotype_snps:
        snps = ''
        if config.genotype_snps.snps:
           snps = config.genotype_snps.snps
        elif config.genotype_snps.reference_version:
            snps_mapping = {
                ('hg19', True): 'https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/GATK/00-All.vcf.gz',
                ('hg19', False): 'https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/00-All.vcf.gz',
                ('hg38', True): 'https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/GATK/00-All.vcf.gz',
                ('hg38', False): 'https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/00-All.vcf.gz'
            }

            if (config.genotype_snps.reference_version, config.genotype_snps.chr_notation) not in snps_mapping:
                raise RuntimeError('Please specify valid values of reference_version and chr_notation')
            else:
                snps = snps_mapping[config.genotype_snps.reference_version, config.genotype_snps.chr_notation]

        os.makedirs(f'{output}/snps', exist_ok=True)
        genotype_snps(
            args=[
                '-N', config.run.normal,
                '-r', config.run.reference,
                '-R', snps,
                '-o', f'{output}/snps/'
            ] + extra_args
        )

    # ----------------------------------------------------

    if config.run.count_alleles:
        os.makedirs(f'{output}/baf', exist_ok=True)
        count_alleles(
            args=[
                '-N', config.run.normal,
                '-T'
            ] + config.run.bams.split() + [
                '-S'
            ] + ('Normal ' + config.run.samples).split() + [
                '-r', config.run.reference,
                '-L'
            ] + glob.glob(f'{output}/snps/*.vcf.gz') + [
                '-O', f'{output}/baf/normal.1bed',
                '-o', f'{output}/baf/tumor.1bed',
                '-l', f'{output}'
            ] + extra_args
        )

    # ----------------------------------------------------

    if config.run.combine_counts:
        _stdout = sys.stdout
        sys.stdout = StringIO()

        combine_counts(args=[
            '-c', f'{output}/rdr/normal.1bed',
            '-C', f'{output}/rdr/tumor.1bed',
            '-B', f'{output}/baf/tumor.1bed',
            '-t', f'{output}/rdr/total.tsv'
        ])

        out = sys.stdout.getvalue()
        sys.stdout.close()
        sys.stdout = _stdout

        os.makedirs(f'{output}/bb', exist_ok=True)
        with open(f'{output}/bb/bulk.bb', 'w') as f:
            f.write(out)

    # ----------------------------------------------------

    if config.run.cluster_bins:
        os.makedirs(f'{output}/bbc', exist_ok=True)
        cluster_bins(args=[
            f'{output}/bb/bulk.bb',
            '-o', f'{output}/bbc/bulk.seg',
            '-O', f'{output}/bbc/bulk.bbc'
        ])

    # ----------------------------------------------------

    if config.run.plot_bins:
        os.makedirs(f'{output}/plots', exist_ok=True)
        plot_bins(args=[
            f'{output}/bbc/bulk.bbc',
            '--rundir', f'{output}/plots'
        ])

    # ----------------------------------------------------

    if config.run.compute_cn:
        os.makedirs(f'{output}/results', exist_ok=True)
        hatchet_main(
            args=[
                '-x', f'{output}/results',
                '-i', f'{output}/bbc/bulk'
            ] + extra_args
        )

    # ----------------------------------------------------

    if config.run.plot_cn:
        os.makedirs(f'{output}/summary', exist_ok=True)
        plot_cn(args=[
            f'{output}/results/best.bbc.ucn',
            '--rundir', f'{output}/summary'
        ])


if __name__ == '__main__':
    main()
