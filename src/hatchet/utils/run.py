import sys
from io import StringIO
import os.path
import glob
import argparse

import hatchet
from hatchet import config
from hatchet.utils.binBAM import main as binBAM
from hatchet.utils.SNPCaller import main as SNPCaller
from hatchet.utils.deBAF import main as deBAF
from hatchet.utils.comBBo import main as comBBo
from hatchet.utils.cluBB import main as cluBB
from hatchet.utils.BBot import main as BBot
from hatchet.bin.HATCHet import main as solve
from hatchet.utils.BBeval import main as BBeval

solve_binary = os.path.join(os.path.dirname(hatchet.__file__), 'solve')


def main(args=None):

    parser = argparse.ArgumentParser(prog='hatchet run', description='Run HATCHet pipeline')
    parser.add_argument("inifile", help=".ini file for run configuration")
    args = parser.parse_args(args)

    config.read(args.inifile)

    output = config.paths.output
    output.rstrip('/')
    os.makedirs(output, exist_ok=True)

    # ----------------------------------------------------

    if config.bin.enabled:
        os.makedirs(f'{output}/rdr', exist_ok=True)
        binBAM(
            args=[
                '-N', config.paths.normal,
                '-T'
            ] + config.paths.bams.split() + [
                '-b', config.bin.size,
                '-S'
            ] + ('Normal ' + config.bin.samples).split() + [
                '-O', f'{output}/rdr/normal.1bed',
                '-o', f'{output}/rdr/tumor.1bed',
                '-t', f'{output}/rdr/total.tsv'
            ]
        )

    # ----------------------------------------------------

    if config.snp.enabled:
        snps = ''
        if config.snp.reference_version and not config.snp.snps:
            snps_mapping = {
                ('hg19', True): 'https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/GATK/00-All.vcf.gz',
                ('hg19', False): 'https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/00-All.vcf.gz',
                ('hg38', True): 'https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/GATK/00-All.vcf.gz',
                ('hg38', False): 'https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/00-All.vcf.gz'
            }

            if (config.snp.reference_version, config.snp.chr_notation) not in snps_mapping:
                raise RuntimeError('Please specify valid values of reference_version and chr_notation')
            else:
                snps = snps_mapping[config.snp.reference_version, config.snp.chr_notation]

        os.makedirs(f'{output}/snps', exist_ok=True)
        SNPCaller(
            args=[
                '-N', config.paths.normal,
                '-r', config.paths.reference,
                '-R', snps,
                '-o', f'{output}/snps/'
            ]
        )

    # ----------------------------------------------------

    if config.baf.enabled:
        os.makedirs('baf', exist_ok=True)
        deBAF(
            args=[
                '-N', config.paths.normal,
                '-T'
            ] + config.paths.bams.split() + [
                '-S'
            ] + ('Normal ' + config.bin.samples).split() + [
                '-r', config.paths.reference,
                '-L'
            ] + glob.glob(f'{output}/snps/*.vcf.gz') + [
                '-O', f'{output}/baf/normal.1bed',
                '-o', f'{output}/baf/tumor.1bed'
            ]
        )

    # ----------------------------------------------------

    if config.combbo.enabled:
        _stdout = sys.stdout
        sys.stdout = StringIO()

        comBBo(args=[
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

    if config.clubb.enabled:
        os.makedirs(f'{output}/bbc', exist_ok=True)
        cluBB(args=[
            f'{output}/bb/bulk.bb',
            '-o', f'{output}/bbc/bulk.seg',
            '-O', f'{output}/bbc/bulk.bbc'
        ])

    # ----------------------------------------------------

    if config.bbot.enabled:
        os.makedirs(f'{output}/plots', exist_ok=True)
        BBot(args=[
            f'{output}/bbc/bulk.bbc',
            '--rundir', 'plots'
        ])

    # ----------------------------------------------------

    if config.solver.enabled:
        os.makedirs(f'{output}/results', exist_ok=True)
        solve(args=[
            solve_binary,
            '-x', f'{output}/results',
            '-i', f'{output}/bbc/bulk',
            '-g', '0.35',
            '-l', '0.6'
        ])

    # ----------------------------------------------------

    if config.bbeval.enabled:
        os.makedirs(f'{output}/summary', exist_ok=True)
        BBeval(args=[
            f'{output}/results/best.bbc.ucn',
            '--rundir', 'summary'
        ])


if __name__ == '__main__':
    main()
