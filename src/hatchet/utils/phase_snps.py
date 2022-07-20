import os.path
import subprocess as pr
import glob
import shlex
import shutil

import hatchet.utils.ArgParsing as ap
from hatchet.utils.Supporting import log, logArgs, error, ensure, run
from hatchet.utils.multiprocessing import Worker


def main(args=None):
    log(msg='# log notes\n', level='STEP')
    args = ap.parse_phase_snps_arguments(args)
    logArgs(args, 80)

    bcftools = args['bcftools']
    shapeit = args['shapeit']
    picard = args['picard']
    bgzip = args['bgzip']
    outdir = args['outdir']

    rpd = args['refpaneldir']
    panel = os.path.join(rpd, '1000GP_Phase3')

    # path to hg19, 1000GP in hg19 coords, potentially needed for liftover
    hg19_path = ''
    # chain files for liftover, chains['hg38_hg19']=path, chains['hg19_hg38']=path
    chains = ''
    # file for renaming chrs with bcftools, rename_files[0] for removing 'chr, rename_files[1] for adding 'chr'
    rename_files = ''

    if args['refvers'] == 'hg38':
        # Download reference panel genome and chain files
        hg19_path = os.path.join(rpd, 'hg19_no_chr.fa')
        if args['chrnot']:
            chains = {
                'hg38_hg19': os.path.join(rpd, 'hg38ToHg19.chr.chain'),
                'hg19_hg38': os.path.join(rpd, 'hg19ToHg38.chr.chain'),
            }
        else:
            chains = {
                'hg38_hg19': os.path.join(rpd, 'hg38ToHg19.no_chr.chain'),
                'hg19_hg38': os.path.join(rpd, 'hg19ToHg38.no_chr.chain'),
            }

        ensure(
            os.path.isfile(chains['hg38_hg19']) and os.path.isfile(chains['hg19_hg38']),
            'The appropriate liftover chain files could not be located! Please run the download-panel '
            'command that downloads these',
        )

    elif args['refvers'] == 'hg19' and args['chrnot']:
        rename_files = [os.path.join(rpd, f'rename_chrs{i}.txt') for i in range(1, 3)]

    # liftover VCFs, phase, liftover again to original coordinates
    os.makedirs(outdir, exist_ok=True)

    chromosomes = []
    for chro in args['chromosomes']:
        if chro.endswith('X') or chro.endswith('Y'):
            log(
                msg=f'Skipping chromosome {chro} (because it ends with X or Y)\n',
                level='WARN',
            )
        else:
            chromosomes.append(chro)

    n_instances = min(max(1, int(args['j'])), len(chromosomes))
    phaser = Phaser(
        panel,
        outdir=outdir,
        hg19=hg19_path,
        ref=args['refgenome'],
        chains=chains,
        rename=rename_files,
        refvers=args['refvers'],
        chrnot=args['chrnot'],
        verbose=False,
        bcftools=bcftools,
        shapeit=shapeit,
        picard=picard,
        bgzip=bgzip,
    )

    _snplist = args['snps']
    phaser_args = [(_snplist[chro], chro) for chro in chromosomes]
    vcfs = phaser.run(work=phaser_args, n_instances=n_instances)
    concat(vcfs, outdir=outdir, bcftools=bcftools)

    # read shapeit output, print fraction of phased snps per chromosome
    print_log(path=outdir, chromosomes=chromosomes)
    cleanup(outdir)


def cleanup(outdir):
    f = []
    # shapeit logs
    [f.extend(glob.glob(f'shapeit*{ext}')) for ext in ['.log', '.mm', '.hh']]
    # intermediate files
    exts = [
        '_phased.vcf.gz',
        '_phased.vcf.gz.csi',
        '_filtered.vcf.gz',
        '_rejected.vcf.gz',
        '_lifted.vcf.gz',
        '_lifted.vcf.gz.tbi',
        '_toFilter.vcf.gz',
        '_toConcat.vcf.gz',
        '_toConcat.vcf.gz.csi',
        '.haps',
        '.sample',
        '_alignments.snp.strand',
        '_alignments.snp.strand.exclude',
    ]
    [f.append(os.path.join(outdir, f'{c}{e}')) for c in range(1, 23) for e in exts]
    [os.remove(i) for i in f if os.path.isfile(i)]


def print_log(path, chromosomes):
    out = open(os.path.join(path, 'phased.log'), 'w')
    print(
        'chrom',
        'phased_snps',
        'original_snps',
        'proportion',
        file=out,
        sep='\t',
    )
    for c in chromosomes:
        for l in open(os.path.join(path, f'{c}_alignments.log'), 'r'):
            if 'SNPs included' in l:
                snps = int(l.split()[1])
            elif 'reference panel sites included' in l:
                phased_snps = int(l.split()[1])
        print(c, phased_snps, snps, float(phased_snps / snps), file=out, sep='\t')


def concat(vcfs, outdir, bcftools):
    infiles = ' '.join(vcfs)
    outfile = os.path.join(outdir, 'phased.vcf.gz')
    run(
        f'{bcftools} concat --output-type z --output {outfile} {infiles}',
        stdouterr_filepath=os.path.join(outdir, 'concat.log'),
    )
    return outfile


class Phaser(Worker):
    def __init__(
        self,
        panel,
        outdir,
        hg19,
        ref,
        chains,
        rename,
        refvers,
        chrnot,
        verbose,
        bcftools,
        shapeit,
        picard,
        bgzip,
    ):
        self.panel = panel
        self.outdir = outdir
        self.hg19 = hg19
        self.ref = ref
        self.chains = chains
        self.rename = rename
        self.refvers = refvers
        self.chrnot = chrnot
        self.verbose = verbose
        self.bcftools = bcftools
        self.shapeit = shapeit
        self.picard = picard
        self.bgzip = bgzip

    def work(self, *args):
        vcf, chromosome = args

        # (1) PREPROCESS
        if self.refvers == 'hg19':
            # no need for liftover, just deal with chr annotation
            if self.chrnot:
                vcf_toFilter = self.change_chr(
                    infile=vcf,
                    chromosome=chromosome,
                    outname='toFilter',
                    rename=self.rename[0],
                )
            else:
                # just copy files, vcfs already appropriately formatted
                vcf_toFilter = self.stage_vcfs(infile=vcf, chromosome=chromosome)
        else:
            # liftover
            vcf_toFilter = self.liftover(
                infile=vcf,
                chromosome=chromosome,
                outname='toFilter',
                chain=self.chains['hg38_hg19'],
                refgen=self.hg19,
                ch=False,
            )

        # (2) FILTERING AND PHASING
        # filter out multi-allelic sites and indels
        vcf_filtered = self.biallelic(infile=vcf_toFilter, chromosome=chromosome)
        vcf_phased = self.run_shapeit(infile=vcf_filtered, chromosome=chromosome)  # phase

        # (3) POSTPROCESS
        if self.refvers == 'hg19':
            if self.chrnot:
                vcf_to_concat = self.change_chr(
                    infile=vcf_phased,
                    chromosome=chromosome,
                    outname='toConcat',
                    rename=self.rename[1],
                )
                # re-index with bcftools after renaming
                self.index(infile=vcf_to_concat, chromosome=chromosome)
            else:
                vcf_to_concat = vcf_phased  # do nothing; vcfs already in original format
        else:
            vcf_to_concat = self.liftover(
                infile=vcf_phased,
                chromosome=chromosome,
                outname='toConcat',
                chain=self.chains['hg19_hg38'],
                refgen=self.ref,
                ch=self.chrnot,
            )

        return vcf_to_concat

    def liftover(self, infile, chromosome, outname, chain, refgen, ch):
        errname = os.path.join(self.outdir, f'{chromosome}_picard.log')
        # output from picard liftover, to be filtered
        tmpfile = os.path.join(self.outdir, f'{chromosome}_lifted.vcf.gz')
        # filtered with bcftools
        outfile = os.path.join(self.outdir, f'{chromosome}_{outname}.vcf.gz')
        # required file of SNPs that didn't liftover
        rejfile = os.path.join(self.outdir, f'{chromosome}_rejected.vcf.gz')

        # --WARN_ON_MISSING_CONTIG true: throws out liftovers to contigs not present in target reference,
        # e.g. small contigs variably present among the assemblies

        cmd1 = (
            f'{self.picard} LiftoverVcf -I={infile} -O={tmpfile} -CHAIN={chain} -R={refgen} -REJECT={rejfile} '
            '--WARN_ON_MISSING_CONTIG=true'
        )
        # need to change 'chr' notation depending on liftover direction
        c = chromosome if not ch else f'chr{chromosome}'
        # filter out mapping to other chromosomes/contigs!
        cmd2 = f'{self.bcftools} filter --output-type z --regions {c} {tmpfile}'
        # remove duplicate sites from liftover
        cmd3 = f'{self.bcftools} norm --remove-duplicates --output {outfile}'

        with open(errname, 'w') as err:
            pic = pr.run(cmd1.split(), stdout=err, stderr=err, universal_newlines=True)
            filt = pr.Popen(
                shlex.split(cmd2),
                stdout=pr.PIPE,
                stderr=err,
                universal_newlines=True,
            )
            norm = pr.Popen(
                shlex.split(cmd3),
                stdin=filt.stdout,
                stdout=err,
                stderr=err,
                universal_newlines=True,
            )
            codes = list(map(lambda p: p.wait(), [filt, norm]))
        if any(c != 0 for c in codes) or pic.returncode != 0:
            log(msg=cmd1, level='ERROR')
            error(
                f'Failed to liftover chromosomes with picard on {infile}.',
                raise_exception=True,
            )
        else:
            os.remove(errname)
        return outfile

    def change_chr(self, infile, chromosome, outname, rename):
        # use bcftools to rename chromosomes
        outfile = os.path.join(self.outdir, f'{chromosome}_{outname}.vcf.gz')
        errname = os.path.join(self.outdir, f'{chromosome}_bcftools.log')

        run(
            f'{self.bcftools} annotate --rename-chrs {rename} --output-type z --output {outfile} {infile}',
            stdouterr_filepath=errname,
            error_msg=f'Failed to re-annotate chromosomes with bcftools on {infile}.',
        )
        return outfile

    def stage_vcfs(self, infile, chromosome):
        # copy file, so that phasing takes same input file name regardless of conditions, unzip
        outfile = os.path.join(self.outdir, f'{chromosome}_toFilter.vcf.gz')
        shutil.copyfile(infile, outfile)
        return outfile

    def biallelic(self, infile, chromosome):
        # use bcftools to discard multi-allelic sites and indels
        outfile = os.path.join(self.outdir, f'{chromosome}_filtered.vcf.gz')
        run(
            (
                f'{self.bcftools} view --max-alleles 2 --exclude-types indels --output-type z --output-file {outfile} '
                f'{infile}'
            ),
            stdouterr_filepath=os.path.join(self.outdir, f'{chromosome}_bcftools.log'),
            error_msg=f'Biallelic sites filtering failed on {infile}.',
        )
        return outfile

    def run_shapeit(self, infile, chromosome):
        # use shapeit with reference panel to phase vcf files
        errname = os.path.join(self.outdir, f'{chromosome}_shapeit.log')

        # define params used across shapeit functions
        inmap = f'{self.panel}/genetic_map_chr{chromosome}_combined_b37.txt'
        inref = (
            f'{self.panel}/1000GP_Phase3_chr{chromosome}.hap.gz '
            f'{self.panel}/1000GP_Phase3_chr{chromosome}.legend.gz '
            f'{self.panel}/1000GP_Phase3.sample'
        )

        # check data with shapeit -check; get list of sites to exclude,
        # such as sites in target VCF that are not present in reference panel
        cmd_check = (
            f'{self.shapeit} -check --input-vcf {self.outdir}/{chromosome}_filtered.vcf.gz '
            f'--input-map {inmap} '
            f'--input-ref {inref} '
            f'--output-log {self.outdir}/{chromosome}_alignments'
        )

        cmd_phase = (
            f'{self.shapeit} --input-vcf {self.outdir}/{chromosome}_filtered.vcf.gz '
            f'--input-map {inmap} '
            f'--input-ref {inref} '
            f'--exclude-snp {self.outdir}/{chromosome}_alignments.snp.strand.exclude '
            f'--output-max {self.outdir}/{chromosome}.haps {self.outdir}/{chromosome}.sample '
            '--chrX --no-mcmc '
            '--seed 0'
        )

        cmd_convert = (
            f'{self.shapeit} -convert --input-haps {self.outdir}/{chromosome} '
            f'--output-vcf {self.outdir}/{chromosome}_phased.vcf'
        )

        cmd_compress = f'{self.bgzip} -f {self.outdir}/{chromosome}_phased.vcf'

        cmd_index = f'{self.bcftools} index -f {self.outdir}/{chromosome}_phased.vcf.gz'

        run(cmd_check, stdouterr_filepath=errname, check_return_codes=False)  # May return 1
        run(
            [cmd_phase, cmd_convert, cmd_compress, cmd_index],
            stdouterr_filepath=errname,
            error_msg=f'Phasing failed on {infile}.',
        )

        return os.path.join(self.outdir, f'{chromosome}_phased.vcf.gz')

    def index(self, infile, chromosome):
        run(
            f'{self.bcftools} index -f {infile}',
            stdouterr_filepath=os.path.join(self.outdir, f'{chromosome}_bcftools.log'),
            error_msg=f'Failed to index {infile} with bcftools.',
        )


if __name__ == '__main__':
    main()
