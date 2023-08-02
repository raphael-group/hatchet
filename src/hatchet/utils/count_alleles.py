import sys
import os
import os.path
import shlex
import subprocess as pr
from scipy.stats import beta
import tempfile

import hatchet.utils.ArgParsing as ap
from hatchet.utils.Supporting import log, logArgs, error, close
from hatchet.utils.multiprocessing import Worker


def main(args=None):
    log(
        msg=(
            '# Parsing the input arguments, checking the consistency of given files, and extracting required '
            'information\n'
        ),
        level='STEP',
    )
    args = ap.parse_count_alleles_arguments(args)
    logArgs(args, 80)

    log(msg='# Extracting heterozygous SNPs\n', level='STEP')

    hetSNPs = selectHetSNPs(
        bcftools=args['bcftools'],
        normal=args['normal'],
        chromosomes=args['chromosomes'],
        num_workers=args['j'],
        snplist=args['snps'],
        mincov=args['mincov'],
        dp=args['maxcov'],
        verbose=args['verbose'],
        outdir=args['outputSnps'],
    )

    if not hetSNPs:
        close('No heterozygous SNPs found in the selected regions of the normal!\n')

    log(
        msg='# Writing the list of selected SNPs, covered and heterozygous in the normal sample\n',
        level='STEP',
    )
    # put lists of SNPs in uniquely named temporary dir to prevent concurrent hatchet runs from overwriting
    # each others temp files
    with tempfile.TemporaryDirectory(dir=args['outputSnps']) as tmpdirname:
        hetsnpsfiles = {}
        for chro in args['chromosomes']:
            hetsnpsfiles[chro] = os.path.join(tmpdirname, 'TMP_{}.tsv'.format(chro))
            with open(hetsnpsfiles[chro], 'w') as f:
                if (args['normal'][1], chro) in hetSNPs:
                    for snp in sorted(hetSNPs[args['normal'][1], chro]):
                        f.write('{}\t{}\n'.format(chro, snp))

        log(
            msg='# Writing the allele counts of the normal sample for selected SNPs\n',
            level='STEP',
        )
        handle = open(args['outputNormal'], 'w') if args['outputNormal'] is not None else sys.stdout
        for chro in args['chromosomes']:
            if (args['normal'][1], chro) in hetSNPs:
                for snp in sorted(hetSNPs[args['normal'][1], chro]):
                    count = hetSNPs[args['normal'][1], chro][snp]
                    handle.write(
                        '{}\t{}\t{}\t{}\t{}\n'.format(
                            chro,
                            snp,
                            args['normal'][1],
                            count[0][1],
                            count[1][1],
                        )
                    )
        if handle is not sys.stdout:
            handle.close()

        log(msg='# Counting SNPs alleles from tumour samples\n', level='STEP')
        rcounts = counting(
            bcftools=args['bcftools'],
            reference=args['reference'],
            samples=args['samples'],
            chromosomes=args['chromosomes'],
            num_workers=args['j'],
            snplist=hetsnpsfiles,
            q=args['q'],
            Q=args['Q'],
            mincov=args['mincov'],
            dp=args['maxcov'],
            E=args['E'],
            verbose=args['verbose'],
            outdir=args['outputSnps'],
        )
        if not rcounts:
            close('The selected SNPs are not covered in the tumors!\n')
        rcounts = {c: dict(map(lambda r: (int(r[2]), dict(r[3])), rcounts[c])) for c in rcounts}
        het = lambda chro: hetSNPs[args['normal'][1], chro]
        form = lambda REF, ALT, T: (
            (REF, T[REF] if REF in T else 0),
            (ALT, T[ALT] if ALT in T else 0),
        )
        counts = {
            c: {o: form(het(c[1])[o][0][0], het(c[1])[o][1][0], rcounts[c][o]) for o in rcounts[c]} for c in rcounts
        }

    log(
        msg='# Writing the allele counts of tumor samples for selected SNPs\n',
        level='STEP',
    )
    handle = open(args['outputTumors'], 'w') if args['outputTumors'] is not None else sys.stdout
    for sample in args['samples']:
        for chro in args['chromosomes']:
            if (sample[1], chro) in counts:
                for snp in counts[sample[1], chro]:
                    count = counts[sample[1], chro][snp]
                    handle.write('{}\t{}\t{}\t{}\t{}\n'.format(chro, snp, sample[1], count[0][1], count[1][1]))
    if handle is not sys.stdout:
        handle.close()


def selectHetSNPs(
    bcftools,
    normal,
    chromosomes,
    num_workers,
    snplist,
    mincov,
    dp,
    verbose,
    outdir,
):
    work = []
    for chro in chromosomes:
        work.append((normal[1], chro))

    worker = HetSNPSelector(
        bcftools,
        mincov,
        dp,
        snplist,
        verbose,
        outdir,
    )

    results = worker.run(work=work, n_instances=num_workers)
    results = {(v[0][0], v[0][1]): v for v in results if v}
    results = {
        c: [(r[0], r[1], r[2], r[3][0], max(r[3][1:], key=(lambda x: x[1]))) for r in results[c] if len(r[3]) >= 2]
        for c in results
    }
    return {c: {int(r[2]): (r[3], r[4]) for r in reversed(results[c])} for c in results if len(results[c]) > 0}

class HetSNPSelector(Worker):
    def __init__(
        self,
        bcftools,
        mincov,
        dp,
        snplist,
        verbose,
        outdir,
    ):
        self.bcftools = bcftools
        self.mincov = mincov
        self.dp = dp
        self.snplist = snplist
        self.verbose = verbose
        self.outdir = outdir

    def work(self, samplename, chromosome):
        return self.selectHets(samplename=samplename, chromosome=chromosome)

    def selectHets(self, samplename, chromosome):
        # query the genotype file for heterozygous SNPs. no need for recounting alleles for normal
        vcffile = self.snplist[chromosome]
        cmd_query = "{} query -f '%CHROM\\t%POS\\t%REF,%ALT\\t%AD\\n' -i 'SUM(FMT/AD)<={} & SUM(FMT/AD)>={} & FMT/GT[0]==\"0/1\"' {}".format(
            self.bcftools, self.dp, self.mincov, vcffile
        )
        errname = os.path.join(self.outdir, '{}_{}_bcftools.log'.format(samplename, chromosome))
        with open(errname, 'w') as err:
            query = pr.Popen(
                shlex.split(cmd_query),
                stdout=pr.PIPE,
                stderr=err,
                universal_newlines=True,
            )
            stdout, _ = query.communicate()
            codes = map(lambda p: p.wait(), [query])
        if any(c != 0 for c in codes):
            raise ValueError(
                error(
                    'Allele counting failed on {} of {}, please check errors in {}!'.format(
                        chromosome, samplename, errname
                    )
                )
            )
        else:
            os.remove(errname)
        return count_alleles_from_stdout(stdout, samplename)

def counting(
    bcftools,
    reference,
    samples,
    chromosomes,
    num_workers,
    snplist,
    q,
    Q,
    mincov,
    dp,
    E,
    verbose,
    outdir,
):

    work = []
    for bam in samples:
        for chro in chromosomes:
            work.append((bam[0], bam[1], chro))

    worker = AlleleCounter(
        bcftools,
        reference,
        q,
        Q,
        mincov,
        dp,
        E,
        snplist,
        verbose,
        outdir,
    )

    results = worker.run(work=work, n_instances=num_workers)
    results = {(v[0][0], v[0][1]): v for v in results if v}
    return results


class AlleleCounter(Worker):
    def __init__(
        self,
        bcftools,
        reference,
        q,
        Q,
        mincov,
        dp,
        E,
        snplist,
        verbose,
        outdir,
    ):
        self.bcftools = bcftools
        self.reference = reference
        self.q = q
        self.Q = Q
        self.mincov = mincov
        self.dp = dp
        self.E = E
        self.snplist = snplist
        self.verbose = verbose
        self.outdir = outdir

    def work(self, bamfile, samplename, chromosome):
        return self.countAlleles(bamfile=bamfile, samplename=samplename, chromosome=chromosome)

    def countAlleles(self, bamfile, samplename, chromosome):
        cmd_mpileup = '{} mpileup {} -Ou -f {} --skip-indels -a INFO/AD -q {} -Q {} -d {} -T {}'.format(
            self.bcftools,
            bamfile,
            self.reference,
            self.q,
            self.Q,
            self.dp,
            self.snplist[chromosome],
        )
        cmd_query = "{} query -f '%CHROM\\t%POS\\t%REF,%ALT\\t%AD\\n' -i 'SUM(AD)<={} & SUM(AD)>={}'".format(
            self.bcftools, self.dp, self.mincov
        )
        if self.E:
            cmd_mpileup += ' -E'
        errname = os.path.join(self.outdir, '{}_{}_bcftools.log'.format(samplename, chromosome))

        with open(errname, 'w') as err:
            mpileup = pr.Popen(
                shlex.split(cmd_mpileup),
                stdout=pr.PIPE,
                stderr=err,
                universal_newlines=True,
            )
            query = pr.Popen(
                shlex.split(cmd_query),
                stdin=mpileup.stdout,
                stdout=pr.PIPE,
                stderr=err,
                universal_newlines=True,
            )
            stdout, _ = query.communicate()
            codes = map(lambda p: p.wait(), [mpileup, query])
        if any(c != 0 for c in codes):
            raise ValueError(
                error(
                    'Allele counting failed on {} of {}, please check errors in {}!'.format(
                        chromosome, samplename, errname
                    )
                )
            )
        else:
            os.remove(errname)
        return count_alleles_from_stdout(stdout, samplename)

def count_alleles_from_stdout(stdout, samplename):
        alleles = {'A', 'C', 'G', 'T'}
        mkcounts = lambda p, q: list(
            map(
                lambda y: (y[0], int(y[1])),
                filter(lambda x: x[0] in alleles, zip(p, q)),
            )
        )
        form = lambda p: (
            samplename,
            p[0],
            p[1],
            mkcounts(p[2].split(','), p[3].split(',')),
        )
        return [form(line.strip().split()) for line in stdout.strip().split('\n') if line != '']


if __name__ == '__main__':
    main()
