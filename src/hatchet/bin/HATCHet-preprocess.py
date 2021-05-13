#!/usr/bin/env python3



import sys, os
import argparse
import subprocess as sp
import multiprocessing as mp
import shlex
import re

utils = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, 'utils')
if not os.path.isdir(utils):
    raise ValueError("utils directory not found in HATCHet's home directory i.e. {}, is anything been moved?".format(utils))
sys.path.append(utils)
from Supporting import *
from hatchet import config, __version__


def parse_args():
    description = "This command automatically runs the HATCHet's preprocessing pipeline, which is composed of three steps: (1) count-reads, (2) count-alleles, and (3) combine-counts."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-t","--tumor", required=True, type=str, help="White-space separated list of input tumor BAM files, corresponding to multiple samples from the same patient (list must be within quotes)")
    parser.add_argument("-n","--normal", required=True, type=str, help="Matched-normal BAM file")
    parser.add_argument("-r","--reference", type=str, required=True, help="Reference genome")
    parser.add_argument("-s","--samplenames", required=False, type=str, default=config.preprocess.samplenames, help="Tumor sample names in a white-space separated list in the same order as the corresponding BAM files (default: file names are used as names)")
    parser.add_argument("-b","--size", type=str, required=False, default=config.preprocess.size, help="Bin size, with or without \"kb\" or \"Mb\" (default: 250kb)")
    parser.add_argument("-c","--minreads", type=int, required=False, default=config.preprocess.minreads, help="Minimum read counts for heterozygous germline SNPs (default: 8)")
    parser.add_argument("-C","--maxreads", type=int, required=False, default=config.preprocess.maxreads, help="Maximum read counts for heterozygous germline SNPs (default: 1000)")
    parser.add_argument("-p","--phred", type=int, required=False, default=config.preprocess.phred, help="Phred quality score (default: 11)")
    parser.add_argument("-x","--rundir", required=False, default=config.preprocess.rundir, type=str, help="Running directory (default: current directory)")
    parser.add_argument("--bcftools", required=False, default=config.paths.bcftools, type=str, help="Path to the directory to \"bcftools\" executable, required in default mode (default: bcftools is directly called as it is in user $PATH)")
    parser.add_argument("--samtools", required=False, default=config.paths.samtools, type=str, help="Path to the directory to \"samtools\" executable, required in default mode (default: samtools is directly called as it is in user $PATH)")
    parser.add_argument("--seed", required=False, type=int, default=config.preprocess.seed, help="Random seed for replication (default: None)")
    parser.add_argument("-j","--jobs", required=False, type=int, default=config.preprocess.jobs, help="Number of parallele jobs to use (default: equal to number of available processors)")
    parser.add_argument("-V", "--version", action='version', version=f'%(prog)s {__version__}')
    args = parser.parse_args()

    tumor = set(t for t in args.tumor.split())
    for t in tumor:
        if not os.path.isfile(t):
            raise ValueError("The following BAM file does not exist: {}".format(t))
    if args.samplenames is None:
        names = set(os.path.splitext(os.path.basename(t))[0] for t in tumor)
        if len(names) != len(tumor):
            names = tumor
    else:
        names = set(t for t in args.samplenames.split())
        if len(names) != len(tumor):
            raise ValueError("A different number of samples names has been provided compared to the number of BAM files, remember to add the list within quotes!")
        
    tumor = set(os.path.abspath(t) for t in tumor)
    if not os.path.isdir(args.rundir):
        raise ValueError("Running directory does not exists: {}".format(args.rundir))
    if not os.path.isfile(args.normal):
        raise ValueError("Matched-normal BAM file does not exist: {}".format(args.normal))
    if not os.path.isfile(args.reference):
        raise ValueError("Reference genome file does not exist: {}".format(args.reference))
    if args.seed and args.seed < 1:
        raise ValueError("The random seed  must be positive!")
    if args.minreads < 1:
        raise ValueError("The minimum number of reads must be positive!")
    if args.maxreads < 1:
        raise ValueError("The maximum number of reads must be positive!")

    size = 0
    try:
        if args.size[-2:] == "kb":
            size = int(args.size[:-2]) * 1000
        elif args.size[-2:] == "Mb":
            size = int(args.size[:-2]) * 1000000
        else:
            size = int(args.size)
    except:
        raise ValueError("Size must be a number, optionally ending with either \"kb\" or \"Mb\"!")

    if not args.jobs:
        args.jobs = mp.cpu_count()
    if args.jobs < 1:
        raise ValueError("The number of jobs must be positive!")

    return {
        "tumor" : list(tumor),
        "normal" : os.path.abspath(args.normal),
        "ref" : os.path.abspath(args.reference),
        "names" : list(names),
        "size" : size,
        "samtools" : args.samtools,
        "bcftools" : args.bcftools,
        "J" : args.jobs,
        "minreads" : args.minreads,
        "maxreads" : args.maxreads,
        "phred" : args.phred,
        "rundir" : os.path.abspath(args.rundir),
        "seed" : args.seed
    }


def main():
    log('Parsing and checking arguments\n', level='PROGRESS')
    args = parse_args()
    log('\n'.join(['Arguments:'] + ['\t{} : {}'.format(a, args[a]) for a in args]) + '\n', level='INFO')

    log('Setting directories\n', level='PROGRESS')
    dbaf, drdr, dbb = setup(args)
    def get_comp(name):
        comp = os.path.join(utils, name)
        if not os.path.isfile(comp):
            raise ValueError("{} not found in utils directory of HATCHet's home i.e. {}, is anything been moved?".format(name, utils))
        return comp

    log('Computing RDRs\n', level='PROGRESS')
    cmd = 'python3 {} -N {} -T {} -S {} -b {} -g {} -j {} -q {} -O {} -o {}'
    nbin = os.path.join(drdr, 'normal.1bed')
    tbin = os.path.join(drdr, 'bulk.1bed')
    cmd = cmd.format(get_comp('count_reads.py'), args['normal'], ' '.join(args['tumor']), 'normal ' + ' '.join(args['names']), args['size'], args['ref'], args['J'], args['phred'], nbin, tbin)
    if args['samtools'] is not None:
        cmd += " --samtools {}".format(args['samtools'])
    runcmd(cmd, drdr, log="bins.log", rundir=args['rundir'])

    log('Computing BAFs\n', level='PROGRESS')
    cmd = 'python3 {} -N {} -T {} -S {} -r {} -j {} -q {} -Q {} -U {} -c {} -C {} -O {} -o {}'
    nbaf = os.path.join(dbaf, 'normal.1bed')
    tbaf = os.path.join(dbaf, 'bulk.1bed')
    cmd = cmd.format(get_comp('count_alleles.py'), args['normal'], ' '.join(args['tumor']), 'normal ' + ' '.join(args['names']), args['ref'], args['J'], args['phred'], args['phred'], args['phred'], args['minreads'], args['maxreads'], nbaf, tbaf)
    if args['samtools'] is not None:
        cmd += " --samtools {}".format(args['samtools'])
    if args['bcftools'] is not None:
        cmd += " --bcftools {}".format(args['bcftools'])
    runcmd(cmd, dbaf, log="bafs.log", rundir=args['rundir'])

    log('Combining RDRs and BAFs\n', level='PROGRESS')
    ctot = os.path.join(args['rundir'], 'total_read.counts')
    cmd = 'python3 {} -c {} -C {} -B {} -m MIRROR -t {}'
    cmd = cmd.format(get_comp('combine_counts.py'), nbin, tbin, tbaf, ctot)
    if args['seed'] is not None:
        cmd += " -e {}".format(args['seed'])
    runcmd(cmd, dbb, out='bulk.bb', log="combo.log", rundir=args['rundir'])


def setup(args):
    dbaf = os.path.join(args['rundir'], 'baf')
    if not os.path.isdir(dbaf):
        os.mkdir(dbaf)

    drdr = os.path.join(args['rundir'], 'rdr')
    if not os.path.isdir(drdr):
        os.mkdir(drdr)

    dbb = os.path.join(args['rundir'], 'bb')
    if not os.path.isdir(dbb):
        os.mkdir(dbb)

    return dbaf, drdr, dbb


def runcmd(cmd, xdir, out=None, log="log", rundir=None):
    j = os.path.join
    tmp = log + '_TMP'
    sout = open(j(xdir, out), 'w') if out is not None else sp.PIPE
    with open(j(xdir, tmp), 'w') as serr:
        proc = sp.Popen(shlex.split(cmd), stdout=sout, stderr=sp.PIPE, cwd=rundir, universal_newlines=True)
        for line in iter(lambda : proc.stderr.read(1), ''):
            sys.stderr.write(line)
            serr.write(line)
    if out is not None:
        sout.flush()
        sout.close()

    with open(j(xdir, tmp), 'r') as i:
        with open(j(xdir, log), 'w') as o:
            for l in i:
                if 'Progress' not in l:
                    o.write(re.sub(r'\033\[[0-9]*m', '', l))
    os.remove(j(xdir, tmp))


if __name__ == '__main__':
    main()
