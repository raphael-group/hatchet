#!/usr/bin/env python3



import sys, os
import argparse
import subprocess as sp
import multiprocessing as mp
import shlex
import re
import pathlib

from .Supporting import *
from hatchet import config
from.ArgParsing import parse_preprocess_args

def main(args=None):
    log('Parsing and checking arguments\n', level='PROGRESS')
    args = parse_preprocess_args(args)
    log('\n'.join(['Arguments:'] + ['\t{} : {}'.format(a, args[a]) for a in args]) + '\n', level='INFO')

    log('Setting directories\n', level='PROGRESS')
    dbaf, drdr, dbb, dsnps, dcounts, drisk = generate_subdirectories(args)

    log("Getting read counts at risk allele sites\n", level="PROGRESS")
    cmd = '{} mpileup {} -f {} -l {} -o {} -s'
    samtools = os.path.join(args["samtools"], "samtools")
    all_samples =  [args['normal']] + args['tumor']
    all_names = ['normal'] + args['names']
    alleles_list = os.path.join(pathlib.Path(__file__).parent.parent.absolute(), 'resources', 'risk_alleles.pos')
    for i in range(len(all_samples)):
        log("Counting sample\"{}\"\n".format(all_names[i]), level="INFO")
        
        my_cmd = cmd.format(samtools, all_samples[i], args['ref'], alleles_list, 
                            os.path.join(drisk, all_names[i] + '.pileup'))
        print(my_cmd)
        runcmd(my_cmd, drisk, log = 'riskalleles.log', rundir = args['rundir'])


    log('Computing RDRs\n', level='PROGRESS')    l;
    cmd = 'python3 -m hatchet binBAM -N {} -T {} -S {} -b {} -g {} -j {} -q {} -O {} -o {}'
    nbin = os.path.join(drdr, 'normal.1bed')
    tbin = os.path.join(drdr, 'bulk.1bed')
    cmd = cmd.format(args['normal'], ' '.join(args['tumor']), 'normal ' + ' '.join(args['names']), args['size'], args['ref'], args['J'], args['phred'], nbin, tbin)
    if args['samtools'] is not None and len(args['samtools']) > 0:
        cmd += " --samtools {}".format(args['samtools'])
    print(cmd)
    runcmd(cmd, drdr, log="bins.log", rundir=args['rundir'])

    log('Calling SNPs\n', level='PROGRESS')
    cmd =  'python3 -m hatchet SNPCaller -N {} -r {} -j {} -c {} -C {} -o {}'
    # TODO: include -R reference SNPs list?
    nbin = os.path.join(drdr, 'normal.1bed')
    tbin = os.path.join(drdr, 'bulk.1bed')
    cmd = cmd.format(args['normal'], args['ref'], args['J'], args['minreads'], args['maxreads'], dsnps)
    if args['samtools'] is not None and len(args['samtools']) > 0:
        cmd += " --samtools {}".format(args['samtools'])
    print(cmd)
    runcmd(cmd, dsnps, log="snps.log", rundir=args['rundir'])

    log('Computing BAFs\n', level='PROGRESS')
    cmd = 'python3 -m hatchet deBAF -N {} -T {} -S {} -r {} -j {} -q {} -Q {} -U {} -c {} -C {} -O {} -o {} -L {}'
    nbaf = os.path.join(dbaf, 'normal.1bed')
    tbaf = os.path.join(dbaf, 'bulk.1bed')
    vcfs = [os.path.join(dsnps, f) for f in os.listdir(dsnps) if f.endswith('.vcf.gz')]
    cmd = cmd.format(args['normal'], ' '.join(args['tumor']), 'normal ' + ' '.join(args['names']), 
                     args['ref'], args['J'], args['phred'], args['phred'], args['phred'], 
                     args['minreads'], args['maxreads'], nbaf, tbaf, " ".join(vcfs))
    if args['samtools'] is not None and len(args['samtools']) > 0:
        cmd += " --samtools {}".format(args['samtools'])
    if args['bcftools'] is not None and len(args['bcftools']) > 0:
        cmd += " --bcftools {}".format(args['bcftools'])
    print(cmd)
    runcmd(cmd, dbaf, log="bafs.log", rundir=args['rundir'])

    log('Combining RDRs and BAFs\n', level='PROGRESS')
    ctot = os.path.join(args['rundir'], config.bin.outputtotal)
    cmd = 'python3 -m hatchet comBBo -c {} -C {} -B {} -t {}'
    cmd = cmd.format(nbin, tbin, tbaf, ctot)
    if args['seed'] is not None:
        cmd += " -e {}".format(args['seed'])
    print(cmd)
    runcmd(cmd, dbb, out='bulk.bb', log="combo.log", rundir=args['rundir'])

    log('Counting reads at each position\n', level='PROGRESS')
    cmd = 'python3 -m hatchet countPos -N {} -T {} -O {} -j {}'
    cmd = cmd.format(args['normal'], ' '.join(args['tumor']), dcounts, args['J'])
    if args['samtools'] is not None and len(args['samtools']) > 0:
        cmd += " --samtools {}".format(args['samtools'])
    print(cmd)
    runcmd(cmd, dcounts, log="counts.log", rundir=args['rundir'])

    log('Preparing gzip file for transfer\n', level='PROGRESS')


def generate_subdirectories(args):
    dbaf = os.path.join(args['rundir'], 'baf')
    if not os.path.isdir(dbaf):
        os.mkdir(dbaf)

    drdr = os.path.join(args['rundir'], 'rdr')
    if not os.path.isdir(drdr):
        os.mkdir(drdr)

    dbb = os.path.join(args['rundir'], 'bb')
    if not os.path.isdir(dbb):
        os.mkdir(dbb)
        
    dsnps = os.path.join(args['rundir'], 'snps')
    if not os.path.isdir(dsnps):
        os.mkdir(dsnps)

    dcounts = os.path.join(args['rundir'], 'counts')
    if not os.path.isdir(dcounts):
        os.mkdir(dcounts)

    drisk = os.path.join(args['rundir'], 'risk')
    if not os.path.isdir(drisk):
        os.mkdir(drisk)

    return dbaf, drdr, dbb, dsnps, dcounts, drisk

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
