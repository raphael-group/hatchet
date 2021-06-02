#!/usr/bin/env python3

import sys, os
import argparse
import subprocess as sp
import multiprocessing as mp
import shlex
import re
import pathlib
import shutil

from .Supporting import *
from hatchet import config
from.ArgParsing import extractChromosomes, parse_preprocess_args

def main(args=None):
    log('Parsing and checking arguments\n', level='PROGRESS')
    args = parse_preprocess_args(args)
    log('\n'.join(['Arguments:'] + ['\t{} : {}'.format(a, args[a]) for a in args]) + '\n', level='INFO')

    log('Setting directories\n', level='PROGRESS')
    dbaf, drdr, dbb, dsnps, dcounts, dabin = generate_subdirectories(args)
    
    samtools = os.path.join(args["samtools"], "samtools")
    all_names = ['normal'] + args['names']
    all_samples =  [args['normal']] + args['tumor']    
    chrs = extractChromosomes(samtools,  [args["normal"], "normal"], [(x, "") for x in args["tumor"]])


    ctot = os.path.join(args['rundir'], config.bin.outputtotal)
    nbaf = os.path.join(dbaf, 'normal.1bed')
    tbaf = os.path.join(dbaf, 'bulk.1bed')
    nbin = os.path.join(drdr, 'normal.1bed')
    tbin = os.path.join(drdr, 'bulk.1bed')
    
    if len(missing_snps(dsnps, chrs)) > 0:
        log('Calling SNPs\n', level='PROGRESS')
        cmd =  'python3 -m hatchet SNPCaller -N {} -r {} -j {} -c {} -C {} -o {}'
        # TODO: expose -R reference SNPs list?
        cmd = cmd.format(args['normal'], args['ref'], args['J'], args['minreads'], args['maxreads'], dsnps)
        if args['samtools'] is not None and len(args['samtools']) > 0:
            cmd += " --samtools {}".format(args['samtools'])
        runcmd(cmd, dsnps, log="snps.log", rundir=args['rundir'])
    else:
        log('Found all SNPCaller output, skipping step.\n', level='PROGRESS')
     
    if len(missing_baf(dbaf)) > 0:
        log('Computing BAFs\n', level='PROGRESS')
        cmd = 'python3 -m hatchet deBAF -N {} -T {} -S {} -r {} -j {} -q {} -Q {} -U {} -c {} -C {} -O {} -o {} -L {} -l {}'
        vcfs = [os.path.join(dsnps, f) for f in os.listdir(dsnps) if f.endswith('.vcf.gz')]
        cmd = cmd.format(args['normal'], ' '.join(args['tumor']), 'normal ' + ' '.join(args['names']), 
                        args['ref'], args['J'], args['phred'], args['phred'], args['phred'], 
                        args['minreads'], args['maxreads'], nbaf, tbaf, " ".join(vcfs), dbaf)
        if args['samtools'] is not None and len(args['samtools']) > 0:
            cmd += " --samtools {}".format(args['samtools'])
        if args['bcftools'] is not None and len(args['bcftools']) > 0:
            cmd += " --bcftools {}".format(args['bcftools'])
        runcmd(cmd, dbaf, log="bafs.log", rundir=args['rundir'])
    else:
        log('Found all deBAF output, skipping step.\n', level='PROGRESS')
    
    if len(missing_rdr(drdr, ctot)) > 0:
        log('Computing RDRs\n', level='PROGRESS')
        cmd = 'python3 -m hatchet binBAM -N {} -T {} -S {} -b {} -g {} -j {} -q {} -O {} -o {}'
        cmd = cmd.format(args['normal'], ' '.join(args['tumor']), 'normal ' + ' '.join(args['names']), args['size'], args['ref'], args['J'], args['phred'], nbin, tbin)
        if args['samtools'] is not None and len(args['samtools']) > 0:
            cmd += " --samtools {}".format(args['samtools'])
        runcmd(cmd, drdr, log="bins.log", rundir=args['rundir'])
    else:
        log('Found all bimBAM output, skipping step.\n', level='PROGRESS')
       
    if len(missing_bb(dbb)) > 0:
        log('Combining RDRs and BAFs\n', level='PROGRESS')
        cmd = 'python3 -m hatchet comBBo -c {} -C {} -B {} -t {}'
        cmd = cmd.format(nbin, tbin, tbaf, ctot)
        if args['seed'] is not None:
            cmd += " -e {}".format(args['seed'])
        runcmd(cmd, dbb, out='bulk.bb', log="combo.log", rundir=args['rundir'])
    else:
        log('Found all comBBo output, skipping step.\n', level='PROGRESS')
    
    if len(missing_counts(dcounts, chrs, all_names)) > 0:
        log('Counting reads at each position\n', level='PROGRESS')
        cmd = 'python3 -m hatchet countPos -B {} -O {} -j {} -S {}'
        cmd = cmd.format(' '.join([args['normal']] + args['tumor']), dcounts, args['J'], ' '.join(all_names))
        if args['samtools'] is not None and len(args['samtools']) > 0:
            cmd += " --samtools {}".format(args['samtools'])
        runcmd(cmd, dcounts, log="counts.log", rundir=args['rundir'])
    else:
        log('Found all countPos output, skipping step.\n', level='PROGRESS')

    if len(missing_arrays(dabin, chrs)) > 0:
        log('Forming intermediate arrays from count files\n', level='PROGRESS')
        log('USING CENTROMERES FOR VERSION hg38 REFERENCE GENOME\n', level='WARN')

        # NOTE: replace with a different centromeres file if reference genome is not HG38
        centromeres_file = os.path.join(pathlib.Path(__file__).parent.parent.absolute(), 'resources', 'hg38.centromeres.txt')
        if not os.path.exists(centromeres_file):
            raise ValueError(error(f"Could not find centromeres files {centromeres_file}!"))
        cmd = f"python3 -m hatchet array -s {args['rundir']} -o {os.path.join(dabin, 'intermediate')} -j {args['J']} -C {centromeres_file}"
        runcmd(cmd, dabin, log="formArray.log", rundir=args['rundir'])
        
        if len(missing_arrays(dabin, chrs)) > 0:
            log('formArray step failed, time to debug :(\n', level='PROGRESS')
        else:
            log('formArray step succeeded, removing counts files\n', level='PROGRESS')
            shutil.rmtree(dcounts)

    else:
        log('Found all formArray output, skipping step.\n', level='PROGRESS')

    log('Checking for output files\n', level='PROGRESS')
    missing = []
    missing.extend(missing_snps(dsnps, chrs))
    missing.extend(missing_baf(dbaf))
    missing.extend(missing_rdr(drdr, ctot))
    missing.extend(missing_bb(dbb))
    missing.extend(missing_arrays(dabin, chrs))
    
    with open(os.path.join(args['rundir'], 'missing_files.log'), 'w') as f:
        if len(missing) == 0:
            log("No output files missing.\n", level="INFO")
            # leave log file empty
        else:
            log("Found missing output files (see missing_files.log).\n", level="INFO")
            f.write('\n'.join(missing))
            f.write('\n')

    if args['zip']:
        log('Collecting and compressing output files\n', level='PROGRESS')
        cmd = f"tar -czvf {args['output'] + '.tar.gz'} {dbaf} {drdr} {dbb} {dsnps} {ctot} {dabin}"
        sp.run(cmd.split())
        
    log('Done\n', level='PROGRESS')

def missing_snps(dsnps, chrs):
    missing = []
    for chr in chrs:
        fname = os.path.join(dsnps, chr + '.vcf.gz')
        if not os.path.exists(fname):
            missing.append("SNPS: Missing file {}".format(fname))
    return missing

def missing_baf(dbaf):
    # deBAF output (baf)
    missing = []
    beds = ['bulk.1bed', 'normal.1bed']
    for file in beds:
        fname = os.path.join(dbaf, file)
        if not os.path.exists(fname):
            missing.append("BAF: Missing file {}".format(fname))
    return missing

def missing_rdr(drdr, ctot):
    # binBAM output (rdr)
    missing = []
    beds = ['bulk.1bed', 'normal.1bed']
    for file in beds:
        fname = os.path.join(drdr, file)
        if not os.path.exists(fname):
            missing.append("RDR: Missing file {}".format(fname))
    if not os.path.exists(ctot):
        missing.append("RDR: Missing file {}".format(ctot))
    return missing

def missing_bb(dbb):
    missing = []
    # comBBo (bb)
    fname = os.path.join(dbb, 'bulk.bb')
    if not os.path.exists(fname):
        missing.append("BB: Missing file {}".format(fname))
    return missing
  
def missing_counts(dcounts, chrs, all_names):
    missing = []
    ### position counting (counts)
    # samtools output
    counts_files = [a for a in os.listdir(dcounts) if a.endswith('starts.gz')]
    if len(counts_files) <= 1:
        missing.append("COUNTS: Missing all count files from directory {}".format(dcounts))
    else:
        for chr in chrs:
            for name in all_names:
                fname = os.path.join(dcounts, '.'.join([name, chr, 'starts.gz']))
                if not os.path.exists(fname):
                    missing.append("COUNTS: Missing samtools output file {}".format(fname))
        # mosdepth output
        for name in all_names:
            fname = os.path.join(dcounts, name + '.mosdepth.global.dist.txt')
            if not os.path.exists(fname):
                missing.append("COUNTS: Missing mosdepth output file {}".format(fname))   
            fname = os.path.join(dcounts, name + '.mosdepth.summary.txt')
            if not os.path.exists(fname):
                missing.append("COUNTS: Missing mosdepth output file {}".format(fname))   
            fname = os.path.join(dcounts, name + '.per-base.bed.gz')
            if not os.path.exists(fname):
                missing.append("COUNTS: Missing mosdepth output file {}".format(fname))   
            fname = os.path.join(dcounts, name + '.per-base.bed.gz.csi')
            if not os.path.exists(fname):
                missing.append("COUNTS: Missing mosdepth output file {}".format(fname))   
    return missing

def missing_arrays(dabin, chrs):
    missing = []
    # formArray (abin/intermediate.* files))

    fname = os.path.join(dabin, f'intermediate.samples')
    if not os.path.exists(fname):
        missing.append(fname)     
    
    for ch in chrs:
        fname = os.path.join(dabin, f'intermediate.{ch}.total')
        if not os.path.exists(fname):
            missing.append(fname)
        fname = os.path.join(dabin, f'intermediate.{ch}.thresholds')
        if not os.path.exists(fname):
            missing.append(fname)
    
    return missing

def pileup_wrapper(params):
    cmd, name = params
    log("Counting sample \"{}\"\n".format(name), level="INFO")
    worker = sp.run(cmd.split())
    worker.check_returncode()
    log("Done counting sample \"{}\"\n".format(name), level="INFO")

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

    dabin = os.path.join(args['rundir'], 'abin')
    if not os.path.isdir(dabin):
        os.mkdir(dabin)

    return dbaf, drdr, dbb, dsnps, dcounts, dabin

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
