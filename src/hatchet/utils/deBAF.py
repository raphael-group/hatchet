#!/usr/bin/python2

import sys
import os.path
import argparse
from multiprocessing import Process, Queue, JoinableQueue, Lock, Value
from scipy.stats import beta
# from statsmodels.stats.proportion import *


import SNPCalling
import AlleleCounting
import ArgParsing as ap
from Supporting import *
import Supporting as sp


def main():
    log(msg="# Parsing the input arguments, checking the consistency of given files, and extracting required information\n", level="STEP")
    args = ap.parse_baf_arguments()
    logArgs(args, 80)

    if args["reference"] is not None:
        defaultMode(args)
    else:
        naiveMode(args)


def defaultMode(args):
    log(msg="# Inferring SNPs from the normal sample\n", level="STEP")
    snps = SNPCalling.call(samtools=args["samtools"], bcftools=args["bcftools"], reference=args["reference"], samples=[args["normal"]],
                           chromosomes=args["chromosomes"], num_workers=args["j"], snplist=args["snps"], q=args["q"], Q=args["Q"],
                           qual=args["qual"], mincov=args["mincov"], dp=args["maxcov"], E=args["E"], regions=args["regions"], verbose=args["verbose"])
    print "BLAGG"
    if not snps: sp.close("No SNPs found in the normal!\n")

    log(msg="# Selecting heterozygous SNPs\n", level="STEP")
    hetSNPs = selectHetSNPs(counts=snps, gamma=args["gamma"], maxshift=args["maxshift"], verbose=args["verbose"])
    if not hetSNPs: sp.close("No heterozygous SNPs found in the selected regions of the normal!\n")

    log(msg="# Writing the list of selected SNPs, covered and heterozygous in the normal sample\n", level="STEP")
    with open(args["outputSnps"], 'w') as f:
        for chro in args["chromosomes"]:
            if (args["normal"][1], chro) in hetSNPs:
                for snp in hetSNPs[args["normal"][1], chro]: f.write("{}\t{}\n".format(snp[1], snp[2]))

    log(msg="# Writing the allele counts of the normal sample for selected SNPs\n", level="STEP")
    if args["outputNormal"] is not None:
        with open(args["outputNormal"], 'w') as f:
            for chro in args["chromosomes"]:
                if (args["normal"][1], chro) in hetSNPs:
                    for count in hetSNPs[args["normal"][1], chro]: f.write("{}\t{}\t{}\t{}\t{}\n".format(count[0], count[1], count[2], count[3], count[4]))
    else:
        for chro in args["chromosomes"]:
            if (args["normal"][1], chro) in hetSNPs:
                for count in hetSNPs[args["normal"][1], chro]: sys.stdout.write("{}\t{}\t{}\t{}\t{}\n".format(count[0], count[1], count[2], count[3], count[4]))

    log(msg="# Counting the alleles of tumor samples for selected SNPs\n", level="STEP")
    counts = AlleleCounting.count(samtools=args["samtools"], bcftools=args["bcftools"], reference=args["reference"], samples=args["samples"],
                                  chromosomes=args["chromosomes"], num_workers=args["j"], snplist=args["outputSnps"], q=args["q"], Q=args["Q"],
                                  mincov=args["mincov"], dp=args["maxcov"], E=args["E"], verbose=args["verbose"])
    if not counts: sp.close("The selected SNPs are not covered in the tumors!\n")

    log(msg="# Writing the allele counts of tumor samples for selected SNPs\n", level="STEP")
    if args["outputTumors"] is not None:
        with open(args["outputTumors"], 'w') as f:
            for sample in args["samples"]:
                for chro in args["chromosomes"]:
                    if (sample[1], chro) in counts:
                        for count in counts[sample[1], chro]: f.write("{}\t{}\t{}\t{}\t{}\n".format(count[0], count[1], count[2], count[3], count[4]))
    else:
        for sample in args["samples"]:
            for chro in args["chromosomes"]:
                if (sample[1], chro) in counts:
                    for count in counts[sample[1], chro]: sys.stdout.write("{}\t{}\t{}\t{}\t{}\n".format(count[0], count[1], count[2], count[3], count[4]))


def naiveMode(args):
    log(msg="# Inferring SNPs from the normal sample\n", level="STEP")
    snps = AlleleCounting.naiveCount(samtools=args["samtools"], samples=[args["normal"]], chromosomes=args["chromosomes"], num_workers=args["j"],
                                     snplist=args["snps"], q=args["q"], Q=args["Q"], E=args["E"], regions=args["regions"], verbose=args["verbose"])
    if not snps: sp.close("No SNPs found in the normal!\n")

    log("# Selecting heterozygous SNPs\n", level="STEP")
    hetSNPs = selectNaiveHetSNPs(counts=snps, chromosomes=args["chromosomes"], mincov=args["mincov"], maxcov=args["maxcov"],
                                 gamma=args["gamma"], maxshift=args["maxshift"], verbose=args["verbose"])
    if not hetSNPs: sp.close("No heterozygous SNPs found in the selected regions of the normal!\n")

    log("# Writing the list of selected SNPs, covered and heterozygous in the normal sample\n", level="STEP")
    with open(args["outputSnps"], 'w') as f:
        for chro in args["chromosomes"]:
            if (args["normal"][1], chro) in hetSNPs:
                for snp in hetSNPs[args["normal"][1], chro]: f.write("{}\t{}\n".format(snp[1], snp[2]))

    log("# Writing the allele counts of the normal sample for selected SNPs\n", level="STEP")
    if args["outputNormal"] is not None:
        with open(args["outputNormal"], 'w') as f:
            for chro in args["chromosomes"]:
                if (args["normal"][1], chro) in hetSNPs:
                    for count in hetSNPs[args["normal"][1], chro]:
                        f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(count[0], count[1], count[2], count[4], count[5], count[6], count[7]))
    else:
        for chro in args["chromosomes"]:
            if (args["normal"][1], chro) in hetSNPs:
                for count in hetSNPs[agrs["normal"][1], chro]:
                    sys.stdout.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(count[0], count[1], count[2], count[4], count[5], count[6], count[7]))

    log("# Counting the alleles of tumor samples for selected SNPs\n", level="STEP")
    snps = AlleleCounting.naiveCount(samtools=args["samtools"], samples=args["samples"], chromosomes=args["chromosomes"], num_workers=args["j"],
                                     snplist=args["outputSnps"], q=args["q"], Q=args["Q"], E=args["E"], verbose=args["verbose"])
    if not counts: sp.close("The selected SNPs are not covered in the tumors!\n")

    log("# Writing the allele counts of tumor samples for selected SNPs\n", level="STEP")
    if args["outputTumors"] is not None:
        with open(args["outputTumors"], 'w') as f:
            for sample in args["samples"]:
                for chro in args["chromosomes"]:
                    if (sample[1], chro) in snps:
                        for count in snps[sample[1], chro]:
                            f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(count[0], count[1], count[2], count[4], count[5], count[6], count[7]))
    else:
        for sample in args["samples"]:
            for chro in args["chromosomes"]:
                for count in hetSNPs[snps, chro]:
                    if (sample[1], chro) in snps:
                        sys.stdout.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(count[0], count[1], count[2], count[4], count[5], count[6], count[7]))



def selectHetSNPs(counts, gamma, maxshift, verbose):
    hetSNPs = {}
    for key in counts:
        hetSNPs[key] = [count for count in counts[key] if isHet(count[3], count[4], gamma) and checkShift(count[3], count[4], maxshift)]
    return hetSNPs


def selectNaiveHetSNPs(counts, chromosomes, mincov, maxcov, gamma, maxshift, verbose):
    hetSNPs = {}
    for key in counts:
        hetSNPs[key] = []
        append = hetSNPs[key].append
        for count in counts[key]:
            if mincov <= count[3] <= maxcov:
                alleles = sorted([count[4], count[5], count[6], count[7]], reverse=True)
                if sum(x > 0 for x in alleles) == 2 and isHet(alleles[0], alleles[1], gamma) and checkShift(alleles[0], alleles[1], maxshift):
                    append(count)
    return hetSNPs


def filterSNPs(snps, chromosomes, regions):
    reg = ap.parseRegions(regions, chromosomes)
    res = {}
    for key in snps:
        res[key] = [count for count in snps[key] if sum(l <= count[2] <= r for (l, r) in reg[count[1]])]
    return res


def isHet(countA, countB, gamma):
    p_lower = gamma / 2.0
    p_upper = 1.0 - p_lower
    [c_lower, c_upper] = beta.ppf([p_lower, p_upper], countA + 1, countB + 1)
    return c_lower <= 0.5 <= c_upper

# def isHet(countA, countB, gamma):
#     lb, ub = proportion_confint(min(countA, countB), countA+countB, alpha=gamma, method='beta')
#     return lb <= 0.5 <= ub

# def isHet(countA, countB, gamma):
#     lb, ub = proportion_confint(min(countA, countB), countA+countB, alpha=gamma, method='jeffreys')
#     return lb <= 0.5 <= ub


def checkShift(countA, countB, maxshift):
    return (0.5 - (float(min(countA, countB)) / float(countA + countB)) ) <= maxshift


def logArgs(args, width):
    text = "\n"
    for key in args:
        text += "\t{}: {}\n".format(key, args[key])
    log(msg=text, level="INFO")


if __name__ == '__main__':
    main()
