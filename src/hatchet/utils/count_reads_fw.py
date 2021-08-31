#!/usr/bin/python3

import sys
import os.path
import argparse

from . import BAMBinning as bb
from . import TotalCounting as tc
from . import ArgParsing as ap
from .Supporting import *



def main(args=None):
    log(msg="# Parsing and checking input arguments\n", level="STEP")
    args = ap.parse_count_reads_fw_arguments(args)
    logArgs(args, 80)

    if args["regions"] is None:
        log(msg="# Retrieving genomic regions to consider from maximum chromosome length\n", level="STEP")
        regions = knownRegions(args["refdict"], args["chromosomes"])
    else:
        log(msg="# Checking the consistency of the given regions\n", level="STEP")
        regions = ap.parseRegions(args["regions"], args["chromosomes"])

    if args["verbose"]:
        msg = "regions: "
        for c in args["chromosomes"]:
            msg += " {}: {}".format(c, regions[c])
        msg += "\n"
        log(msg=msg, level="INFO")

    log(msg="# Binning and counting the normal sample\n", level="STEP")
    normal_bins = bb.bin(samtools=args["samtools"], samples=[args["normal"]], chromosomes=args["chromosomes"],
                         num_workers=args["j"], q=args["q"], size=args["size"], regions=regions, verbose=args["verbose"])
    if not normal_bins: close("No bins in the normal sample!\n")

    log(msg="# Writing the read counts for bins of normal sample\n", level="STEP")
    handle = open(args['outputNormal'], 'w') if args['outputNormal'] is not None else sys.stdout
    for c in args["chromosomes"]:
        for (samplename, chromosome, start, stop, n_reads) in normal_bins[args["normal"][1], c]:
            handle.write('\t'.join([str(chromosome), str(start), str(stop), samplename, str(n_reads)]) + '\n')
    if handle is not sys.stdout:
        handle.close()

    log(msg="# Binning and counting the tumor samples\n", level="STEP")
    tumor_bins = bb.bin(samtools=args["samtools"], samples=args["samples"], chromosomes=args["chromosomes"],
                         num_workers=args["j"], q=args["q"], size=args["size"], regions=regions, verbose=args["verbose"])
    if not tumor_bins: close("No bins in the tumor samples!\n")

    log(msg="# Writing the read counts for bins of tumor samples\n", level="STEP")
    handle = open(args['outputTumors'], 'w') if args['outputTumors'] is not None else sys.stdout
    for sample in sorted(args["samples"]):
        for c in args["chromosomes"]:
            for (samplename, chromosome, start, stop, n_reads) in tumor_bins[sample[1], c]:
                handle.write('\t'.join([str(chromosome), str(start), str(stop), samplename, str(n_reads)]) + '\n')
    if handle is not sys.stdout:
        handle.close()

    log(msg="# Counting total number of reads for normal and tumor samples\n", level="STEP")
    total_counts = tc.tcount(samtools=args["samtools"], samples=({args["normal"]}|args["samples"]), chromosomes=args["chromosomes"],
                             num_workers=args["j"], q=args["q"], verbose=args["verbose"])

    try:
        total = {sample[1] : sum(total_counts[sample[1], chromosome] for chromosome in args["chromosomes"]) for sample in args["samples"]}
        total[args["normal"][1]] = sum(total_counts[args["normal"][1], chromosome] for chromosome in args["chromosomes"])
    except:
        raise KeyError(error("Either a chromosome or a sample has not been considered in the total counting!"))

    log(msg="# Writing the total read counts for all samples in {}\n".format(args["outputTotal"]), level="STEP")
    with open(args["outputTotal"], 'w') as f:
        f.write("{}\t{}\n".format(args["normal"][1], total[args["normal"][1]]))
        for sample in sorted(args["samples"]):
            f.write("{}\t{}\n".format(sample[1], total[sample[1]]))


def knownRegions(refdict, chromosomes):
    ends = {c : None for c in chromosomes}
    assert os.path.isfile(refdict)
    with open(refdict, 'r') as i:
        for l in i:
            if '@SQ' in l:
                assert 'SN:' in l and 'LN:' in l
                c = l.split('SN:')[1].split()[0]
                if c in chromosomes:
                    end = int(l.split('LN:')[1].split()[0])
                    ends[c] = end
    if None in ends.values():
        log(msg="The following chromosomes have not been found in the dictionary of the reference genome: \n\t{}".format(','.join([c for c in ends if ends[c] == None])), level="WARN")

    res = {}
    for c in chromosomes:
        if ends[c] is None:
            raise ValueError("Length of chromosome {} could not be determined. Are you using the correct reference genome?".format(c))
        res[c] = [(0, ends[c])]
        
    return res


def logArgs(args, width):
    text = "\n"
    for key in args:
        text += "\t{}: {}\n".format(key, args[key])
    log(msg=text, level="INFO")


if __name__ == '__main__':
    main()
