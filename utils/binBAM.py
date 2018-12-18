#!/usr/bin/python2

import sys
import os.path
import argparse

import BAMBinning as bb
import TotalCounting as tc
import ArgParsing as ap
from Supporting import *



def main():
    log(msg="# Parsing and checking input arguments\n", level="STEP")
    args = ap.parse_bin_arguments()
    logArgs(args, 80)

    if args["regions"] is None:
        log(msg="# Retrieving genomic regions to consider from maximum chromosome length\n", level="STEP")
        regions = knownRegions(args["referencename"], args["chromosomes"])
    else:
        log(msg="# Checking the consistency of the given regions\n", level="STEP")
        regions = parseRegions(args["regions"], args["chromosomes"])

    if args["verbose"]:
        msg = " Regions:\n"
        for c in args["chromosomes"]:
            msg += "\t{}: {}\n".format(c, regions[c])
    else:
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
    if args["outputNormal"] is not None:
        with open(args["outputNormal"], 'w') as f:
            for c in args["chromosomes"]:
                for count in normal_bins[args["normal"][1], c]:
                    f.write("{}\t{}\t{}\t{}\t{}\n".format(count[0], count[1], count[2], count[3], count[4]))
    else:
        for c in args["chromosomes"]:
            for count in normal_bins[args["normal"][1], c]:
                sys.stdout.write("{}\t{}\t{}\t{}\t{}\n".format(count[0], count[1], count[2], count[3], count[4]))

    log(msg="# Binning and counting the tumor samples\n", level="STEP")
    tumor_bins = bb.bin(samtools=args["samtools"], samples=args["samples"], chromosomes=args["chromosomes"],
                         num_workers=args["j"], q=args["q"], size=args["size"], regions=regions, verbose=args["verbose"])
    if not tumor_bins: close("No bins in the tumor samples!\n")

    log(msg="# Writing the read counts for bins of tumor samples\n", level="STEP")
    if args["outputTumors"] is not None:
        with open(args["outputTumors"], 'w') as f:
            for sample in args["samples"]:
                for c in args["chromosomes"]:
                    for count in tumor_bins[sample[1], c]:
                        f.write("{}\t{}\t{}\t{}\t{}\n".format(count[0], count[1], count[2], count[3], count[4]))
    else:
        for sample in args["samples"]:
            for c in args["chromosomes"]:
                for count in tumor_bins[sample[1], c]:
                    sys.stdout.write("{}\t{}\t{}\t{}\t{}\n".format(count[0], count[1], count[2], count[3], count[4]))

    log(msg="# Counting total number of reads for normal and tumor samples\n", level="STEP")
    total_counts = tc.tcount(samtools=args["samtools"], samples=({args["normal"]}|args["samples"]), chromosomes=args["chromosomes"],
                             num_workers=args["j"], q=args["q"], verbose=args["verbose"])

    try:
        total = {sample[1] : sum(total_counts[sample[1], chromosome] for chromosome in args["chromosomes"]) for sample in args["samples"]}
        total[args["normal"][1]] = sum(total_counts[args["normal"][1], chromosome] for chromosome in args["chromosomes"])
    except:
        raise KyeError(sp.error("Either a chromosome or a sample has not been considered in the total counting!"))

    log(msg="# Writing the total read counts for all samples in {}\n".format(args["outputTotal"]), level="STEP")
    with open(args["outputTotal"], 'w') as f:
        f.write("{}\t{}\n".format(args["normal"][1], total[args["normal"][1]]))
        for sample in args["samples"]:
            f.write("{}\t{}\n".format(sample[1], total[sample[1]]))


def knownRegions(referencename, chromosomes):
    ends = []
    if referencename == "hg18":
        ends = [-1, 247249719, 242951149, 199501827, 191273063, 180857866, 170899992, 158821424, 146274826, 140273252, 135374737, 134452384, 132349534, 114142980, 106368585, 100338915, 88827254, 78774742, 76117153, 63811651, 62435964, 46944323, 49691432]
    elif referencename == "hg19":
        ends = [-1, 249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566]
    elif referencename == "hg38":
        ends = [-1, 248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636, 138394717, 133797422, 135086622, 133275309, 114364328, 107043718, 101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983, 50818468]
    else:
        raise ValueError(error("The given reference name is not recognized!"))
    res = {}
    for c in chromosomes:
        res[c] = [(0, ends[int(digits(c))])]
    return res



def logArgs(args, width):
    text = "\n"
    for key in args:
        text += "\t{}: {}\n".format(key, args[key])
    log(msg=text, level="INFO")


if __name__ == '__main__':
    main()
