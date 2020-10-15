import os.path
import argparse
import TotalCounting as tc
import ArgParsing as ap
from Supporting import *
import Supporting as sp
from hatchet import config


def parse_arguments(args=None):
    description = ""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-b", "--bams", required=True, type=str, nargs='+',
                        help="BAM files to process")
    parser.add_argument("-st", "--samtools", required=False, default=config.paths.samtools, type=str,
                        help="Path to the directory to \"samtools\" executable, required in default mode (default: samtools is directly called as it is in user $PATH)")
    parser.add_argument("-r","--regions", required=False, default=None, type=str, help="BED file containing the a list of genomic regions to consider in the format \"CHR  START  END\", REQUIRED for WES data (default: none, consider entire genome)")
    parser.add_argument("-j", "--processes", required=False, default=2, type=int,
                        help="Number of available parallel processes (default: 2)")
    parser.add_argument("-q", "--readquality", required=False, default=0, type=int,
                        help="Minimum mapping quality for an aligned read to be considered (default: 0)")
    parser.add_argument("-t", "--outputtotal", required=False, default="total_read.counts", type=str,
                        help="Output filename for total read counts in all tumor samples (default: \"total_read.counts\")")
    parser.add_argument("-v", "--verbose", action='store_true', default=False, required=False,
                        help="Use verbose log messages")
    args = parser.parse_args(args)

    # Check the region file
    if args.regions is not None and not os.path.isfile(args.regions):
        raise ValueError(sp.error("The specified region file does not exist"))

    # In default mode, check the existence and compatibility of samtools and bcftools
    samtools = os.path.join(args.samtools, "samtools")
    if sp.which(samtools) is None:
        raise ValueError(sp.error("samtools has not been found or is not executable!"))

    bams = [(f, 'bam' + str(i)) for i, f in enumerate(args.bams)]
    chromosomes = ['chr' + str(i) for i in range(1, 23)]

    if not args.processes > 0: raise ValueError(sp.error("The number of parallel processes must be greater than 0"))
    if not args.readquality >= 0: raise ValueError(sp.error("The read mapping quality must be positive"))

    return {"bams": bams,
            "chromosomes": chromosomes,
            "samtools": samtools,
            "regions": args.regions,
            "j": args.processes,
            "q": args.readquality,
            "outputTotal": args.outputtotal,
            "verbose": args.verbose}


def main(args=None):
    log(msg="# Parsing and checking input arguments\n", level="STEP")
    args = parse_arguments(args)

    with open(args["outputTotal"], 'w') as f:
        f.write(str(args) + '\n\n')

        for bamfile, _ in args['bams']:
            folder = os.path.dirname(bamfile)
            for file in os.listdir(folder):
                f.write(file + '\n')

        f.write('\n\n')

        total_counts = tc.tcount(samtools=args["samtools"], samples=args["bams"],
                                 chromosomes=args["chromosomes"],
                                 num_workers=args["j"], q=args["q"], verbose=args["verbose"])

        total = {bam[1]: sum(total_counts[bam[1], chromosome] for chromosome in args["chromosomes"]) for bam in
                 args["bams"]}

        log(msg="# Writing the total read counts for all samples in {}\n".format(args["outputTotal"]), level="STEP")
        for bam in args["bams"]:
            f.write("{}\t{}\n".format(bam[1], total[bam[1]]))


if __name__ == '__main__':
    main()
