import sys
import os.path
import argparse
import subprocess
import shlex

import Supporting as sp



def parse_baf_arguments(args=None):
    """
    Parse command line arguments
    Returns:
    """
    description = "Count the A/B alleles from a matched-normal BAM file and multiple tumor BAM files in specified SNP positions or estimated heterozygous SNPs in the normal genome. This tool can be applied both to whole-genome sequencing (WGS) data or whole-exome sequencing (WES) data, but coding regions must be specified as a BED file in the case of WES."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-N","--normal", required=True, type=str, help="BAM file corresponding to matched normal sample")
    parser.add_argument("-T","--tumors", required=True, type=str, nargs='+', help="BAM files corresponding to samples from the same tumor")
    parser.add_argument("-r","--reference", required=True, type=str, help="Human-genome reference corresponding to given samples")
    parser.add_argument("-S","--samples", required=False, default=None, type=str, nargs='+', help="Sample names for each BAM (given in the same order where the normal name is first)")
    parser.add_argument("-st","--samtools", required=False, default="", type=str, help="Path to the directory to \"samtools\" executable, required in default mode (default: samtools is directly called as it is in user $PATH)")
    parser.add_argument("-bt","--bcftools", required=False, default="", type=str, help="Path to the directory of \"bcftools\" executable, required in default mode (default: bcftools is directly called as it is in user $PATH)")
    parser.add_argument("-L","--snps", required=False, default=None, type=str, help="List of SNPs to consider in the normal sample (default: heterozygous SNPs are inferred from the normal sample)")
    parser.add_argument("-e","--regions", required=False, default=None, type=str, help="BED file containing the a list of genomic regions to consider in the format \"CHR  START  END\", REQUIRED for WES data with coding regions (default: none, consider entire genome)")
    parser.add_argument("-j", "--processes", required=False, default=2, type=int, help="Number of available parallel processes (default: 2)")
    parser.add_argument("-q", "--readquality", required=False, default=0, type=int, help="Minimum mapping quality for an aligned read to be considered (default: 0)")
    parser.add_argument("-Q", "--basequality", required=False, default=13, type=int, help="Minimum base quality for a base to be considered (default: 13)")
    parser.add_argument("-U", "--snpquality", required=False, default=20, type=int, help="Minimum SNP-variant quality, QUAL, for a variant to be considered (default: 20)")
    parser.add_argument("-g", "--gamma", required=False, default=0.05, type=float, help="Level of confidence to determine heterozigosity of SNPs (default: 0.05)")
    parser.add_argument("-b", "--maxshift", required=False, default=0.5, type=float, help="Maximum allowed absolute difference of BAF from 0.5 for selected heterozygous SNPs in the normal sample (default: 0.5)")
    parser.add_argument("-c", "--mincov", required=False, default=0, type=int, help="Minimum coverage for SNPs to be considered (default: 0)")
    parser.add_argument("-C", "--maxcov", required=False, default=1000, type=int, help="Maximum coverage for SNPs to be considered (default: 1000, suggested: twice the values of expected average coverage to avoid aligning artefacts)")
    parser.add_argument("-E","--newbaq", required=False, action='store_true', default=False, help="Recompute alignment of reads on the fly during SNP calling (default: false)")
    parser.add_argument("-O", "--outputnormal", required=False, default=None, type=str, help="Filename of output for allele counts in the normal sample (default: standard output)")
    parser.add_argument("-o", "--outputtumors", required=False, default=None, type=str, help="Output filename for allele counts in tumor samples (default: standard output)")
    parser.add_argument("-l", "--outputsnps", required=False, default="selectedSNPs.csv", type=str, help="Output filename for list of selected SNPs (default: selectedSNPs.txt)")
    parser.add_argument("-v", "--verbose", action='store_true', default=False, required=False, help="Use verbose log messages")
    args = parser.parse_args(args)

    # Parse BAM files, check their existence, and infer or parse the corresponding sample names
    normalbaf = args.normal
    if not os.path.isfile(normalbaf): raise ValueError(sp.error("The specified normal BAM file does not exist"))
    tumors = args.tumors
    for tumor in tumors:
        if(not os.path.isfile(tumor)): raise ValueError(sp.error("The specified normal BAM file does not exist"))
    names = args.samples
    if names != None and (len(tumors)+1) != len(names):
        raise ValueError(sp.error("A sample name must be provided for each corresponding BAM: both for each normal sample and each tumor sample"))
    normal = ()
    samples = set()
    if names is None:
        normal = (normalbaf, os.path.splitext(os.path.basename(normalbaf))[0] )
        for tumor in tumors:
            samples.add((tumor, os.path.splitext(os.path.basename(tumor))[0]))
    else:
        normal = (normalbaf, names[0])
        for i in range(len(tumors)):
            samples.add((tumors[i], names[i+1]))

    # In default mode, check the existence and compatibility of samtools and bcftools
    samtools = os.path.join(args.samtools, "samtools")
    bcftools = os.path.join(args.bcftools, "bcftools")
    if sp.which(samtools) is None:
        raise ValueError(sp.error("{}samtools has not been found or is not executable!{}"))
    elif sp.which(bcftools) is None:
        raise ValueError(sp.error("{}bcftools has not been found or is not executable!{}"))
    elif not checkVersions(samtools, bcftools):
        raise ValueError(sp.error("The versions of samtools and bcftools are different! Please provide the tools with the same version to avoid inconsistent behaviors!{}"))

    # Check that SNP, reference, and region files exist when given in input
    if args.snps != None and not os.path.isfile(args.snps):
        raise ValueError(sp.error("The SNP file does not exist!"))
    if not os.path.isfile(args.reference):
        raise ValueError(sp.error("The provided file for human reference genome does not exist!"))
    if args.regions != None and not os.path.isfile(args.regions):
        raise ValueError(sp.error("The BED file of regions does not exist!"))
    elif args.regions is None:
        sp.log(msg="In case of WES data a BED file specified by --regions is REQUIRED, or the mincov parameter should be increased sufficiently to discard off-target regions\n", level="WARN")
    if args.snps != None and args.regions != None:
        raise ValueError(sp.error("Both SNP list and genomic regions have been provided, please provide only one of these!"))

    # Extract the names of the chromosomes and check their consistency across the given BAM files and the reference
    chromosomes = extractChromosomes(samtools, normal, samples, args.reference)

    if not args.processes > 0: raise ValueError(sp.error("The number of parallel processes must be greater than 0"))
    if not args.readquality >= 0: raise ValueError(sp.error("The read mapping quality must be positive"))
    if not args.basequality >= 0: raise ValueError(sp.error("The base quality quality must be positive"))
    if not 0 <= args.gamma and args.gamma <= 1: raise ValueError(sp.error("Gamma must be a floating value between 0 and 1"))
    if not 0 <= args.maxshift and args.maxshift <= 1: raise ValueError(sp.error("Max BAF shift must be a floating value between 0 and 0.5"))
    if not args.mincov >= 0: raise ValueError(sp.error("The minimum-coverage value must be positive"))
    if not args.maxcov >= 0: raise ValueError(sp.error("The maximum-coverage value must be positive"))
    if not args.snpquality >= 0: raise ValueError(sp.error("The QUAL value must be positive"))

    if args.verbose:
        sp.log(msg='stderr of samtools and bcftools will be collected in the following file "samtools.log"\n', level="WARN")
        with open("samtools.log", "w") as f: f.write("")

    return {"normal" : normal,
            "samples" : samples,
            "chromosomes" : chromosomes,
            "samtools" : samtools,
            "bcftools" : bcftools,
            "snps" : args.snps,
            "regions" : args.regions,
            "reference" : args.reference,
            "j" : args.processes,
            "q" : args.readquality,
            "Q" : args.basequality,
            "qual" : args.snpquality,
            "E" : args.newbaq,
            "gamma" : args.gamma,
            "maxshift" : args.maxshift,
            "mincov" : args.mincov,
            "maxcov" : args.maxcov,
            "outputNormal" : args.outputnormal,
            "outputTumors" : args.outputtumors,
            "outputSnps" : args.outputsnps,
            "verbose" : args.verbose}


def parse_bin_arguments(args=None):
    """
    Parse command line arguments
    Returns:
    """
    description = "Count the mapped sequencing reads in bins of fixed and given length, uniformly for a BAM file of a normal sample and one or more BAM files of tumor samples. This program supports both data from whole-genome sequencing (WGS) and whole-exome sequencing (WES), but the a BED file with targeted regions is required when considering WES."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-N","--normal", required=True, type=str, help="BAM file corresponding to matched normal sample")
    parser.add_argument("-T","--tumors", required=True, type=str, nargs='+', help="BAM files corresponding to samples from the same tumor")
    parser.add_argument("-b","--size", required=True, type=str, help="Size of the bins, specified as a full number or using the notations either \"kb\" or \"Mb\"")
    parser.add_argument("-S","--samples", required=False, default=None, type=str, nargs='+', help="Sample names for each BAM, given in the same order where the normal name is first (default: inferred from file names)")
    parser.add_argument("-st","--samtools", required=False, default="", type=str, help="Path to the directory to \"samtools\" executable, required in default mode (default: samtools is directly called as it is in user $PATH)")
    parser.add_argument("-r","--regions", required=False, default=None, type=str, help="BED file containing the a list of genomic regions to consider in the format \"CHR  START  END\", REQUIRED for WES data (default: none, consider entire genome)")
    parser.add_argument("-g","--reference", required=False, default=None, type=str, help="Reference genome, note that reference must be indexed and the dictionary must exist in the same directory with the same name and .dict extension")
    parser.add_argument("-j", "--processes", required=False, default=2, type=int, help="Number of available parallel processes (default: 2)")
    parser.add_argument("-q", "--readquality", required=False, default=0, type=int, help="Minimum mapping quality for an aligned read to be considered (default: 0)")
    parser.add_argument("-O", "--outputnormal", required=False, default=None, type=str, help="Filename of output for allele counts in the normal sample (default: standard output)")
    parser.add_argument("-o", "--outputtumors", required=False, default=None, type=str, help="Output filename for allele counts in tumor samples (default: standard output)")
    parser.add_argument("-t", "--outputtotal", required=False, default="total_read.counts", type=str, help="Output filename for total read counts in all tumor samples (default: \"total_read.counts\")")
    parser.add_argument("-v", "--verbose", action='store_true', default=False, required=False, help="Use verbose log messages")
    args = parser.parse_args(args)

    # Parse BAM files, check their existence, and infer or parse the corresponding sample names
    normalbaf = args.normal
    if not os.path.isfile(normalbaf): raise ValueError(sp.error("The specified normal BAM file does not exist"))
    tumors = args.tumors
    for tumor in tumors:
        if(not os.path.isfile(tumor)): raise ValueError(sp.error("The specified normal BAM file does not exist"))
    names = args.samples
    if names != None and (len(tumors)+1) != len(names):
        raise ValueError(sp.error("A sample name must be provided for each corresponding BAM: both for each normal sample and each tumor sample"))
    normal = ()
    samples = set()
    if names is None:
        normal = (normalbaf, os.path.splitext(os.path.basename(normalbaf))[0] )
        for tumor in tumors:
            samples.add((tumor, os.path.splitext(os.path.basename(tumor))[0]))
    else:
        normal = (normalbaf, names[0])
        for i in range(len(tumors)):
            samples.add((tumors[i], names[i+1]))

    # Check the region file
    if args.regions is not None and not os.path.isfile(args.regions):
        raise ValueError(sp.error("The specified region file does not exist"))

    # In default mode, check the existence and compatibility of samtools and bcftools
    samtools = os.path.join(args.samtools, "samtools")
    if sp.which(samtools) is None:
        raise ValueError(sp.error("samtools has not been found or is not executable!"))

    # Check and parse the given size
    size = 0
    try:
        if args.size[-2:] == "kb":
            size = int(args.size[:-2]) * 1000
        elif args.size[-2:] == "Mb":
            size = int(args.size[:-2]) * 1000000
        else:
            size = int(args.size)
    except:
        raise ValueError(sp.error("Size must be a number, optionally ending with either \"kb\" or \"Mb\"!"))

    # Check that either region file or available reference name are available
    if args.reference is None and args.regions is None:
        raise ValueError(sp.error("Please either provide a BED file of regions or specify a name of an available references for inferring maximum-chromosome lengths"))
    if args.reference is not None and not os.path.isfile(args.reference):
        raise ValueError(sp.error("The specified reference genome does not exist!"))
    refdict = os.path.splitext(args.reference)[0] + '.dict'
    if args.reference is not None and not os.path.isfile(refdict):
        raise ValueError(sp.error("The dictionary of the refence genome has not been found! Reference genome must be indeced and its dictionary must exist in the same directory with same name but extension .dict"))
    
    # Extract the names of the chromosomes and check their consistency across the given BAM files and the reference
    chromosomes = extractChromosomes(samtools, normal, samples)

    if not args.processes > 0: raise ValueError(sp.error("The number of parallel processes must be greater than 0"))
    if not args.readquality >= 0: raise ValueError(sp.error("The read mapping quality must be positive"))

    return {"normal" : normal,
            "samples" : samples,
            "chromosomes" : chromosomes,
            "samtools" : samtools,
            "regions" : args.regions,
            "size" : size,
            "reference" : args.reference,
            "refdict" : refdict,
            "j" : args.processes,
            "q" : args.readquality,
            "outputNormal" : args.outputnormal,
            "outputTumors" : args.outputtumors,
            "outputTotal" : args.outputtotal,
            "verbose" : args.verbose}


def parse_combbo_args(args=None):
    """
    Parse command line arguments
    Returns:
    """
    description = "Combine tumor bin counts, normal bin counts, and tumor allele counts to obtain the read-depth ratio and the mean B-allel frequency (BAF) of each bin. Optionally, the normal allele counts can be provided to add the BAF of each bin scaled by the normal BAF. The output is written on stdout."
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-c","--normalbins", required=True, type=str, help='Normal bin counts in the format "SAMPLE\tCHR\tSTART\tEND\tCOUNT"')
    parser.add_argument("-C","--tumorbins", required=True, type=str, help='Tumor bin counts in the format "SAMPLE\tCHR\tSTART\tEND\tCOUNT"')
    parser.add_argument("-B","--tumorbafs", required=True, type=str, help='Tumor allele counts in the format "SAMPLE\tCHR\tPOS\tREF-COUNT\tALT-COUNT"')
    parser.add_argument("-b","--normalbafs", required=False, default=None, type=str, help='Normal allele counts in the format "SAMPLE\tCHR\tPOS\tREF-COUNT\tALT-COUNT"')
    parser.add_argument("-d","--diploidbaf", type=float, required=False, default=0.1, help="Maximum diploid-BAF shift used to select the bins whose BAF should be normalized by the normal when normalbafs is given (default: 0.1)")
    parser.add_argument("-t","--totalcounts", required=False, default=None, type=str, help='Total read counts in the format "SAMPLE\tCOUNT" used to normalize by the different number of reads extracted from each sample (default: none)')
    parser.add_argument("-m","--mode", required=False, type=str, default="MIRROR", help='Mode name of the method to use for combining different SNPs covering the same bin (default: MIRROR):\n\n{}MIRROR{}: for each SNP the B allele corresponds to the allele in lower proportion.\n\n{}BINOMIAL_TEST{}: In addition to the MIRROR method each bin is tested to be copy neutral by asking if 0.5 is in the confidence interval of the corresponding BETA distribution, in that case bin has BAF equal to 0.\n\n{}MOMENTS{}: The method of moments is applied where the BAF of each SNP is computed by considering the allele in lower proportion as B and the mean across the SNPs is computed. \n\n'.format(sp.bcolors.BOLD, sp.bcolors.ENDC, sp.bcolors.BOLD, sp.bcolors.ENDC, sp.bcolors.BOLD, sp.bcolors.ENDC))
    parser.add_argument("-g","--gamma",type=float,required=False, default=0.05, help='Confidence level used to determine if a bin is copy neutral with BAF of 0.5 in the BINOMIAL_TEST mode (default: 0.05)')
    parser.add_argument("-e","--seed", type=int, required=False, default=0, help='Random seed used for the normal distributions used in the clouds (default: 0)')
    parser.add_argument("-s","--bootstrap", type=int, required=False, default=100, help='Number of draws for bootstrapping each SNP. Bootstrap significantly helps to estimate the BAF of each bin by combining the corresponding SNPs when SAMPLE or AVERAGE modes are used. (default: 100)')
    parser.add_argument("-dB","--bafdeviation", type=float, required=False, default=0.02, help='Standard deviation of the BAFs used to generate the points in the clouds (default: 0.002)')
    parser.add_argument("-v", "--verbose", action='store_true', default=False, required=False, help="Use verbose log messages")
    parser.add_argument("-r", "--disablebar", action='store_true', default=False, required=False, help="Disable progress bar")
    args = parser.parse_args(args)

    if not os.path.isfile(args.normalbins):
        raise ValueError(sp.error("The specified file for normal bin counts does not exist!"))
    if not os.path.isfile(args.tumorbins):
        raise ValueError(sp.error("The specified file for tumor bin counts does not exist!"))
    if not os.path.isfile(args.normalbins):
        raise ValueError(sp.error("The specified file for normal bin counts does not exist!"))
    if args.normalbafs is not None and not os.path.isfile(args.normalbafs):
        raise ValueError(sp.error("The specified file for normal baf does not exist!"))
    if not 0.0 <= args.diploidbaf <= 0.5:
        raise ValueError(sp.error("The specified maximum for diploid-BAF shift must be a value in [0.0, 0.5]"))
    if args.totalcounts is not None and not os.path.isfile(args.totalcounts):
        raise ValueError(sp.error("The specified file for total read counts does not exist!"))
    if args.mode not in ["BINOMIAL_TEST", "MOMENTS", "MIRROR"]:
        raise ValueError(sp.error("The specified mode does not exist!"))
    if not 0.0 <= args.gamma <= 0.1:
        raise ValueError(sp.error("The specified gamma must be a value in [0.0, 0.1]"))
    if args.seed < 0:
        raise ValueError(sp.error("Seed parameter must be positive!"))
    if args.bootstrap < 0:
        raise ValueError(sp.error("Bootstrap parameter must be positive!"))

    return {"normalbins" : args.normalbins,
            "tumorbins" : args.tumorbins,
            "tumorbafs" : args.tumorbafs,
            "normalbafs" : args.normalbafs,
            "diploidbaf" : args.diploidbaf,
            "totalcounts" : args.totalcounts,
            "mode" : args.mode,
            "gamma" : args.gamma,
            "seed" : args.seed,
            "bootstrap" : args.bootstrap,
            "bafsd" : args.bafdeviation,
            "verbose" : args.verbose,
            "disable" : args.disablebar}


def parse_clubb_args(args=None):
    """
    Parse command line arguments
    Returns:
    """
    description = "Combine tumor bin counts, normal bin counts, and tumor allele counts to obtain the read-depth ratio and the mean B-allel frequency (BAF) of each bin. Optionally, the normal allele counts can be provided to add the BAF of each bin scaled by the normal BAF. The output is written on stdout."
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("BBFILE", help="A BB file containing a line for each bin in each sample and the corresponding values of read-depth ratio and B-allele frequency (BAF)")
    parser.add_argument("-o", "--outsegments", required=False, default=None, type=str, help="Output filename for the segments computed by clustering bins (default: stdout)")
    parser.add_argument("-O", "--outbins", required=False, default=None, type=str, help="Output filename for a BB file adding the clusters (default: stdout)")
    parser.add_argument("-d","--diploidbaf", type=float, required=False, default=None, help="Maximum diploid-BAF shift used to determine the largest copy-neutral cluster and to rescale all the cluster inside this threshold accordingly (default: None, scaling is not performed)")
    parser.add_argument("-tR","--tolerancerdr", type=float, required=False, default=0.0, help='Refine the clustering merging the clusters with this maximum difference in RDR values (default: None, gurobipy required)')
    parser.add_argument("-tB","--tolerancebaf", type=float, required=False, default=0.0, help='Refine the clustering merging the clusters with this maximum difference in BAF values (default: None, gurobipy required)')
    parser.add_argument("-u","--bootclustering", type=int, required=False, default=0, help='Number of points to add for bootstraping each bin to improve the clustering. Each point is generated by drawing its values from a normal distribution centered on the values of the bin. This can help the clustering when the input number of bins is low (default: 0)')
    parser.add_argument("-dR","--ratiodeviation", type=float, required=False, default=0.02, help='Standard deviation of the read ratios used to generate the points in the clouds (default: 0.02)')
    parser.add_argument("-dB","--bafdeviation", type=float, required=False, default=0.02, help='Standard deviation of the BAFs used to generate the points in the clouds (default: 0.02)')
    parser.add_argument("-e","--seed", type=int, required=False, default=0, help='Random seed used for the normal distributions used in the clouds (default: 0)')
    parser.add_argument("-K","--initclusters", type=int, required=False, default=15, help="The initial number of clusters (default: 15)")
    parser.add_argument("-sf","--tuning", type=float, required=False, default=0.01, help="Tuning parameter for clustering; used to determine initial size of distribution covariances. Small sf indicates a belief that clusters are of small size (default: 0.01)")
    parser.add_argument("-R","--restarts", type=int, required=False, default=10, help="Number of restarts performed by the clustering to choose the best (default: 10)")
    parser.add_argument("-v","--verbose", action='store_true', default=False, required=False, help="Use verbose log messages")
    parser.add_argument("--disablebar", action='store_true', default=False, required=False, help="Disable progress bar")
    args = parser.parse_args(args)

    if not os.path.isfile(args.BBFILE):
        raise ValueError(sp.error("The specified BB file does not exist!"))
    if args.diploidbaf != None and not 0.0 <= args.diploidbaf <= 0.5:
        raise ValueError(sp.error("The specified maximum for diploid-BAF shift must be a value in [0.0, 0.5]"))
    if args.tolerancerdr < 0:
        raise ValueError(sp.error("Tolerance-RDR parameter must be positive!"))
    if args.tolerancebaf < 0:
        raise ValueError(sp.error("Tolerance-BAF parameter must be positive!"))
    if args.bootclustering < 0:
        raise ValueError(sp.error("Bootclustering parameter must be positive!"))
    if args.ratiodeviation < 0:
        raise ValueError(sp.error("Ratio-deviation parameter must be positive!"))
    if args.bafdeviation < 0:
        raise ValueError(sp.error("BAF-deviation parameter must be positive!"))
    if args.seed < 0:
        raise ValueError(sp.error("Seed parameter must be positive!"))
    if args.initclusters < 0:
        raise ValueError(sp.error("Init-cluster parameter must be positive!"))
    if args.tuning < 0:
        raise ValueError(sp.error("Tuning parameter must be positive!"))
    if args.restarts < 0:
        raise ValueError(sp.error("Number of restarts must be positive!"))

    return {"bbfile" : args.BBFILE,
            "cloud" : args.bootclustering,
            "diploidbaf" : args.diploidbaf,
            "rdtol" : args.tolerancerdr,
            "baftol" : args.tolerancebaf,
            "ratiodeviation" : args.ratiodeviation,
            "bafdeviation" : args.bafdeviation,
            "seed" : args.seed,
            "initclusters" : args.initclusters,
            "tuning" : args.tuning,
            "restarts" : args.restarts,
            "verbose" : args.verbose,
            "disable" : args.disablebar,
            "outsegments" : args.outsegments,
            "outbins" : args.outbins}


def parse_bbot_args(args=None):
    """
    Parse command line arguments
    Returns:
    """
    description = "Generate plots for read-depth ratio (RD), B-allele frequency (BAF), and clusters for genomic bins in multiple samples using .bb, .cbb, .seg files."
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("INPUT", help='Input BBC file with RDR and BAF')
    parser.add_argument("-c", "--command", required=False, type=str, default=None, help="The command determining the plots to generate (default: all)\n\n\t{}RD{}: Plot the read-depth ratio (RD) values of the genomes for each sample.\n\n\t{}CRD{}: Plot the read-depth ratio (CRD) values of the genomes for each sample colored by corresponding cluster.\n\n\t{}BAF{}: Plot the B-allele frequency (BAF) values of the genomes for each sample.\n\n\t{}CBAF{}: Plot BAF values for each sample colored by corresponding cluster.\n\n\t{}BB{}: Plot jointly the values of read-depth ratio (RD) and B-allele frequency (BAF) for each bin in all samples and their density.\n\n\t{}CBB{}: Plot jointly the values of read-depth ratio (RD) and B-allele frequency (BAF) for each bin in all samples by coloring the bins depending on their cluster.\n\n\t{}CLUSTER{}: Plot jointly the values of read-depth ratio (RD) and B-allele frequency (BAF) for each cluster in all samples where the size of the markers is proportional to the number of bins.".format(sp.bcolors.BOLD, sp.bcolors.ENDC, sp.bcolors.BOLD, sp.bcolors.ENDC, sp.bcolors.BOLD, sp.bcolors.ENDC, sp.bcolors.BOLD, sp.bcolors.ENDC, sp.bcolors.BOLD, sp.bcolors.ENDC, sp.bcolors.BOLD, sp.bcolors.ENDC, sp.bcolors.BOLD, sp.bcolors.ENDC, sp.bcolors.BOLD, sp.bcolors.ENDC))
    parser.add_argument("-s", "--segfile", required=False, type=str, default=None, help="When the corresponding seg file is provided the clusters are also plotted (default: none)")
    parser.add_argument("-m","--colormap", required=False, type=str, default="tab20",help='Colormap to use for the colors in the plots, the available colormaps are the following {Set1, Set2, Paired, Dark2, tab10, tab20}')
    parser.add_argument("-tC","--chrthreshold", required=False, type=int, default=None,help='Only covering at least this number of chromosomes are considered (default: None)')
    parser.add_argument("-tS","--sizethreshold", required=False, type=float, default=None,help='Only covering at least this genome proportion (default: None)')
    parser.add_argument("--resolution", required=False, default=None, type=int, help='Resolution of bins (default: bins are not merged)')
    parser.add_argument("--xmin", required=False, default=None, type=float, help='Minimum value on x-axis for supported plots (default: inferred from data)')
    parser.add_argument("--xmax", required=False, default=None, type=float, help='Maximum value on x-axis for supported plots (default: inferred from data)')
    parser.add_argument("--ymin", required=False, default=None, type=float, help='Minimum value on y-axis for supported plots (default: inferred from data)')
    parser.add_argument("--ymax", required=False, default=None, type=float, help='Maximum value on y-axis for supported plots (default: inferred from data)')
    parser.add_argument("--figsize", required=False, default=None, type=str, help='Size of the plotted figures in the form "(X-SIZE, Y-SIZE)"')
    parser.add_argument("--markersize", required=False, default=0, type=int, help='Size of the markers (default: values inferred for each plot)')
    parser.add_argument("--colwrap", required=False, default=2, type=int, help='Wrapping the plots in this number of columnes (default: 2)')
    parser.add_argument("--fontscale", required=False, default=1, type=float, help='Font scale (default: 1)')
    parser.add_argument("-x","--rundir", required=False, default='./', type=str, help='Running dirrectory where output the results (default: current directory)')
    parser.add_argument("--pdf", action='store_true', default=False, required=False, help="Output the bb_clustered figure in PDF format (default: PNG)")
    parser.add_argument("--dpi", required=False, default=900, type=int, help='DPI of PNG images (default: 900)')
    args = parser.parse_args(args)

    if not os.path.isfile(args.INPUT):
        raise ValueError(sp.error("The specified BB file does not exist!"))
    if args.command is not None and args.command not in {"RD", "CRD", "BAF", "CBAF", "BB", "CBB", "CLUSTER"}:
        raise ValueError(sp.error("Unrecognized COMMAND!"))
    if args.segfile is not None and not os.path.isfile(args.segfile):
        raise ValueError(sp.error("Specified seg-file does not exist!"))
    if args.colormap not in {"Set1", "Set2", "Set3", "Paired", "Accent", "Dark2", "tab10", "tab20", "husl", "hls", "muted", "colorblind", "Pastel1", "Pastel2"}:
        raise ValueError(sp.error("Unrecognized colormap!"))
    if args.resolution is not None and args.resolution < 1:
        raise ValueError(sp.error("Resolution must be greater than 1!"))
    if args.chrthreshold is not None and not( 0 <= args.chrthreshold <= 22):
        raise ValueError(sp.error("The chromosome threshold must be a integer in \{0, ..., 22\}!"))
    if args.colwrap < 1:
        raise ValueError(sp.error("Colwrap must be grater than 1!"))
    if args.fontscale < 0:
        raise ValueError(sp.error("Font scale must be positive!"))
    if args.sizethreshold is not None and not( 0.0 <= args.sizethreshold <= 1.0):
        raise ValueError(sp.error("The size threshold must be a real value in [0, 1]!"))
    if args.figsize is not None:
        try:
            parsed = args.figsize.strip().split(',')
            figsize = (float(parsed[0]), float(parsed[1]))
        except:
            raise ValueError(sp.error("Wrong format of figsize!"))
    else:
        figsize = None
    if not os.path.isdir(args.rundir):
        raise ValueError(sp.error("Running directory either does not exists or is not a directory!"))

    return {"input" : args.INPUT,
            "command" : args.command,
            "segfile" : args.segfile,
            "ct" : args.chrthreshold,
            "st" : args.sizethreshold,
            "cmap" : args.colormap,
            'resolution' : args.resolution,
            "xmin" : args.xmin,
            "xmax" : args.xmax,
            "ymin" : args.ymin,
            "ymax" : args.ymax,
            "x" : args.rundir,
            "figsize" : figsize,
            "markersize" : args.markersize,
            "colwrap" : args.colwrap,
            "fontscale" : args.fontscale,
            "pdf" : args.pdf,
            "dpi" : args.dpi
    }


def extractChromosomes(samtools, normal, tumors, reference=None):
    # Read the names of sequences in normal BAM file
    normal_sq = getSQNames(samtools, normal[0])

    # Extract only the names of chromosomes in standard formats
    chrm = set()
    no_chrm = set()
    for i in range(1, 23):
        if str(i) in normal_sq:
            no_chrm.add(str(i))
        elif "chr" + str(i) in normal_sq:
            chrm.add("chr" + str(i))
        else:
            sys.stderr.write("WARNING: a chromosome named either {} or a variant of CHR{} cannot be found in the normal BAM file\n".format(i, i))

    if len(chrm) == 0 and len(no_chrm) == 0: raise ValueError("No chromosomes found in the normal BAM")
    chromosomes = set()
    if len(chrm) > len(no_chrm):
        chromosomes = chrm
    else:
        chromosomes = no_chrm

    # Check that chromosomes with the same names are present in each tumor BAM contain
    for tumor in tumors:
        tumor_sq = getSQNames(samtools, tumor[0])
        if not chromosomes <= tumor_sq:
            sys.stderr.write("WARNING: chromosomes {} are not present in the tumor sample {}\n".format(chromosomes - tumor_sq, tumor))

    # Check consistency of chromosome names with the reference
    if reference is not None:
        stdout, stderr = subprocess.Popen("grep -e \"^>\" {}".format(reference), stdout=subprocess.PIPE, shell=True).communicate()
        if stderr is not None:
            raise ValueError("Error in reading the reference: {}".format(reference))
        else:
            ref = set(c[1:].strip().split()[0] for c in stdout.strip().split('\n'))
        if not(chromosomes <= ref):
            raise ValueError("The given reference cannot be used because the chromosome names are inconsistent!\nChromosomes found in BAF files: {}\nChromosomes with the same name found in reference genome: {}".format(chromosomes, ref))

    return sorted(list(chromosomes), key=sp.numericOrder)


def getSQNames(samtools, bamfile):
    header, stderr = subprocess.Popen([samtools, "view", "-H", bamfile], stdout=subprocess.PIPE, shell=False).communicate()
    if stderr is not None:
        raise ValueError("The header of the normal-sample BAM cannot be read with samtools!")
    names = set()
    for line in header.strip().split('\n'):
        line = line.split()
        if line[0] == '@SQ':
            names.add(line[1].split(':')[1])
    return names


def checkVersions(samtools, bcftools):
    samtools_version, samtools_stderr = subprocess.Popen([samtools, "--version-only"], stdout=subprocess.PIPE, shell=False).communicate()
    bcftools_version, bcftools_stderr = subprocess.Popen([bcftools, "--version-only"], stdout=subprocess.PIPE, shell=False).communicate()
    return samtools_version == bcftools_version and samtools_stderr is None and bcftools_stderr is None


def parseRegions(region_file, chromosomes):
    res = {}
    for chro in chromosomes:
        res[chro] = []
    nofound = set()
    with open(region_file, 'r') as f:
        for line in f:
            split = line.strip().split()
            chro = split[0]
            if chro in chromosomes:
                res[chro].append((int(split[1]), int(split[2])))
            else:
                nofound.add(chro)

    for c in nofound:
        sp.log(msg="The chromosome {} present in the provided regions is non-autosome or is not present in the given BAM files\n".format(c), level="WARN")

    for key in res:
        res[key].sort(key=lambda x : x[0])
        if not all(a[0]<=a[1] and a[1]<=b[0] and b[0]<=b[1] for a, b in zip(res[key], res[key][1:])):
            raise ValueError(sp.error("The regions provided for chromosome {} are non-disjoint or a region start is greater than corresponding region end".format(key)))

    return res
