import sys
import os.path
import argparse
import subprocess

from . import Supporting as sp
from hatchet import config, __version__


def parse_genotype_snps_arguments(args=None):
    """
    Parse command line arguments
    Returns:
    """
    description = "Genotype and call SNPs in a matched-normal sample."
    parser = argparse.ArgumentParser(prog='hatchet genotype-snps', description=description)
    parser.add_argument("-N","--normal", required=True, type=str, help="BAM file corresponding to matched normal sample")
    parser.add_argument("-r","--reference", required=True, type=str, help="Human reference genome of BAMs")
    parser.add_argument("-st","--samtools", required=False, default=config.paths.samtools, type=str, help="Path to the directory to \"samtools\" executable, required in default mode (default: samtools is directly called as it is in user $PATH)")
    parser.add_argument("-bt","--bcftools", required=False, default=config.paths.bcftools, type=str, help="Path to the directory of \"bcftools\" executable, required in default mode (default: bcftools is directly called as it is in user $PATH)")
    parser.add_argument("-R","--snps", required=False, default=config.genotype_snps.snps, type=str, help="List of SNPs to consider in the normal sample (default: heterozygous SNPs are inferred from the normal sample)")
    parser.add_argument("-j", "--processes", required=False, default=config.genotype_snps.processes, type=int, help="Number of available parallel processes (default: 2)")
    parser.add_argument("-q", "--readquality", required=False, default=config.genotype_snps.readquality, type=int, help="Minimum mapping quality for an aligned read to be considered (default: 0)")
    parser.add_argument("-Q", "--basequality", required=False, default=config.genotype_snps.basequality, type=int, help="Minimum base quality for a base to be considered (default: 11)")
    parser.add_argument("-c", "--mincov", required=False, default=config.genotype_snps.mincov, type=int, help="Minimum coverage for SNPs to be considered (default: 0)")
    parser.add_argument("-C", "--maxcov", required=False, default=config.genotype_snps.maxcov, type=int, help="Maximum coverage for SNPs to be considered (default: 1000, suggested: twice the values of expected average coverage to avoid aligning artefacts)")
    parser.add_argument("-E","--newbaq", required=False, action='store_true', default=config.genotype_snps.newbaq, help="Recompute alignment of reads on the fly during SNP calling (default: false)")
    parser.add_argument("-o", "--outputsnps", required=False, default=config.genotype_snps.outputsnps, type=str, help="Output folder for SNPs separated by chromosome (default: ./)")
    parser.add_argument("-v", "--verbose", action='store_true', default=config.genotype_snps.verbose, required=False, help="Use verbose log messages")
    parser.add_argument("-V", "--version", action='version', version=f'%(prog)s {__version__}')
    args = parser.parse_args(args)

    # Parse BAM files, check their existence, and infer or parse the corresponding sample names
    normalbaf = args.normal
    if not os.path.isfile(os.path.abspath(normalbaf)):
        raise ValueError(sp.error("The specified normal BAM file does not exist"))
    normal = (os.path.abspath(normalbaf), 'Normal')

    # In default mode, check the existence and compatibility of samtools and bcftools
    samtools = os.path.join(args.samtools, "samtools")
    bcftools = os.path.join(args.bcftools, "bcftools")
    if sp.which(samtools) is None:
        raise ValueError(sp.error("{}samtools has not been found or is not executable!{}"))
    elif sp.which(bcftools) is None:
        raise ValueError(sp.error("{}bcftools has not been found or is not executable!{}"))

    # Check that SNP, reference, and region files exist when given in input
    if not os.path.isfile(args.reference):
        raise ValueError(sp.error("The provided file for human reference genome does not exist!"))
    if args.outputsnps != None and not os.path.isdir(args.outputsnps):
        raise ValueError(sp.error("The folder for output SNPs does not exist!"))
    if args.snps != None and len(args.snps) < 2:
        args.snps = None
    if args.snps != None and not (os.path.isfile(args.snps) or sp.urlexists(args.snps)):
        raise ValueError(sp.error("The provided list of SNPs does not exist!"))

    # Extract the names of the chromosomes and check their consistency across the given BAM files and the reference
    chromosomes = extractChromosomes(samtools, normal, [], args.reference)

    if not args.processes > 0: raise ValueError(sp.error("The number of parallel processes must be greater than 0"))
    if not args.readquality >= 0: raise ValueError(sp.error("The read mapping quality must be positive"))
    if not args.basequality >= 0: raise ValueError(sp.error("The base quality quality must be positive"))
    if not args.mincov >= 0: raise ValueError(sp.error("The minimum-coverage value must be positive"))
    if not args.maxcov >= 0: raise ValueError(sp.error("The maximum-coverage value must be positive"))

    if args.verbose:
        sp.log(msg='stderr of samtools and bcftools will be collected in the following file "samtools.log"\n', level="WARN")
        with open("samtools.log", "w") as f: f.write("")

    return {"normal" : normal,
            "chromosomes" : chromosomes,
            "samtools" : samtools,
            "bcftools" : bcftools,
            "snps" : args.snps,
            "reference" : args.reference,
            "j" : args.processes,
            "q" : args.readquality,
            "Q" : args.basequality,
            "E" : args.newbaq,
            "mincov" : args.mincov,
            "maxcov" : args.maxcov,
            "outputsnps" : os.path.abspath(args.outputsnps),
            "verbose" : args.verbose}


def parse_count_alleles_arguments(args=None):
    """
    Parse command line arguments
    Returns:
    """
    description = "Count the A/B alleles from a matched-normal BAM file and multiple tumor BAM files in specified SNP positions or estimated heterozygous SNPs in the normal genome. This tool can be applied both to whole-genome sequencing (WGS) data or whole-exome sequencing (WES) data, but coding regions must be specified as a BED file in the case of WES."
    parser = argparse.ArgumentParser(prog='hatchet count-alleles', description=description)
    parser.add_argument("-N","--normal", required=True, type=str, help="BAM file corresponding to matched normal sample")
    parser.add_argument("-T","--tumors", required=True, type=str, nargs='+', help="BAM files corresponding to samples from the same tumor")
    parser.add_argument("-r","--reference", required=True, type=str, help="Human reference genome of BAMs")
    parser.add_argument("-L","--snps", required=True, type=str, nargs='+', help="List of SNPs to consider in the normal sample")
    parser.add_argument("-S","--samples", required=False, default=config.count_alleles.samples, type=str, nargs='+', help="Sample names for each BAM (given in the same order where the normal name is first)")
    parser.add_argument("-st","--samtools", required=False, default=config.paths.samtools, type=str, help="Path to the directory to \"samtools\" executable, required in default mode (default: samtools is directly called as it is in user $PATH)")
    parser.add_argument("-bt","--bcftools", required=False, default=config.paths.bcftools, type=str, help="Path to the directory of \"bcftools\" executable, required in default mode (default: bcftools is directly called as it is in user $PATH)")
    parser.add_argument("-e","--regions", required=False, default=config.count_alleles.regions, type=str, help="BED file containing the a list of genomic regions to consider in the format \"CHR  START  END\", REQUIRED for WES data with coding regions (default: none, consider entire genome)")
    parser.add_argument("-j", "--processes", required=False, default=config.count_alleles.processes, type=int, help="Number of available parallel processes (default: 2)")
    parser.add_argument("-q", "--readquality", required=False, default=config.count_alleles.readquality, type=int, help="Minimum mapping quality for an aligned read to be considered (default: 0)")
    parser.add_argument("-Q", "--basequality", required=False, default=config.count_alleles.basequality, type=int, help="Minimum base quality for a base to be considered (default: 11)")
    parser.add_argument("-U", "--snpquality", required=False, default=config.count_alleles.snpquality, type=int, help="Minimum SNP-variant quality, QUAL, for a variant to be considered (default: 11)")
    parser.add_argument("-g", "--gamma", required=False, default=config.count_alleles.gamma, type=float, help="Level of confidence to determine heterozigosity of SNPs (default: 0.05)")
    parser.add_argument("-b", "--maxshift", required=False, default=config.count_alleles.maxshift, type=float, help="Maximum allowed absolute difference of BAF from 0.5 for selected heterozygous SNPs in the normal sample (default: 0.5)")
    parser.add_argument("-c", "--mincov", required=False, default=config.count_alleles.mincov, type=int, help="Minimum coverage for SNPs to be considered (default: 0)")
    parser.add_argument("-C", "--maxcov", required=False, default=config.count_alleles.maxcov, type=int, help="Maximum coverage for SNPs to be considered (default: 1000, suggested: twice the values of expected average coverage to avoid aligning artefacts)")
    parser.add_argument("-E","--newbaq", required=False, action='store_true', default=config.count_alleles.newbaq, help="Recompute alignment of reads on the fly during SNP calling (default: false)")
    parser.add_argument("-O", "--outputnormal", required=False, default=config.count_alleles.outputnormal, type=str, help="Filename of output for allele counts in the normal sample (default: standard output)")
    parser.add_argument("-o", "--outputtumors", required=False, default=config.count_alleles.outputtumors, type=str, help="Output filename for allele counts in tumor samples (default: standard output)")
    parser.add_argument("-l", "--outputsnps", required=False, default=config.count_alleles.outputsnps, type=str, help="Output directory for lists of selected SNPs (default: output dir of patient)")
    parser.add_argument("-v", "--verbose", action='store_true', default=config.count_alleles.verbose, required=False, help="Use verbose log messages")
    parser.add_argument("-V", "--version", action='version', version=f'%(prog)s {__version__}')
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

    # Check that SNP, reference, and region files exist when given in input
    snplists = {}
    for f in args.snps:
        if not os.path.isfile(f):
            raise ValueError(sp.error("The specified SNP file {} does not exist!".format(f)))
        snplists = {os.path.basename(f).split('.')[0] : f for f in args.snps}
    if not os.path.isfile(args.reference):
        raise ValueError(sp.error("The provided file for human reference genome does not exist!"))
    if args.regions != None and not os.path.isfile(args.regions):
        raise ValueError(sp.error("The BED file of regions does not exist!"))
    #elif args.regions is None:
    #    sp.log(msg="In case of WES data a BED file specified by --regions is REQUIRED, or the mincov parameter should be increased sufficiently to discard off-target regions\n", level="WARN")
    if args.snps != None and args.regions != None:
        raise ValueError(sp.error("Both SNP list and genomic regions have been provided, please provide only one of these!"))

    # Extract the names of the chromosomes and check their consistency across the given BAM files and the reference
    chromosomes = extractChromosomes(samtools, normal, samples, args.reference)
    #for c in chromosomes:
    #    if c not in snplists:
    #        raise ValueError(sp.error('The SNP file for analyzed chromosome {} was expected with name {}.* but not found in the provided list!'.format(c, c)))
    snplists = {c : snplists.get(c, []) for c in chromosomes}

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
            "snps" : snplists,
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


def parse_count_reads_arguments(args=None):
    """
    Parse command line arguments
    Returns:
    """
    description = "Count the mapped sequencing reads in bins of fixed and given length, uniformly for a BAM file of a normal sample and one or more BAM files of tumor samples. This program supports both data from whole-genome sequencing (WGS) and whole-exome sequencing (WES), but the a BED file with targeted regions is required when considering WES."
    parser = argparse.ArgumentParser(prog='hatchet count-reads', description=description)
    parser.add_argument("-N","--normal", required=True, type=str, help="BAM file corresponding to matched normal sample")
    parser.add_argument("-T","--tumors", required=True, type=str, nargs='+', help="BAM files corresponding to samples from the same tumor")
    parser.add_argument("-b","--size", required=True, type=str, help="Size of the bins, specified as a full number or using the notations either \"kb\" or \"Mb\"")
    parser.add_argument("-S","--samples", required=False, default=config.count_reads.samples, type=str, nargs='+', help="Sample names for each BAM, given in the same order where the normal name is first (default: inferred from file names)")
    parser.add_argument("-st","--samtools", required=False, default=config.paths.samtools, type=str, help="Path to the directory to \"samtools\" executable, required in default mode (default: samtools is directly called as it is in user $PATH)")
    parser.add_argument("-r","--regions", required=False, default=config.count_reads.regions, type=str, help="BED file containing the a list of genomic regions to consider in the format \"CHR  START  END\", REQUIRED for WES data (default: none, consider entire genome)")
    parser.add_argument("-g","--reference", required=False, default=config.paths.reference, type=str, help="Reference genome, note that reference must be indexed and the dictionary must exist in the same directory with the same name and .dict extension")
    parser.add_argument("-j", "--processes", required=False, default=config.count_reads.processes, type=int, help="Number of available parallel processes (default: 2)")
    parser.add_argument("-q", "--readquality", required=False, default=config.count_reads.readquality, type=int, help="Minimum mapping quality for an aligned read to be considered (default: 0)")
    parser.add_argument("-O", "--outputnormal", required=False, default=config.count_reads.outputnormal, type=str, help="Filename of output for allele counts in the normal sample (default: standard output)")
    parser.add_argument("-o", "--outputtumors", required=False, default=config.count_reads.outputtumors, type=str, help="Output filename for allele counts in tumor samples (default: standard output)")
    parser.add_argument("-t", "--outputtotal", required=False, default=config.count_reads.outputtotal, type=str, help="Output filename for total read counts in all tumor samples (default: \"total_read.counts\")")
    parser.add_argument("-v", "--verbose", action='store_true', default=config.count_reads.verbose, required=False, help="Use verbose log messages")
    parser.add_argument("-V", "--version", action='version', version=f'%(prog)s {__version__}')
    args = parser.parse_args(args)

    # Parse BAM files, check their existence, and infer or parse the corresponding sample names
    normalbaf = args.normal
    if not os.path.isfile(normalbaf): raise ValueError(sp.error("The specified normal BAM file does not exist"))
    tumors = args.tumors
    for tumor in tumors:
        if(not os.path.isfile(tumor)): raise ValueError(sp.error(f"The specified tumor BAM file {tumor} does not exist"))
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


def parse_combine_counts_args(args=None):
    """
    Parse command line arguments
    Returns:
    """
    description = "Combine tumor bin counts, normal bin counts, and tumor allele counts to obtain the read-depth ratio and the mean B-allel frequency (BAF) of each bin. Optionally, the normal allele counts can be provided to add the BAF of each bin scaled by the normal BAF. The output is written on stdout."
    parser = argparse.ArgumentParser(prog='hatchet combine-counts', description=description, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-c","--normalbins", required=True, type=str, help='Normal bin counts in the format "SAMPLE\tCHR\tSTART\tEND\tCOUNT"')
    parser.add_argument("-C","--tumorbins", required=True, type=str, help='Tumor bin counts in the format "SAMPLE\tCHR\tSTART\tEND\tCOUNT"')
    parser.add_argument("-B","--tumorbafs", required=True, type=str, help='Tumor allele counts in the format "SAMPLE\tCHR\tPOS\tREF-COUNT\tALT-COUNT"')
    parser.add_argument("-p","--phase", required=False, default=config.combine_counts.phase, type=str, help='Phasing of heterozygous germline SNPs in the format "CHR\tPOS\t<string containing 0|1 or 1|0>"')
    parser.add_argument("-l","--blocklength", required=False, default=config.combine_counts.blocklength, type=str, help="Size of the haplotype blocks, specified as a full number or using the notations either \"kb\" or \"Mb\" (default: 50kb)")
    parser.add_argument("-d","--diploidbaf", type=float, required=False, default=config.combine_counts.diploidbaf, help="Maximum diploid-BAF shift used to select the bins whose BAF should be normalized by the normal when normalbafs is given (default: 0.1)")
    parser.add_argument("-t","--totalcounts", required=False, default=config.combine_counts.totalcounts, type=str, help='Total read counts in the format "SAMPLE\tCOUNT" used to normalize by the different number of reads extracted from each sample (default: none)')
    parser.add_argument("-g","--gamma",type=float,required=False, default=config.combine_counts.gamma, help='Confidence level used to determine if a bin is copy neutral with BAF of 0.5 in the BINOMIAL_TEST mode (default: 0.05)')
    parser.add_argument("-e","--seed", type=int, required=False, default=config.combine_counts.seed, help='Random seed used for the normal distributions used in the clouds (default: 0)')
    parser.add_argument("-v", "--verbose", action='store_true', default=config.combine_counts.verbose, required=False, help="Use verbose log messages")
    parser.add_argument("-r", "--disablebar", action='store_true', default=config.combine_counts.disablebar, required=False, help="Disable progress bar")
    parser.add_argument("-V", "--version", action='version', version=f'%(prog)s {__version__}')
    args = parser.parse_args(args)

    if not os.path.isfile(args.normalbins):
        raise ValueError(sp.error("The specified file for normal bin counts does not exist!"))
    if not os.path.isfile(args.tumorbins):
        raise ValueError(sp.error("The specified file for tumor bin counts does not exist!"))
    if not os.path.isfile(args.normalbins):
        raise ValueError(sp.error("The specified file for normal bin counts does not exist!"))
    if args.phase == 'None':
        args.phase = None
    if args.phase is not None and not os.path.isfile(args.phase):
        raise ValueError(sp.error("The specified file for phase does not exist!"))
    if not 0.0 <= args.diploidbaf <= 0.5:
        raise ValueError(sp.error("The specified maximum for diploid-BAF shift must be a value in [0.0, 0.5]"))
    if args.totalcounts is not None and not os.path.isfile(args.totalcounts):
        raise ValueError(sp.error("The specified file for total read counts does not exist!"))
    if not 0.0 <= args.gamma <= 0.1:
        raise ValueError(sp.error("The specified gamma must be a value in [0.0, 0.1]"))
    if args.seed < 0:
        raise ValueError(sp.error("Seed parameter must be positive!"))

    size = 0
    try:
        if args.blocklength[-2:] == "kb":
            size = int(args.blocklength[:-2]) * 1000
        elif args.blocklength[-2:] == "Mb":
            size = int(args.blocklength[:-2]) * 1000000
        else:
            size = int(args.blocklength)
    except:
        raise ValueError(sp.error("Size must be a number, optionally ending with either \"kb\" or \"Mb\"!"))

    return {"normalbins" : args.normalbins,
            "tumorbins" : args.tumorbins,
            "tumorbafs" : args.tumorbafs,
            "phase" : args.phase,
            "block" : size,
            "diploidbaf" : args.diploidbaf,
            "totalcounts" : args.totalcounts,
            "gamma" : args.gamma,
            "seed" : args.seed,
            "verbose" : args.verbose,
            "disable" : args.disablebar}


def parse_cluster_bins_args(args=None):
    """
    Parse command line arguments
    Returns:
    """
    description = "Combine tumor bin counts, normal bin counts, and tumor allele counts to obtain the read-depth ratio and the mean B-allel frequency (BAF) of each bin. Optionally, the normal allele counts can be provided to add the BAF of each bin scaled by the normal BAF. The output is written on stdout."
    parser = argparse.ArgumentParser(prog='hatchet cluster-bins', description=description, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("BBFILE", help="A BB file containing a line for each bin in each sample and the corresponding values of read-depth ratio and B-allele frequency (BAF)")
    parser.add_argument("-o", "--outsegments", required=False, default=config.cluster_bins.outsegments, type=str, help="Output filename for the segments computed by clustering bins (default: stdout)")
    parser.add_argument("-O", "--outbins", required=False, default=config.cluster_bins.outbins, type=str, help="Output filename for a BB file adding the clusters (default: stdout)")
    parser.add_argument("-d","--diploidbaf", type=float, required=False, default=config.cluster_bins.diploidbaf, help="Maximum diploid-BAF shift used to determine the largest copy-neutral cluster and to rescale all the cluster inside this threshold accordingly (default: None, scaling is not performed)")
    parser.add_argument("-tR","--tolerancerdr", type=float, required=False, default=config.cluster_bins.tolerancerdr, help='Refine the clustering merging the clusters with this maximum difference in RDR values (default: None, gurobipy required)')
    parser.add_argument("-tB","--tolerancebaf", type=float, required=False, default=config.cluster_bins.tolerancebaf, help='Refine the clustering merging the clusters with this maximum difference in BAF values (default: None, gurobipy required)')
    parser.add_argument("-u","--bootclustering", type=int, required=False, default=config.cluster_bins.bootclustering, help='Number of points to add for bootstraping each bin to improve the clustering. Each point is generated by drawing its values from a normal distribution centered on the values of the bin. This can help the clustering when the input number of bins is low (default: 0)')
    parser.add_argument("-dR","--ratiodeviation", type=float, required=False, default=config.cluster_bins.ratiodeviation, help='Standard deviation of the read ratios used to generate the points in the clouds (default: 0.02)')
    parser.add_argument("-dB","--bafdeviation", type=float, required=False, default=config.cluster_bins.bafdeviation, help='Standard deviation of the BAFs used to generate the points in the clouds (default: 0.02)')
    parser.add_argument("-e","--seed", type=int, required=False, default=config.cluster_bins.seed, help='Random seed used for clustering AND the normal distributions used in the clouds (default: 0)')
    parser.add_argument("-K","--initclusters", type=int, required=False, default=config.cluster_bins.initclusters, help="The maximum number of clusters to infer (default: 50)")
    parser.add_argument("-c","--concentration", type=float, required=False, default=config.cluster_bins.concentration, help="Tuning parameter for clustering (concentration parameter for Dirichlet process prior). Higher favors more clusters, lower favors fewer clusters (default 0.02 = 1/K).")
    parser.add_argument("-R","--restarts", type=int, required=False, default=config.cluster_bins.restarts, help="Number of restarts performed by the clustering to choose the best (default: 10)")
    parser.add_argument("-v","--verbose", action='store_true', default=config.cluster_bins.verbose, required=False, help="Use verbose log messages")
    parser.add_argument("--disablebar", action='store_true', default=config.cluster_bins.disablebar, required=False, help="Disable progress bar")
    parser.add_argument("-V", "--version", action='version', version=f'%(prog)s {__version__}')
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
    if args.concentration < 0:
        raise ValueError(sp.error("Concentration parameter must be positive!"))
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
            "concentration" : args.concentration,
            "restarts" : args.restarts,
            "verbose" : args.verbose,
            "disable" : args.disablebar,
            "outsegments" : args.outsegments,
            "outbins" : args.outbins}


def parse_plot_bins_args(args=None):
    """
    Parse command line arguments
    Returns:
    """
    description = "Generate plots for read-depth ratio (RD), B-allele frequency (BAF), and clusters for genomic bins in multiple samples using .bb, .cbb, .seg files."
    parser = argparse.ArgumentParser(prog='hatchet plot-bins', description=description, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("INPUT", help='Input BBC file with RDR and BAF')
    parser.add_argument("-c", "--command", required=False, type=str, default=config.plot_bins.command, help="The command determining the plots to generate (default: all)\n\n\t{}RD{}: Plot the read-depth ratio (RD) values of the genomes for each sample.\n\n\t{}CRD{}: Plot the read-depth ratio (CRD) values of the genomes for each sample colored by corresponding cluster.\n\n\t{}BAF{}: Plot the B-allele frequency (BAF) values of the genomes for each sample.\n\n\t{}CBAF{}: Plot BAF values for each sample colored by corresponding cluster.\n\n\t{}BB{}: Plot jointly the values of read-depth ratio (RD) and B-allele frequency (BAF) for each bin in all samples and their density.\n\n\t{}CBB{}: Plot jointly the values of read-depth ratio (RD) and B-allele frequency (BAF) for each bin in all samples by coloring the bins depending on their cluster.\n\n\t{}CLUSTER{}: Plot jointly the values of read-depth ratio (RD) and B-allele frequency (BAF) for each cluster in all samples where the size of the markers is proportional to the number of bins.".format(sp.bcolors.BOLD, sp.bcolors.ENDC, sp.bcolors.BOLD, sp.bcolors.ENDC, sp.bcolors.BOLD, sp.bcolors.ENDC, sp.bcolors.BOLD, sp.bcolors.ENDC, sp.bcolors.BOLD, sp.bcolors.ENDC, sp.bcolors.BOLD, sp.bcolors.ENDC, sp.bcolors.BOLD, sp.bcolors.ENDC, sp.bcolors.BOLD, sp.bcolors.ENDC))
    parser.add_argument("-s", "--segfile", required=False, type=str, default=config.plot_bins.segfile, help="When the corresponding seg file is provided the clusters are also plotted (default: none)")
    parser.add_argument("-m","--colormap", required=False, type=str, default=config.plot_bins.colormap, help='Colormap to use for the colors in the plots, the available colormaps are the following {Set1, Set2, Paired, Dark2, tab10, tab20}')
    parser.add_argument("-tC","--chrthreshold", required=False, type=int, default=config.plot_bins.chrthreshold, help='Only covering at least this number of chromosomes are considered (default: None)')
    parser.add_argument("-tS","--sizethreshold", required=False, type=float, default=config.plot_bins.sizethreshold, help='Only covering at least this genome proportion (default: None)')
    parser.add_argument("--resolution", required=False, default=config.plot_bins.resolution, type=int, help='Resolution of bins (default: bins are not merged)')
    parser.add_argument("--xmin", required=False, default=config.plot_bins.xmin, type=float, help='Minimum value on x-axis for supported plots (default: inferred from data)')
    parser.add_argument("--xmax", required=False, default=config.plot_bins.xmax, type=float, help='Maximum value on x-axis for supported plots (default: inferred from data)')
    parser.add_argument("--ymin", required=False, default=config.plot_bins.ymin, type=float, help='Minimum value on y-axis for supported plots (default: inferred from data)')
    parser.add_argument("--ymax", required=False, default=config.plot_bins.ymax, type=float, help='Maximum value on y-axis for supported plots (default: inferred from data)')
    parser.add_argument("--figsize", required=False, default=config.plot_bins.figsize, type=str, help='Size of the plotted figures in the form "(X-SIZE, Y-SIZE)"')
    parser.add_argument("--markersize", required=False, default=config.plot_bins.markersize, type=int, help='Size of the markers (default: values inferred for each plot)')
    parser.add_argument("--colwrap", required=False, default=config.plot_bins.colwrap, type=int, help='Wrapping the plots in this number of columnes (default: 2)')
    parser.add_argument("--fontscale", required=False, default=config.plot_bins.fontscale, type=float, help='Font scale (default: 1)')
    parser.add_argument("-x","--rundir", required=False, default=config.plot_bins.rundir, type=str, help='Running dirrectory where output the results (default: current directory)')
    parser.add_argument("--pdf", action='store_true', default=config.plot_bins.pdf, required=False, help="Output the bb_clustered figure in PDF format (default: PNG)")
    parser.add_argument("--dpi", required=False, default=config.plot_bins.dpi, type=int, help='DPI of PNG images (default: 900)')
    parser.add_argument("-V", "--version", action='version', version=f'%(prog)s {__version__}')
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
        stdout, stderr = subprocess.Popen("grep -e \"^>\" {}".format(reference), stdout=subprocess.PIPE, shell=True, universal_newlines=True).communicate()
        if stderr is not None:
            raise ValueError("Error in reading the reference: {}".format(reference))
        else:
            ref = set(c[1:].strip().split()[0] for c in stdout.strip().split('\n'))
        if not(chromosomes <= ref):
            raise ValueError("The given reference cannot be used because the chromosome names are inconsistent!\nChromosomes found in BAF files: {}\nChromosomes with the same name found in reference genome: {}".format(chromosomes, ref))

    return sorted(list(chromosomes), key=sp.numericOrder)


def getSQNames(samtools, bamfile):
    header, stderr = subprocess.Popen([samtools, "view", "-H", bamfile], stdout=subprocess.PIPE, shell=False, universal_newlines=True).communicate()
    if stderr is not None:
        raise ValueError("The header of the normal-sample BAM cannot be read with samtools!")
    names = set()
    for line in header.strip().split('\n'):
        line = line.split()
        if line[0] == '@SQ':
            names.add(line[1].split(':')[1])
    return names


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
        if len(res[key]) == 0:
            raise ValueError(sp.error("The regions of chromosome {} could not be determined.".format(key)))
        res[key].sort(key=lambda x : x[0])
        if not all(a[0]<=a[1] and a[1]<=b[0] and b[0]<=b[1] for a, b in zip(res[key], res[key][1:])):
            raise ValueError(sp.error("The regions provided for chromosome {} are non-disjoint or a region start is greater than corresponding region end".format(key)))

    return res
