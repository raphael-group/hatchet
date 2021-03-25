#!/usr/bin/env bash

####################################################################################
# Please set up the correct configuration values here below before running HATCHet #
####################################################################################

REF="/path/to/reference.fa" #Please make sure to have produced the reference dictionary /path/to/reference.dict
LIST=""             # If HATCHet has internet, it selects a list of known germline SNPs based on REF_VERS and CHR_NOTATION below. If not, please provide full path to a locally stored list (.vcf.gz) here.
REF_VERS=""         # Reference version used to select list of known germline SNPs; possible values are "hg19" or "hg38", or leave blank "" if you wish for all positions to be genotyped by bcftools
CHR_NOTATION=true   # Does your reference name chromosomes with "chr" prefix?; possible values true/false
# SAM="/path/to/samtools-home/bin/" #Uncomment if samtools is not already in PATH
# BCF="/path/to/bcftools-home/bin/" #Uncomment if bcftools is not already in PATH
XDIR="/path/to/running-dir/"      #Path for output
NORMAL="/path/to/matched-normal.bam"
BAMS="/path/to/tumor-sample1.bam /path/to/tumor-sample2.bam"
NAMES="Primary Met" #Use the same order as the related tumor BAM files in BAMS above
J=$(python -c 'import multiprocessing as mp; print(mp.cpu_count())') #Replace with fixed number if you do not want to use all available cpus
MINREADS=8 #Use 8 for WGS with >30x and 20 for WES with ~100x
MAXREADS=300 #Use 300 for WGS with >30x and Use 1000 for WES with ~100x
BIN="50kb"   #Bin size for calculating RDR and BAF

################################################################################################################################
# To run HATCHet with phasing please do the following:                                                                         #
# 1. Use a phasing algorithm with the SNP VCF files generated in ${SNP}*.vcf.gz (snps folder by default)                       #
# 2. Combine the phased SNPs for all chromosomes in a unique phased file with `CHROM  POS  PHASE` where:                       #
#      - CHROM is the chromosome of the SNP;                                                                                   #
#      - POS is the genomic position of the SNP;                                                                               #
#      - PHASE is any string that contains 0|1 and 1|0 (lines without those will be excluded as well as those starting with #) #
# 3. Provide the path to the phased file in the variable PHASE here below                                                      #
# 4. Choose haplotype block size BLOCK, 50kb is used by default
# Note: a phased VCF file (with phased genotypes 0|1 and 1|0) works and `bcftools concat` can be used to combine chromosomes   #                                         #
# If using reference-phasing algorithm please make sure the ouput VCF are w.r.t. same reference genome, otherwise please       #
# use LiftOver to convert it or bcftools --annotate to add or remove `chr` notation                                            #
################################################################################################################################
PHASE="None"  #Path to phased file; specify "None" to run hatchet without phasing
BLOCK="50kb"  #Haplotype block size  used for combining SNPs


# These specify the subdirectories created and used by HATCHet and do not need to be changed

RDR="rdr/"
SNP="snps/"
BAF="baf/"
BB="bb/"
BBC="bbc/"
PLO="plots/"
RES="results/"
SUM="summary/"


