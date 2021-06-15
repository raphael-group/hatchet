#!/usr/bin/env bash

####################################################################################
# Please set up the correct configuration values here below before running HATCHet #
####################################################################################

REF="/n/fs/ragr-data/datasets/ref-genomes/GRCh37_NCBI/GRCh37.p13.fa" #Please make sure to have produced the reference dictionary /path/to/reference.dict
LIST=""             # If HATCHet has internet, it selects a list of known germline SNPs based on REF_VERS and CHR_NOTATION below. If not, please provide full path to a locally stored list (.vcf.gz) here.
REF_VERS="hg19"         # Reference version used to select list of known germline SNPs; possible values are "hg19" or "hg38", or leave blank "" if you wish for all positions to be genotyped by bcftools
CHR_NOTATION=false   # Does your reference name chromosomes with "chr" prefix?; possible values true/false
# SAM="/path/to/samtools-home/bin/" #Uncomment if samtools is not already in PATH
# BCF="/path/to/bcftools-home/bin/" #Uncomment if bcftools is not already in PATH
XDIR="/opt/ragr/bjarnold/hatchet_shapeit_phasing/hatchet/script/dipg2/SJHGG010325"      #Path for output
NORMAL="/opt/ragr/zaccaria/Gundem_15/A12/b627a3ea-5c38-4cc4-8f5b-75b45995dbe6/95bf9ae68444e5e14976fbbe79dc466a.bam"
BAMS="/opt/ragr/zaccaria/Gundem_15/A12/0e3e44bc-9142-4b02-8d8a-d85504e08ca4/cb5ab4c2dd4074286b87209dfa443944.bam"
NAMES="Primary" #Use the same order as the related tumor BAM files in BAMS above
J=22 #Replace with fixed number if you do not want to use all available cpus
MINREADS=8 #Use 8 for WGS with >30x and 20 for WES with ~100x
MAXREADS=1000 #Use 300 for WGS with >30x and Use 1000 for WES with ~100x
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
PHASE=$XDIR"/phase/phased.vcf.gz"  #Path to phased file; specify "None" to run hatchet without phasing
BLOCK="50kb"  #Haplotype block size  used for combining SNPs
REF_PANEL="1000GP_Phase3" # Currently, only "1000GP_Phase3" is supported 
REF_PANEL_DIR="/opt/ragr/bjarnold/hatchet_shapeit_phasing/hatchet/script/refpanel"



# These specify the subdirectories created and used by HATCHet and do not need to be changed

RDR="rdr/"
SNP="snps/"
BAF="baf/"
BB="bb/"
BBC="bbc/"
PLO="plots/"
RES="results/"
SUM="summary/"


