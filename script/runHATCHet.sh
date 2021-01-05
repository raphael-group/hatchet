#!/usr/bin/env bash

####################################################################################
# Please set up the correct configuration values here below before running HATCHet #
####################################################################################

REF="/path/to/reference.fa" #Please make sure to have produced the reference dictionary /path/to/reference.dict
SAM="/path/to/samtools-home/bin/" #Uncomment if samtools is already in PATH
BCF="/path/to/bcftools-home/bin/" #Uncomment if bcftools is already in PATH
XDIR="/path/to/running-dir/"
NORMAL="/path/to/matched-normal.bam"
BAMS="/path/to/tumor-sample1.bam /path/to/tumor-sample2.bam"
NAMES="Primary Met" #Use the same order as the related tumor BAM files in BAMS above
J=$(python -c 'import multiprocessing as mp; print mp.cpu_count()') #Replace with fixed number if you do not want to use all available cpus 
MINREADS=8 #Use 8 for WGS with >30x and 20 for WES with ~100x
MAXREADS=300 #Use 300 for WGS with >30x and Use 1000 for WES with ~100x
LIST=""
#We reccommend to provide a list of known SNPs, please uncomment the appropriate one below according to the provided reference genome
#If a list is not provided, all genomic positions will be genotyped by BCFtools and will be considered for selection of heterozygous SNPs
#Uncomment the following for reference genome hg19 with `chr` notation
#LIST="https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/GATK/00-All.vcf.gz"
#Uncomment the following for reference genome hg19 without `chr` notation
#LIST="https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/00-All.vcf.gz"
#Uncomment the following for reference genome hg38 with `chr` notation
#LIST="https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/GATK/00-All.vcf.gz"
#Uncomment the following for reference genome hg38 without `chr` notation
#LIST="https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/00-All.vcf.gz"

####################################################################################


##################################################################
# For default run please run the following without changes       #
# Otherwise please follow the related HATCHet's reccommendations #
##################################################################
set -e
set -o xtrace
PS4='\''[\t]'\'
ALLNAMES="Normal ${NAMES}"
export PATH=$PATH:${SAM}
export PATH=$PATH:${BCF}
#source /path/to/virtualenv-python2.7/bin/activate

RDR=${XDIR}rdr/
mkdir -p ${RDR}
SNP=${XDIR}snps/
mkdir -p ${SNP}
BAF=${XDIR}baf/
mkdir -p ${BAF}
BB=${XDIR}bb/
mkdir -p ${BB}
BBC=${XDIR}bbc/
mkdir -p ${BBC}
ANA=${XDIR}analysis/
mkdir -p ${ANA}
RES=${XDIR}results/
mkdir -p ${RES}
EVA=${XDIR}evaluation/
mkdir -p ${EVA}

cd ${XDIR}

python2 -m hatchet binBAM -N ${NORMAL} -T ${BAMS} -S ${ALLNAMES} -b 50kb -g ${REF} -j ${J} -O ${RDR}normal.rdr -o ${BIN}tumor.rdr -v &> ${BIN}bins.log

python2 -m hatchet SNPCaller -N ${NORMAL} -r ${REF} -j ${J} -c ${MINREADS} -C ${MAXREADS} -R ${LIST} -o ${SNP} -v &> ${BAF}bafs.log

python2 -m hatchet deBAF -N ${NORMAL} -T ${BAMS} -S ${ALLNAMES} -r ${REF} -j ${J} -c ${MINREADS} -C ${MAXREADS} -L ${SNP} -O ${BAF}normal.baf -o ${BAF}tumor.baf -v &> ${BAF}bafs.log

python2 -m hatchet comBBo -c ${RDR}normal.rdr -C ${RDR}tumor.rdr -B ${BAF}tumor.baf -e ${RANDOM} > ${BB}bulk.bb

python2 -m hatchet cluBB ${BB}bulk.bb -o ${BBC}bulk.seg -O ${BBC}bulk.bbc -e ${RANDOM} -tB 0.04 -tR 0.15 -d 0.08

cd ${ANA}
python2 -m hatchet BBot -c RD --figsize 6,3 ${BBC}bulk.bbc &
python2 -m hatchet BBot -c CRD --figsize 6,3 ${BBC}bulk.bbc &
python2 -m hatchet BBot -c BAF --figsize 6,3 ${BBC}bulk.bbc &
python2 -m hatchet BBot -c BB ${BBC}bulk.bbc &
python2 -m hatchet BBot -c CBB ${BBC}bulk.bbc -tS 0.01 &
wait

cd ${RES}
python2 -m hatchet solve -i ${BBC}bulk -n2,8 -p 400 -u 0.03 -r ${RANDOM} -j ${J} -eD 6 -eT 12 -g 0.35 -l 0.6 &> >(tee >(grep -v Progress > hatchet.log))

## Increase -l to 0.6 to decrease the sensitivity in high-variance or noisy samples, and decrease it to -l 0.3 in low-variance samples to increase the sensitivity and explore multiple solutions with more clones.
## Increase -u if solutions have clone proportions equal to the minimum threshold -u
## Decrease the number of restarts to 200 or 100 for fast runs, as well as user can decrease the number of clones to -n 2,6 when appropriate or when previous runs suggest fewer clones.
## Increase the single-clone confidence to `-c 0.6` to increase the confidence in the presence of a single tumor clone and further increase this value when interested in a single clone.

cd ${EVA}
python -m hatchet BBeval ${RES}/best.bbc.ucn

