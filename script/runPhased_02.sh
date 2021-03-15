#!/usr/bin/env bash

source ./config.sh

##################################################################
# For default run please execute the following without changes   #
# Otherwise please follow the related HATCHet's reccommendations #
# To run HATCHet with phasing of SNPs please see below           #
##################################################################
set -e
set -o xtrace
PS4='\''[\t]'\'
ALLNAMES="Normal ${NAMES}"
export PATH=$PATH:${SAM}
export PATH=$PATH:${BCF}
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1

cd ${XDIR}
mkdir -p ${BB}
mkdir -p ${BBC}
mkdir -p ${PLO}
mkdir -p ${RES}
mkdir -p ${SUM}

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

python3 -m hatchet comBBo -c ${RDR}normal.1bed -C ${RDR}tumor.1bed -B ${BAF}tumor.1bed -t ${RDR}total.tsv -p ${PHASE} -l ${BLOCK} -e ${RANDOM} > ${BB}bulk.bb

python3 -m hatchet cluBB ${BB}bulk.bb -o ${BBC}bulk.seg -O ${BBC}bulk.bbc -e ${RANDOM} -tB 0.04 -tR 0.15 -d 0.08

cd ${PLO}
python3 -m hatchet BBot -c RD --figsize 6,3 ../${BBC}bulk.bbc
python3 -m hatchet BBot -c CRD --figsize 6,3 ../${BBC}bulk.bbc
python3 -m hatchet BBot -c BAF --figsize 6,3 ../${BBC}bulk.bbc
python3 -m hatchet BBot -c BB ../${BBC}bulk.bbc
python3 -m hatchet BBot -c CBB ../${BBC}bulk.bbc -tS 0.01

cd ../${RES}
python3 -m hatchet solve -i ../${BBC}bulk -n2,6 -p 400 -u 0.03 -r ${RANDOM} -j ${J} -eD 6 -eT 12 -g 0.35 -l 0.6 &> >(tee >(grep -v Progress > hatchet.log))

## Increase -l to 0.6 to decrease the sensitivity in high-variance or noisy samples, and decrease it to -l 0.3 in low-variance samples to increase the sensitivity and explore multiple solutions with more clones.
## Increase -u if solutions have clone proportions equal to the minimum threshold -u
## Decrease the number of restarts to 200 or 100 for fast runs, as well as user can increase the number of clones to -n 2,8 when appropriate or when previous runs suggest fewer clones (i.e. OBJ function keeps decreasing).
## Increase the single-clone confidence to `-c 0.6` to increase the confidence in the presence of a single tumor clone and further increase this value when interested in a single clone.

cd ../${SUM}
python3 -m hatchet BBeval ../${RES}/best.bbc.ucn

