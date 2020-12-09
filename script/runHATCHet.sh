#!/usr/bin/env bash

REF="/path/to/reference.fa"
SAM="/path/to/samtools-home/bin/"
BCF="/path/to/bcftools-home/bin/"

XDIR="/path/to/running-dir/"
NORMAL="/path/to/matched-normal.bam"
BAMS="/path/to/tumor-sample1.bam /path/to/tumor-sample2.bam"
ALLNAMES="Normal Primary Met"
NAMES="Primary Met"
J=22

set -e
set -o xtrace
PS4='\''[\t]'\'
export PATH=$PATH:${SAM}
export PATH=$PATH:${BCF}
#source /path/to/virtualenv-python3.8/bin/activate

BIN=${XDIR}bin/
mkdir -p ${BIN}
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

\time -v python3 -m hatchet binBAM -st ${SAM} -N ${NORMAL} -T ${BAMS} -S ${ALLNAMES} -b 50kb -g ${REF} -j ${J} -q 11 -O ${BIN}normal.bin -o ${BIN}bulk.bin -v &> ${BIN}bins.log

\time -v python3 -m hatchet deBAF -bt ${BCF} -st ${SAM} -N ${NORMAL} -T ${BAMS} -S ${ALLNAMES} -r ${REF} -j ${J} -q 11 -Q 11 -U 11 -c 8 -C 300 -O ${BAF}normal.baf -o ${BAF}bulk.baf -v &> ${BAF}bafs.log

\time -v python3 -m hatchet comBBo -c ${BIN}normal.bin -C ${BIN}bulk.bin -B ${BAF}bulk.baf -m MIRROR -e 12 > ${BB}bulk.bb

\time -v python3 -m hatchet cluBB ${BB}bulk.bb -o ${BBC}bulk.seg -O ${BBC}bulk.bbc -e ${RANDOM} -tB 0.04 -tR 0.15 -d 0.08

cd ${ANA}
\time -v python3 -m hatchet BBot -c RD --figsize 6,3 ${BBC}bulk.bbc &
\time -v python3 -m hatchet BBot -c CRD --figsize 6,3 ${BBC}bulk.bbc &
\time -v python3 -m hatchet BBot  -c BAF --figsize 6,3 ${BBC}bulk.bbc &
\time -v python3 -m hatchet BBot  -c BB ${BBC}bulk.bbc &
\time -v python3 -m hatchet BBot  -c CBB ${BBC}bulk.bbc -tS 0.01 &
wait

cd ${RES}
\time -v python3 -m hatchet solve -i ${BBC}bulk -n2,8 -p 400 -v 3 -u 0.03 -r ${RANDOM} -j ${J} -eD 6 -eT 12 -g 0.35 -l 0.6 &> >(tee >(grep -v Progress > hatchet.log))

## Increase -l to 0.6 to decrease the sensitivity in high-variance or noisy samples, and decrease it to -l 0.3 in low-variance samples to increase the sensitivity and explore multiple solutions with more clones.
## Increase -u if solutions have clone proportions equal to the minimum threshold -u
## Decrease the number of restarts to 200 or 100 for fast runs, as well as user can decrease the number of clones to -n 2,6 when appropriate or when previous runs suggest fewer clones.
## Increase the single-clone confidence to `-c 0.6` to increase the confidence in the presence of a single tumor clone and further increase this value when interested in a single clone.

cd ${EVA}
\time -v python -m hatchet BBeval ${RES}/best.bbc.ucn

