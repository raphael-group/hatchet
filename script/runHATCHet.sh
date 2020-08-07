#!/usr/bin/env bash

REF="/media/vineetb/t5-vineetb/raphael-group/hg19.fa"
SAM="/opt/raphael-group/samtools/bin/"
BCF="/opt/raphael-group/bcftools/bin/"
BNPY="/home/vineetb/git_checkouts/bnpy"

HATCHET_HOME="/home/vineetb/git_checkouts/hatchet"
HATCHET="${HATCHET_HOME}/src/hatchet/bin/HATCHet.py"
UTILS="${HATCHET_HOME}/src/hatchet/utils/"
SOLVER="${HATCHET_HOME}/src/hatchet/solve"

XDIR="/media/vineetb/t5-vineetb/raphael-group/hatchet/out/"
NORMAL="/media/vineetb/t5-vineetb/raphael-group/data/hatchet/SRR5906250.sorted.bam"
BAMS="/media/vineetb/t5-vineetb/raphael-group/data/hatchet/SRR5906251.sorted.bam /media/vineetb/t5-vineetb/raphael-group/data/hatchet/SRR5906253.sorted.bam"
ALLNAMES="Normal TumorOP Tumor2"
NAMES="TumorOP Tumor2"
J=8

set -e
set -o xtrace
PS4='\''[\t]'\'
export PATH=$PATH:${SAM}
export PATH=$PATH:${BCF}
#source /path/to/virtualenv-python2.7/bin/activate

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

\time -v python2 ${UTILS}binBAM.py -st ${SAM} -N ${NORMAL} -T ${BAMS} -S ${ALLNAMES} -b 50kb -g ${REF} -j ${J} -q 11 -O ${BIN}normal.bin -o ${BIN}bulk.bin -v &> ${BIN}bins.log

\time -v python2 ${UTILS}deBAF.py -bt ${BCF} -st ${SAM} -N ${NORMAL} -T ${BAMS} -S ${ALLNAMES} -r ${REF} -j ${J} -q 11 -Q 11 -U 11 -c 8 -C 300 -O ${BAF}normal.baf -o ${BAF}bulk.baf -v &> ${BAF}bafs.log

\time -v python2 ${UTILS}comBBo.py -c ${BIN}normal.bin -C ${BIN}bulk.bin -B ${BAF}bulk.baf -m MIRROR -e 12 > ${BB}bulk.bb

\time -v python2 ${UTILS}cluBB.py ${BB}bulk.bb -by ${BNPY} -o ${BBC}bulk.seg -O ${BBC}bulk.bbc -e ${RANDOM} -tB 0.04 -tR 0.15 -d 0.08

cd ${ANA}
\time -v python2 ${UTILS}BBot.py -c RD --figsize 6,3 ${BBC}bulk.bbc &
\time -v python2 ${UTILS}BBot.py -c CRD --figsize 6,3 ${BBC}bulk.bbc &
\time -v python2 ${UTILS}BBot.py -c BAF --figsize 6,3 ${BBC}bulk.bbc &
\time -v python2 ${UTILS}BBot.py -c BB ${BBC}bulk.bbc &
\time -v python2 ${UTILS}BBot.py -c CBB ${BBC}bulk.bbc -tS 0.01 &
wait

cd ${RES}
\time -v python2 ${HATCHET} ${SOLVER} -i ${BBC}bulk -n2,8 -p 400 -v 3 -u 0.03 -r ${RANDOM} -j ${J} -eD 6 -eT 12 -g 0.35 -l 0.6 &> >(tee >(grep -v Progress > hatchet.log))

## Increase -l to 0.6 to decrease the sensitivity in high-variance or noisy samples, and decrease it to -l 0.3 in low-variance samples to increase the sensitivity and explore multiple solutions with more clones.
## Increase -u if solutions have clone proportions equal to the minimum threshold -u
## Decrease the number of restarts to 200 or 100 for fast runs, as well as user can decrease the number of clones to -n 2,6 when appropriate or when previous runs suggest fewer clones.
## Increase the single-clone confidence to `-c 0.6` to increase the confidence in the presence of a single tumor clone and further increase this value when interested in a single clone.

cd ${EVA}
\time -v python ${UTILS}BBeval.py ${RES}/best.bbc.ucn

