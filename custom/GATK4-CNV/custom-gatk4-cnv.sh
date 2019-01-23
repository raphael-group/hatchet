#!/usr/bin/bash

# This is a custom complete pipeline of HATCHet which considers in input segmented files for one or more samples from the same patient, produced by the GATK4 CNV pipeline.

BNPY="/path/to/bnpy-dev/" # Provide the full path to BNPY's repository

HATCHET_HOME="/path/to/hatchet_home" # Provide the full path to HATCHet's repository

HATCHET="${HATCHET_HOME}/bin/hatchet.py"
UTILS="${HATCHET_HOME}/utils/"
SOLVER="${HATCHET_HOME}/build/solve"
CNVTOBB="${HATCHET_HOME}/custom/GATK4-CNV/gatk4cnsToBB.py"

XDIR="/path/to/running-dir/" # Specify running directory where all results will be written
SEGS="/path/to/tumor-sample1.seg /path/to/tumor-sample2.seg" # Specify the full path to segmented files
NAMES="Sample1 Sample2" # Specify the names of samples in the same order as above
J=22 # Specify the number of threads that can be used

set -e
set -o xtrace
PS4='\''[\t]'\'
#source /path/to/virtualenv-python2.7/bin/activate

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

\time -v python2 ${CNVTOBB} "${SEGS}" -- samples "${NAMES}" > ${BB}bulk.bb

\time -v python2 ${UTILS}cluBB.py ${BB}bulk.bb -by ${BNPY} -o ${BBC}bulk.seg -O ${BBC}bulk.bbc -e 12 -tB 0.04 -tR 0.15 -d 0.08

cd ${ANA}
\time -v python2 ${UTILS}BBot.py -c RD --figsize 6,3 ${BBC}bulk.bbc &
\time -v python2 ${UTILS}BBot.py -c CRD --figsize 6,3 ${BBC}bulk.bbc &
\time -v python2 ${UTILS}BBot.py -c BAF --figsize 6,3 ${BBC}bulk.bbc &
\time -v python2 ${UTILS}BBot.py -c BB ${BBC}bulk.bbc &
\time -v python2 ${UTILS}BBot.py -c CBB ${BBC}bulk.bbc &
wait

cd ${RES}
\time -v python2 ${HATCHET} ${SOLVER} -i ${BBC}bulk -n2,8 -p 400 -v 3 -u 0.03 -r 12 -j ${J} -eD 6 -eT 12 -g 0.35 -l 0.6 |& tee hatchet.log

cd ${EVA}
\time -v python ${UTILS}BBeval.py ${RES}/best.bbc.ucn

