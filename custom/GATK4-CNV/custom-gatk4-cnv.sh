#!/usr/bin/bash

# This is a custom complete pipeline of HATCHet2 which considers in input segmented files for one or more samples from the same patient, produced by the GATK4 CNV pipeline.

HATCHET_HOME="/path/to/hatchet_home" # Provide the full path to HATCHet2's repository

CNVTOBB="${HATCHET_HOME}/custom/GATK4-CNV/gatk4cnsToBB.py"

XDIR="/path/to/running-dir/" # Specify running directory where all results will be written
SEGS="/path/to/tumor-sample1.seg /path/to/tumor-sample2.seg" # Specify the full path to segmented files
NAMES="Sample1 Sample2" # Specify the names of samples in the same order as above
J=22 # Specify the number of threads that can be used

set -e
set -o xtrace
PS4='\''[\t]'\'
#source /path/to/virtualenv-python3.8/bin/activate

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

\time -v python3 ${CNVTOBB} "${SEGS}" -- samples "${NAMES}" > ${BB}bulk.bb

\time -v python3 -m hatchet cluster-bins ${BB}bulk.bb -o ${BBC}bulk.seg -O ${BBC}bulk.bbc -e 12 -tB 0.04 -tR 0.15 -d 0.08

cd ${ANA}
\time -v python3 -m hatchet plot-bins -c RD --figsize 6,3 ${BBC}bulk.bbc &
\time -v python3 -m hatchet plot-bins -c CRD --figsize 6,3 ${BBC}bulk.bbc &
\time -v python3 -m hatchet plot-bins -c BAF --figsize 6,3 ${BBC}bulk.bbc &
\time -v python3 -m hatchet plot-bins -c BB ${BBC}bulk.bbc &
\time -v python3 -m hatchet plot-bins -c CBB ${BBC}bulk.bbc &
wait

cd ${RES}
\time -v python3 -m hatchet solve -i ${BBC}bulk -n2,8 -p 400 -v 3 -u 0.03 -r 12 -j ${J} -eD 6 -eT 12 -g 0.35 -l 0.6 |& tee hatchet.log

cd ${EVA}
\time -v python3 -m hatchet plot-cn ${RES}/best.bbc.ucn
