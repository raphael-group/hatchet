#!/usr/bin/env bash

source ./config.sh


cd ${XDIR}
mkdir -p ${PHASE}
python3 -m hatchet Phase -j 2 -R ${REF_PANEL} -L ${SNP}*.vcf.gz -o ${PHASE}

