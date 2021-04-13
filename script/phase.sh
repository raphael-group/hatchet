#!/usr/bin/env bash

source ./config_lo.sh
#source ./config.sh


cd ${XDIR}
mkdir -p ${PHASE}
python3 -m hatchet Phase -j 22 -R ${REF_PANEL} -V ${REF_VERS} -N ${CHR_NOTATION} -L ${SNP}*.vcf.gz -o ${PHASE}




