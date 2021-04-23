#!/usr/bin/env bash

#source ./config_lo.sh
#source ./config_hg38.sh
source ./config.sh

#mkdir -p ${REF_PANEL_DIR}
#python3 -m hatchet PhasePrep -D ${REF_PANEL_DIR} -R ${REF_PANEL} -V ${REF_VERS} -N ${CHR_NOTATION}

cd ${XDIR}
mkdir -p ${PHASE}
python3 -m hatchet Phase -j 22 -g ${REF} -D ${REF_PANEL_DIR} -V ${REF_VERS} -N ${CHR_NOTATION} -L ${SNP}*.vcf.gz -o "phase/" 




