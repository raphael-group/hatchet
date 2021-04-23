#!/usr/bin/env bash

source ./config.sh

set -e
set -o xtrace
PS4='\''[\t]'\'

mkdir -p ${REF_PANEL_DIR}
python3 -m hatchet PhasePrep -D ${REF_PANEL_DIR} -R ${REF_PANEL} -V ${REF_VERS} -N ${CHR_NOTATION}

