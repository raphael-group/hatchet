#!/usr/bin/env bash

set -e
set -o xtrace

# -----------------------------------
# Reference Genome
# -----------------------------------
# Uncomment ONE of the REFERENCE lines below
# Links obtained from
#   https://gatk.broadinstitute.org/hc/en-us/articles/360035890711?id=23390
#   https://cloud.google.com/life-sciences/docs/resources/public-datasets/reference-genomes
# GRCh37
# REFERENCE=ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_genomic.fna.gz
# GRCh37lite
# REFERENCE=https://www.bcgsc.ca/downloads/genomes/9606/hg19/1000genomes/bwa_ind/genome/GRCh37-lite.fa
# hg19
REFERENCE=http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
# b37
# REFERENCE=ftp://ftp.broadinstitute.org/pub/seq/references/Homo_sapiens_assembly19.fasta
# hs37d5
# REFERENCE=ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
# HumanG1Kv37
# REFERENCE=http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz
# GRCh38
# REFERENCE=http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
# GRCh38_Verily_v1
# REFERENCE=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz

REFERENCE_FILENAME=$(basename -- "$REFERENCE")
wget $REFERENCE -P /mnt/data/
if [[ $REFERENCE_FILENAME == *.gz ]]
then
    gunzip /mnt/data/$REFERENCE_FILENAME
    REFERENCE_FILENAME="${REFERENCE_FILENAME%.*}"
fi
REFERENCE_NAME="${REFERENCE_FILENAME%.*}"

$HATCHET_PATHS_SAMTOOLS/samtools faidx /mnt/data/$REFERENCE_FILENAME
$HATCHET_PATHS_SAMTOOLS/samtools dict /mnt/data/$REFERENCE_FILENAME > /mnt/data/$REFERENCE_NAME.dict
export HATCHET_PATHS_REFERENCE=/mnt/data/$REFERENCE_FILENAME
# -----------------------------------

TUMOR_BAMS=""
ALLNAMES="NORMAL"
for BAM_N in "${!TUMORBAM@}"; do
  TUMOR_BAM="${!BAM_N}"
  TUMOR_BAMS+=" ${TUMOR_BAM}"
  TUMOR_BAM_NAME="$(basename "${TUMOR_BAM}")"
  TUMOR_BAM_NAME="${TUMOR_BAM_NAME%.*}"
  ALLNAMES+=" ${TUMOR_BAM_NAME}"
done

BIN=${OUTPUT_FOLDER}/bin/
BAF=${OUTPUT_FOLDER}/baf/
BB=${OUTPUT_FOLDER}/bb/
BBC=${OUTPUT_FOLDER}/bbc/
ANA=${OUTPUT_FOLDER}/analysis/
RES=${OUTPUT_FOLDER}/results/
EVA=${OUTPUT_FOLDER}/evaluation/

mkdir -p ${BIN} ${BAF} ${BB} ${BBC} ${ANA} ${RES} ${EVA}

python -m hatchet count-reads -N ${NORMALBAM} -T ${TUMOR_BAMS} -S ${ALLNAMES} -b 50kb -O ${BIN}normal.1bed -o ${BIN}bulk.1bed &> ${BIN}bins.log
python -m hatchet count-alleles -N ${NORMALBAM} -T ${TUMOR_BAMS} -S ${ALLNAMES} -O ${BAF}normal.1bed -o ${BAF}bulk.1bed &> ${BAF}bafs.log
python -m hatchet combine-counts -c ${BIN}normal.1bed -C ${BIN}bulk.1bed -B ${BAF}bulk.1bed > ${BB}bulk.bb
python -m hatchet cluster-bins ${BB}bulk.bb -o ${BBC}bulk.seg -O ${BBC}bulk.bbc

cd ${ANA}
python -m hatchet plot-bins -c RD ${BBC}bulk.bbc
python -m hatchet plot-bins -c CRD ${BBC}bulk.bbc
python -m hatchet plot-bins -c BAF ${BBC}bulk.bbc
python -m hatchet plot-bins -c BB ${BBC}bulk.bbc
python -m hatchet plot-bins -c CBB ${BBC}bulk.bbc

# ------------------------------------------------------
# Commented out till the solver works inside a container
# ------------------------------------------------------
# cd ${RES}
# python -m hatchet compute-cn -i ${BBC}bulk &> >(tee >(grep -v Progress > hatchet.log))

# cd ${EVA}
# python -m hatchet plot-cn ${RES}/best.bbc.ucn
