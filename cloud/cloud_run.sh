#!/bin/bash

dsub \
  --provider google-v2 \
  --project PROJECT_ID \
  --regions us-east1 \
  --logging "gs://BUCKET_NAME/RUN_NAME/logging" \
  --disk-size 1000 \
  --name "RUN_NAME" \
  --image gcr.io/durable-tracer-294016/hatchet \
  --input NORMALBAM="gs://gdc-tcga-phs000178-controlled/BRCA/DNA/WGS/WUGSC/ILLUMINA/b9774dd35c320f70de8f2b81c15d5a98.bam" \
  --input NORMALBAI="gs://gdc-tcga-phs000178-controlled/BRCA/DNA/WGS/WUGSC/ILLUMINA/b9774dd35c320f70de8f2b81c15d5a98.bam.bai" \
  --input TUMORBAM1="gs://gdc-tcga-phs000178-controlled/BRCA/DNA/WGS/WUGSC/ILLUMINA/2258e57e8e0af9db6969a1da86177ca7.bam" \
  --input TUMORBAI1="gs://gdc-tcga-phs000178-controlled/BRCA/DNA/WGS/WUGSC/ILLUMINA/2258e57e8e0af9db6969a1da86177ca7.bam.bai" \
  --output-recursive OUTPUT_FOLDER="gs://BUCKET_NAME/RUN_NAME/output" \
  --script "_run.sh" \
  --wait
