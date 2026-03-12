# HATCHet3 (Under development)

HATCHet3 infers subclonal and allele-specific integer copy number variations and clone proportions using multi-sample bulk sequencing data from various assay types: short-read/long-read Whole-genome/exome DNA sequencing data.

### Quick Start
```bash
mamba env create -f ./environment.yaml -p /path/to/envs/hatchet_env
conda activate /path/to/envs/hatchet_env
pip install -e .

snakemake -p --cores 4 -s ./Snakefile \
    --configfile config/config.yaml \
    --directory <output_dir>
```
