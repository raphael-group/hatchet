# HATCHet3 (Under development)

HATCHet3 infers clone-specific and allele-specific copy number aberrations from multi-sample bulk DNA sequencing data includes short-read WES/WES and long-read.

### Quick Start
```bash
mamba env create -f ./environment.yaml -p /path/to/envs/hatchet_env
conda activate /path/to/envs/hatchet_env
pip install -e .

# edit config/hatchet.yaml (bb_dir, genome_size, region_bed), then:
snakemake -p --cores 4 -s ./Snakefile \
    --configfile config/hatchet.yaml \
    --directory <output_dir>
```
