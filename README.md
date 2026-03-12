# HATCHet3 (Under development)

HATCHet3 infers subclonal and allele-specific integer copy number variations and clone proportions using multi-sample bulk sequencing data from various assay types: short-read/long-read Whole-genome/exome DNA sequencing data.

### Installation

```bash
mamba env create -f ./environment.yaml -p /path/to/envs/hatchet_env
conda activate /path/to/envs/hatchet_env
pip install -e .
```

### Quick Start (Snakemake)

```bash
snakemake -p --cores 4 -s ./Snakefile \
    --configfile config/config.yaml \
    --directory <output_dir>
```

The pipeline runs three steps in order:

1. **cluster-bins** -- clusters genomic bins via HMM (output: `bbc/`)
2. **compute-cn** -- infers copy number states via ILP/CD solver (output: `results/`)
3. **plot-cn** -- generates 1D CN profiles and 2D RDR-vs-BAF scatter plots (output: `summary/`)

#### Key config options (`config/hatchet.yaml`)

| Section | Parameter | Default | Description |
|---|---|---|---|
| top-level | `bb_dir` | -- | Input directory with NPZ count matrices and `bb.tsv.gz` |
| top-level | `genome_size` | -- | Reference chromosome sizes file |
| top-level | `region_bed` | -- | Reference chromosome BED file |
| `cluster_bins` | `minK` / `maxK` | 3 / 30 | Range of HMM cluster states to search |
| `compute_cn` | `mode` | `ilp` | Solver mode: `ilp`, `cd`, or `both` |
| `compute_cn` | `solver` | `gurobi` | ILP solver: `gurobi` or `cbc` |
| `compute_cn` | `minClone` / `maxClone` | 2 / 4 | Range of tumor clones to solve for |
| `compute_cn` | `diploid` / `tetraploid` | true / true | Enable diploid and/or tetraploid modes |

If you use `mode=ilp` during `compute_cn` step, be awared that runtime is slow when `maxK>3`. You may also look at `<bbc>/plots/` to decide a smaller number of clusters if BIC/ICL model-selection picks too many states.

#### Running individual steps

```bash
# Step 1: cluster bins
hatchet cluster-bins --bb_dir <bb_dir> --bbc_dir <bbc_dir> --genome_size <genome.sizes>

# Step 2: compute copy numbers
hatchet compute-cn --bbc <bbc_dir/bulk.bbc> --seg <bbc_dir/bulk.seg> \
    --result_dir <results> --genome_size <genome.sizes> --region_bed <regions.bed>

# Step 3: plot results
hatchet plot-cn --bbc <results/best.bbc.ucn> --seg <results/best.seg.ucn> \
    --plot_dir <summary> --gamma_file <results/gammas.tsv> \
    --genome_size <genome.sizes> --region_bed <regions.bed> --ploidy diploid
```

Use `hatchet <command> --help` for full argument details.
