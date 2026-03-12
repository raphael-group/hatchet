# HATCHet3 Tutorial

This tutorial walks through running HATCHet3 on bulk DNA sequencing data.

## Prerequisites

- Conda environment with HATCHet3 installed (see README)
- Input genotyping results from universal-genotyping pipeline
- Reference files: chromosome sizes file and segments BED file

## Input Format

HATCHet3 expects pre-processed genotyping output in `bb_dir/`:

```
bb_dir/
  bb.tsv.gz           # bin metadata (chr, start, end, sample, RDR, BAF, counts)
  bb.rdr.npz          # read-depth ratio matrix (bins x samples)
  bb.Aallele.npz      # A-allele count matrix
  bb.Ballele.npz      # B-allele count matrix
  bb.Tallele.npz      # total allele count matrix
  sample_ids.tsv      # sample ID table with columns: SAMPLE, sample_type
```

## Running the Full Pipeline

### 1. Configure

Copy and edit the config template:

```bash
cp config/hatchet.yaml config/my_config.yaml
```

Key settings to update:
- `bb_dir`: absolute path to your genotyping output directory
- `genome_size`: path to chromosome sizes file (e.g., `hg38.chrom.sizes`)
- `region_bed`: path to chromosome segments BED file (e.g., `hg38.segments.bed`)

### 2. Run with Snakemake

```bash
snakemake -p --cores 4 -s ./Snakefile \
    --configfile config/my_config.yaml \
    --directory output/my_sample
```

This runs three steps sequentially:

1. **cluster-bins**: clusters genomic bins using an HMM
2. **compute-cn**: infers integer copy number states via ILP/CD optimization
3. **plot-cn**: generates 1D CN profiles and 2D RDR-vs-BAF scatter plots

Output structure:

```
output/my_sample/
  bbc/                  # cluster-bins output
    bulk.bbc            # per-bin cluster assignments with allele counts
    bulk.seg            # per-cluster summary statistics
  results/              # compute-cn output
    best.bbc.ucn        # model-selected CN solution (per-bin)
    best.seg.ucn        # model-selected CN solution (per-cluster)
    gammas.tsv          # RDR scaling factors per sample
    elbow_curve.png     # model selection elbow plot
    results.diploid.n*.bbc.ucn.tsv   # per-n solutions
  summary/                # plot-cn output
    diploid_n*/         # per-solution plots
      <sample>.1D.png   # 1D genome-wide CN profile
      <sample>.2D.png   # 2D RDR-vs-BAF scatter plot
```

## Running Individual Steps

### Step 1: Cluster Bins

```bash
hatchet cluster-bins \
    --bb_dir /path/to/bb_dir \
    --bbc_dir output/bbc \
    --genome_size /path/to/hg38.chrom.sizes \
    --minK 3 --maxK 30 \
    --restarts 10 --top_restarts 5 \
    --niters 10 \
    --decode_method viterbi \
    --score_method bic \
    --verbosity 1
```

Key parameters:
- `--minK / --maxK`: range of cluster counts to search (default: 3-30)
- `--score_method`: model selection criterion, `bic` or `icl` (default: bic)
- `--decode_method`: HMM decoding, `viterbi` or `map` (default: viterbi)
- `--restarts`: number of random initializations per K (default: 10)
- `--niters`: EM iterations per restart (default: 10)

### Step 2: Compute Copy Numbers

```bash
hatchet compute-cn \
    --bbc output/bbc/bulk.bbc \
    --seg output/bbc/bulk.seg \
    --result_dir output/results \
    --genome_size /path/to/hg38.chrom.sizes \
    --region_bed /path/to/hg38.segments.bed \
    --mode ilp \
    --solver gurobi \
    --minClone 2 --maxClone 4 \
    --diploid --tetraploid \
    --reg_term MAXCN \
    --verbosity 2
```

Key parameters:
- `--mode`: solver mode — `ilp` (exact), `cd` (coordinate descent), or `both` (default: ilp)
- `--solver`: ILP backend — `gurobi` (requires license) or `cbc` (open-source) (default: gurobi)
- `--minClone / --maxClone`: range of tumor clones to solve for (default: 2-4)
- `--diploid / --tetraploid`: enable diploid (cn_max=6) and/or tetraploid (cn_max=12) modes
- `--reg_term`: regularization — `RAW`, `MAXCN`, `DROOT_SUM`, or `DADJ_SUM` (default: MAXCN)

### Step 3: Plot Results

```bash
hatchet plot-cn \
    --bbc output/results/best.bbc.ucn \
    --seg output/results/best.seg.ucn \
    --plot_dir output/summary \
    --gamma_file output/results/gammas.tsv \
    --genome_size /path/to/hg38.chrom.sizes \
    --region_bed /path/to/hg38.segments.bed \
    --ploidy diploid
```

Key parameters:
- `--ploidy`: `diploid` or `tetraploid` (required; selects the correct gamma scaling)
- `--maxlim_fcn`: upper limit for the FCN y-axis (default: 30)
- `--img_type`: output format — `png`, `pdf`, or `svg` (default: png)

## Tips

- **Solver choice**: Gurobi is faster but requires a license. CBC (`--solver cbc`) is open-source and included in the conda environment.
- **Cluster count range**: If you know the approximate number of CN segments, narrow the `minK`/`maxK` range to speed up cluster-bins.
- **Manual K selection**: If cluster-bins picks a suboptimal K, use `compute_cn.k` in the config (or run compute-cn with a specific `bbc/labels/bulk<K>.bbc` file) to override.
- **Coordinate descent**: Use `--mode cd` or `--mode both` for large instances where ILP is slow. CD uses random restarts (`--cd_nseeds 400`) and can run in parallel (`--cd_njobs 8`).

Use `hatchet <command> --help` for full argument details.
