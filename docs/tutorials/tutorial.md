# HATCHet3 Tutorial

This tutorial walks through running HATCHet3 on bulk DNA sequencing data.

## Prerequisites

- Conda environment with HATCHet3 installed (see README)
- Input genotyping results from the [Universal-Genotyping-Pipeline](https://github.com/raphael-group/Universal-Genotyping-Pipeline)
- Reference files: chromosome sizes file and whitelist regions BED file

## Input Format

HATCHet3 expects pre-processed genotyping output in `bb_dir/`:

```
bb_dir/
  bb.tsv.gz           # bin metadata (chr, start, end, region_id, switchprobs)
  bb.rdr.npz          # read-depth ratio matrix (bins x samples)
  bb.depth.npz        # read-depth matrix (bins x samples)
  bb.Aallele.npz      # A-allele count matrix
  bb.Ballele.npz      # B-allele count matrix
  bb.Tallele.npz      # total allele count matrix
  sample_ids.tsv      # sample ID table with a SAMPLE column
```

## Running the Full Pipeline

### 1. Configure

Copy and edit the config template:

```bash
cp config/snakemake-hatchet.yaml config/my_config.yaml
```

Key settings to update:
- `bb_dir`: absolute path to your genotyping output directory
- `genome_size`: path to chromosome sizes file (e.g., `data/reference/hg38.chrom.sizes`)
- `region_bed`: path to whitelist regions BED file (e.g., `data/reference/hg38.regions.bed`)
- `patient_id`: output filename prefix for the per-solution plots

Per-stage options live under the `cluster-bins:` and `compute-cn:` blocks; the Snakefile
forwards only valid flags to each subcommand (see [Reference](../reference.md#parameter-reference)).

### 2. Run with Snakemake

```bash
snakemake -p --cores 4 -s ./Snakefile \
    --configfile config/my_config.yaml \
    --directory output/my_sample
```

The workflow runs two rules in sequence:

1. **cluster-bins**: clusters genomic bins with a phase-aware Gaussian + Beta-Binomial
   factorial HMM, sweeping K and selecting by ICL/BIC.
2. **compute-cn**: infers integer copy numbers and clone proportions via ILP/CD, sweeps a
   regularization path, selects clone number and ploidy, and emits its own CN plots under
   `results/plots/`.

The final targets are `results/best.bbc.ucn` and `results/best.seg.ucn`. Per-rule wall time,
CPU, and peak memory are recorded in `logs/*.benchmark.tsv`. See
[Output Structure](../reference.md#output-structure) for the full directory layout.

### 3. (Optional) Regenerate plots

`compute-cn` already produces plots. To regenerate or customize a specific solution, run the
standalone command:

```bash
hatchet plot-cn \
    --bbc output/my_sample/results/best.bbc.ucn \
    --seg output/my_sample/results/best.seg.ucn \
    --gamma_file output/my_sample/results/gammas.tsv \
    --ploidy diploid \
    --genome_size data/reference/hg38.chrom.sizes \
    --region_bed data/reference/hg38.regions.bed \
    --img_type svg \
    -O output/my_sample/plots_custom
```

## Tips

- **Solver choice**: Gurobi is faster but requires a license. CBC (`--solver cbc`) is open-source and included in the conda environment.
- **Cluster count range**: If you know the approximate number of CN segments, narrow the `minK`/`maxK` range to speed up cluster-bins.
- **Manual K selection**: If cluster-bins picks a suboptimal K, set `compute-cn.k` in the config (or point compute-cn at `bbc/labels/bulk<K>.bbc` / `bulk<K>.seg`) to override.
- **Coordinate descent**: Use `--mode cd` or `--mode both` for large instances where ILP is slow. CD uses random restarts (`--cd_nseeds 400`) and runs in parallel (`--cd_njobs 8`).

Use `hatchet <command> --help` for full argument details.
