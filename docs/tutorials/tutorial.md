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
- `genome_size`: path to chromosome sizes file (e.g., `data/reference/hg38.chrom.sizes`)
- `region_bed`: path to chromosome segments BED file (e.g., `data/reference/hg38.regions.bed`)

### 2. Run with Snakemake

```bash
snakemake -p --cores 4 -s ./Snakefile \
    --configfile config/my_config.yaml \
    --directory output/my_sample
```

This runs three steps sequentially:

1. **cluster-bins**: clusters genomic bins using an HMM
2. **compute-cn**: infers integer copy number states via ILP/CD optimization
3. **plot-cn**: generates 1D CN profiles and 2D RDR-vs-BAF scatter plots for model-selected final result.

See [Output Structure](../reference.md#output-structure) for the full directory layout.

## Tips

- **Solver choice**: Gurobi is faster but requires a license. CBC (`--solver cbc`) is open-source and included in the conda environment.
- **Cluster count range**: If you know the approximate number of CN segments, narrow the `minK`/`maxK` range to speed up cluster-bins.
- **Manual K selection**: If cluster-bins picks a suboptimal K, use `compute_cn.k` in the config (or run compute-cn with a specific `bbc/labels/bulk<K>.bbc` file) to override.
- **Coordinate descent**: Use `--mode cd` or `--mode both` for large instances where ILP is slow. CD uses random restarts (`--cd_nseeds 400`) and can run in parallel (`--cd_njobs 8`).

Use `hatchet <command> --help` for full argument details.
