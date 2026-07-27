# `cluster-bins`
`cluster-bins` performs local-global genome segmentation using a Gaussian RDR + Beta-Binomial BAF multi-sample factorial HMM with phase switch correction.

## Input

The preprocessed bin-by-sample matrices in `--bb_dir` (RDR, phased allele counts, bin metadata, and the sample table), plus the reference `--genome_size` and `--region_bed`. See [reference.md#input](../reference.md#input) for file layouts and formats.

## Usage

```console
$ hatchet cluster-bins --help
usage: hatchet cluster-bins [-h] --bb_dir BB_DIR --bbc_dir BBC_DIR [--force]
                            [--minK MINK] [--maxK MAXK] [-t T]
                            [--restarts RESTARTS]
                            [--top_restarts TOP_RESTARTS]
                            [--n_local_trials N_LOCAL_TRIALS]
                            [--niters NITERS] [--decode_method {viterbi,map}]
                            [--score_method {bic,icl}]
                            [--score_criteria SCORE_CRITERIA]
                            [--init_method {cna_plus_plus,kmeans_plus_plus}]
                            [--free_baf_c0]
                            [--training_method {baum_welch,viterbi}]
                            [--min_tau MIN_TAU] [--max_tau MAX_TAU]
                            [--share_tau | --no-share_tau] [--baf_eps BAF_EPS]
                            [--min_covar MIN_COVAR] [--ig_alpha IG_ALPHA]
                            [--tau_iters TAU_ITERS] [--seed SEED] [--log_rdr]
                            --genome_size GENOME_SIZE --region_bed REGION_BED
                            [--verbosity VERBOSITY]
                            [--bal_lrt_alpha BAL_LRT_ALPHA]
                            [--bal_margin BAL_MARGIN]
                            [--filter_std FILTER_STD] [--min_nbins MIN_NBINS]
                            [--ub_nbins UB_NBINS] [--skip_mhbafs]
```

## Main parameters

Here we describe the main parameters. See [reference.md#cluster-bins](../reference.md#cluster-bins) for the full parameter description.

### HMM inference

- **Number of clusters (`--minK`, `--maxK`).** cluster-bins fits the factorial HMM for every K in the closed interval [`--minK` (default 3), `--maxK` (default 30)] and selects the best K by a model score (below). Widen the range if the chosen K lands at either end; narrow it to save runtime when the expected number of distinct states is known.

- **Off-diagonal transition mass (`-t`, default 1e-6).** The off-diagonal transition cost of the HMM transition matrix balances *global* information (RDR/BAF shared across samples) against *local* information (keeping adjacent bins in the same segment). Smaller `-t` favors local continuity (smoother, more contiguous segments); larger `-t` favors global grouping. Reduce it by orders of magnitude for noisier or low-coverage data.

- **Model selection (`--score_method`).** The number of clusters (`K`) is chosen based on model-selection over `--score_method` (`icl`, default, or `bic`). `icl` penalizes overlapping clusters compared `bic` using additonal posterior entropy cost, so it tends to return fewer, better-separated clusters.

> [!TIP]
> Model selection is a heuristic, not a guarantee. We recommend the user to inspect the model-selection score curve plot (`bbc/plots/model_scores.pdf`) alongside the per-K RDR-BAF 1D/2D clustering plots (`bbc/plots/K<K>.pdf`) to see if a better fitted solution may exist and use that solution instead for `compute-cn` - every solution is written to `bbc/labels/bulk<K>.bbc` / `bbc/labels/bulk<K>.seg`.

### Balanced-cluster detection

After HMM state decoding, each cluster is tested for being allelic balanced (BAF = 0.5) or not with an interval likelihood-ratio test on the raw allele counts. `--bal_margin` (default 0.03) is the half-width of the neutral zone `[0.5-δ, 0.5+δ]` treated as a balanced cluster, and `--bal_lrt_alpha` (default 0.05) is the significance level; a cluster is called *balanced* only if it passes in every sample. These clusters anchor the BAF baseline used downstream by `compute-cn`.

> [!TIP]
> Under strong reference-mapping bias the BAF of truly balanced clusters can be pulled off 0.5, so the automatic test may failed to detect them. In this case, we recommend users to manually inspect the BAF centers using 1D/2D clustering plots and decide the balanced cluster by toggling the values of `is_balanced` column (True/False) under `seg` file.

### Post-cluster filtering

Outlier clusters are removed before output: any cluster with fewer than `--min_nbins` bins (default 10) is dropped, and among small clusters (<= `--ub_nbins` bins, default 50) those whose variance exceeds the mean by more than `--filter_std` standard deviations (default 2.0) are treated as variance outliers and removed. Raise `--min_nbins` or lower `--filter_std` to filter more aggressively.

## Output

Clustering results written to `--bbc_dir`: the model-selected $K$ cluster results (`bulk.bbc`, `bulk.seg`), per-$K$ sweeps under `labels/`, `model_scores.tsv`, and diagnostic plots under `plots/`. See [reference.md#output](../reference.md#output) for the full directory tree.
