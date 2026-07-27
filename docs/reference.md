# Reference

## Table of Contents
1. [Input](#input)
2. [Environment variables](#environment-variables)
3. [Parameters](#parameters)
    - [cluster-bins](#cluster-bins)
    - [compute-cn](#compute-cn)
    - [plot-cn](#plot-cn)
    - [plot-panel](#plot-panel)
4. [Output](#output)

## Input
```
<bb_dir>/
  bb.tsv.gz         # BED format bin meta-informations, one row per bin, aligned with matrix rows.
                    #   columns: #CHR, START, END, region_id, switchprobs, #SNPS
  sample_ids.tsv    # dataset meta-informations, one row per sample, aligned with matrix columns.
                    #   columns: SAMPLE, sample_type (normal|tumor), [assay_type]
  bb.rdr.npz        # (bins, tumor samples) float  - NumPy format, read-depth ratio, tumor samples only
  bb.depth.npz      # (bins, samples) float        - NumPy format, read depth, all samples
  bb.Aallele.npz    # (bins, samples) int          - NumPy format, A-haplotype allele counts
  bb.Ballele.npz    # (bins, samples) int          - NumPy format, B-haplotype allele counts
  bb.Tallele.npz    # (bins, samples) int          - NumPy format, total allele counts
```

`sample_ids.tsv` columns (one row per sample; row order defines the matrix column order):

| Column | Required | Description |
|---|---|---|
| `SAMPLE` | yes | Sample name; also the matrix column label |
| `sample_type` | yes | `normal` or `tumor`; the tumor-only `bb.rdr.npz` holds just the `tumor` columns |
| `assay_type` | no | Optional sequencing-assay label; samples are grouped by assay |

## Environment variables

| Variable | Default | Description |
|---|---|---|
| `HATCHET_DISABLE_CPP` | `0` | When set to `1`, `true`, or `yes` (case-insensitive), `cluster-bins` uses the pure-Python (Numba) HMM backend instead of the compiled C++ extension (`_hmm_cpp`). Useful for debugging or when the extension fails to build. Example: `HATCHET_DISABLE_CPP=1 hatchet cluster-bins ...` |

## Parameters

Defaults below are the CLI defaults from `hatchet_parser.py`. Use
`hatchet <command> --help` for the authoritative list.

### cluster-bins

| Parameter | Default | Description |
|---|---|---|
| `--bb_dir` | *(required)* | Directory containing NPZ count matrices and `bb.tsv.gz` input files |
| `--bbc_dir` | *(required)* | Output directory for BBC/seg TSV files |
| `--genome_size` | *(required)* | Reference chromosome sizes file |
| `--force` | False | Re-run even if results already exist (default: skip existing) |
| `--minK` | 3 | Minimum number of HMM cluster states |
| `--maxK` | 30 | Maximum number of HMM cluster states |
| `-t` | 1e-6 | Initial off-diagonal transition mass |
| `--restarts` | 10 | Number of random restarts per K |
| `--top_restarts` | *(= restarts)* | Number of top-scoring inits to run full EM on |
| `--n_local_trials` | 3 | Candidate bins evaluated per k-means++ seeding step |
| `--niters` | 50 | Number of EM iterations per restart |
| `--decode_method` | `map` | HMM decoding: `viterbi` (most-likely path) or `map` (marginal per-bin posterior) |
| `--score_method` | `icl` | Model selection score: `bic` or `icl` |
| `--score_criteria` | `min` | How to pick K from the score curve: `min`, `elbow`, or `margin-<int>` |
| `--init_method` | `cna_plus_plus` | Init method: `cna_plus_plus` (HMM-aware) or `kmeans_plus_plus` (sklearn KMeans++) |
| `--training_method` | `baum_welch` | HMM training: `baum_welch` (soft EM) or `viterbi` (hard EM) |
| `--free_baf_c0` | False | Allow cluster-0 BAF to update during EM (default: fixed at 0.5) |
| `--min_tau` | 1 | Minimum Beta-Binomial dispersion tau |
| `--max_tau` | 1e6 | Maximum Beta-Binomial dispersion tau |
| `--share_tau` | True | Share BB dispersion tau across clusters within a sample (`--no-share_tau` for per-cluster) |
| `--tau_iters` | 3 | Number of EM iterations during which BAF dispersion tau is updated |
| `--baf_eps` | 1e-3 | BAF mean Brent search bounds `[baf_eps, 1-baf_eps]`; sequencing error floor |
| `--min_covar` | 1e-3 | Minimum RDR variance floor applied after each M-step |
| `--ig_alpha` | 10.0 | Inverse-gamma prior shape parameter for RDR variance updates |
| `--log_rdr` | False | Use log(RDR) instead of raw RDR in the Gaussian emission |
| `--seed` | 42 | Random seed for HMM init step |
| `--bal_lrt_alpha` | 0.05 | Significance level for the balanced-cluster interval LRT |
| `--bal_margin` | 0.03 | Half-width of neutral zone `[0.5-δ, 0.5+δ]` for the balanced-cluster test |
| `--filter_std` | 2.0 | Filter clusters whose variance deviates from mean by `filter_std × std` |
| `--min_nbins` | 10 | Remove clusters with fewer than `min_nbins` bins |
| `--ub_nbins` | 50 | Variance-outlier filtering only applies to clusters with #bins ≤ `ub_nbins` |
| `--skip_mhbafs` | False | Skip minor-haplotype BAF folding after decoding |
| `--verbosity` | 0 | Verbose level: 0, 1, or 2 |

### compute-cn

| Parameter | Default | Description |
|---|---|---|
| `--bbc` | *(required)* | Filename for BBC table (e.g., `bbc/bulk.bbc`) |
| `--seg` | *(required)* | Filename for SEG table (e.g., `bbc/bulk.seg`) |
| `--result_dir` | *(required)* | Output directory for computed CN results |
| `--genome_size` | *(required)* | Reference chromosome sizes file |
| `--region_bed` | *(required)* | Reference chromosome BED file |
| `--mode` | `cd` | Solver mode: `cd`, `ilp`, `both`, or `cnt_cd` (**Experimental**) |
| `--solver` | `gurobi` | ILP solver backend: `gurobi` or `cbc` |
| `--model_select` | `bic` | Clone-number/ploidy selection: `elbow` or `bic` |
| `--force` | False | Re-solve even if results already exist (default: skip existing) |
| `--timelimit` | None | ILP solver time limit in seconds |
| `--obj_type` | `imf` | Fitting objective: `imf` (weighted L1) or `ci` (CI-violation hinge) |
| `--fcn_ci_alpha` | 0.05 | Significance level for the FCN confidence interval (0.05 -> 95% CI) |
| `--min_ci_margin` | 0.1 | Hard minimum CI half-width in FCN space |
| `--minClone` | 2 | Minimum number of tumor clones |
| `--maxClone` | 4 | Maximum number of tumor clones |
| `--diploid` | False | Solve under diploid assumption |
| `--tetraploid` | False | Solve under tetraploid/WGD assumption |
| `--reg_term` | `DBOX_L1` | Regularizer: `RAW`, `MAXCN`, `DBOX_L1`, `DBOX_L0`, `DROOT_SUM`, or `DADJ_SUM` |
| `--reg_steps` | 15 | Number of steps in the regularization path |
| `--reg_bound` | 0.15 | Maximum penalty weight for the regularization path |
| `--fix_cn_dip` | None | Fix diploid cluster CN states, e.g. `6:2|0;8:3|1` |
| `--fix_cn_tet` | None | Fix tetraploid cluster CN states, e.g. `6:4|2` |
| `--zero_cn_thres` | 0.005 | Clusters with weight ≥ this fraction of total cannot take a (0,0) CN state |
| `--no_ampdel` | False | Disable the amp/del symmetry constraint |
| `--num_cnstates` | -1 | Constrain the number of distinct CN states per clone (-1 = unconstrained) |
| `-eD` / `--diploidcmax` | 8 | Max copy number for diploid mode (0 = inferred from scaled FCN) |
| `-eT` / `--tetraploidcmax` | 12 | Max copy number for tetraploid mode (0 = inferred from scaled FCN) |
| `--min_prop` | 0.01 | Minimum clone proportion |
| `--purities` | None | Semicolon-separated `sample:purity` pairs; fixes normal-clone proportion to 1 − purity |
| `--cd_niters` | 10 | CD: max outer iterations per seed |
| `--cd_convergence_iters` | 2 | CD: consecutive convergence iterations required to stop |
| `--cd_tol` | 0.001 | CD: stop when U-step objective changes less than this |
| `--cd_nseeds` | 400 | CD: number of random restarts |
| `--cd_njobs` | 8 | CD: number of parallel worker processes |
| `--cd_seed` | 42 | CD: random seed for reproducibility |
| `--u_init` | `dirichlet` | U initialization: `dirichlet`, `bubble`, or `bin_dir` |
| `--u_dir_alpha` | *(solver default)* | Dirichlet alpha for U initialization; lower = sparser |
| `--u_bin_p` | *(solver default)* | `bin_dir`: per-cell Bernoulli presence probability |
| `--solver_threads` | *(solver default)* | Max threads per solver call (Gurobi); set to 1 for parallel CD workers |
| `--tree_file` | None | **Experimental** (`cnt_cd`): Newick tree file; if omitted, enumerate all unlabeled shapes |
| `--eps_fit` | 0.01 | **Experimental** (`cnt_cd`): fit tolerance for the C-step CNT lexicographic bound |
| `--plot_ascn` | False | Plot CN profile with the allele-CN row scheme |
| `--patient_id` | `panel` | Output filename prefix for per-(ploidy, n) plots |
| `--verbosity` | 0 | Verbose level: 0, 1, or 2 |

### plot-cn

Standalone plotting command. `compute-cn` already emits plots under `results/plots/`;
use `plot-cn` to regenerate or customize a specific solution.

| Parameter | Default | Description |
|---|---|---|
| `--bbc` | *(required)* | BBC UCN table (e.g., `results/best.bbc.ucn`) |
| `--seg` | *(required)* | SEG UCN table (e.g., `results/best.seg.ucn`) |
| `--genome_size` | *(required)* | Reference chromosome sizes file |
| `--region_bed` | *(required)* | Reference chromosome BED file |
| `-g` / `--gamma_file` | *(required)* | Gamma scaling-factor file from `compute-cn` |
| `--ploidy` | *(required)* | `diploid` or `tetraploid` (selects the gamma column) |
| `-O` / `--plot_dir` | *(required)* | Output directory for figures |
| `-s` / `--solfile` | None | Optional solution file to override CN states in the BBC table |
| `--img_type` | `png` | File format: `pdf`, `png`, or `svg` |
| `--dpi` | 500 | Image resolution |
| `--transparent` | False | Transparent background |
| `--keep_gap` | False | Keep gap regions in the plot |
| `--tail_alpha` | 0.8 | Transparency of the tail region per CN state |
| `--center_alpha` | 1.0 | Transparency of the center region per CN state |
| `--onetail_area` | 0.025 | Per-tail area per CN state used to set transparency |
| `--maxlim_fcn` | 30 | Figure axis limit for FCN |
| `--plot_ascn` | False | Plot CN profile with the allele-CN row scheme |
| `--patient_id` | *(none)* | Output filename prefix for combined plots |

### plot-panel

Composes a multi-sample CN panel from several per-sample BBC UCN solutions listed in
a panel TSV (`--panel_file`).

| Parameter | Default | Description |
|---|---|---|
| `--panel_file` | *(required)* | Panel TSV listing per-sample BBC UCN paths |
| `--genome_size` | *(required)* | Reference chromosome sizes file |
| `--region_bed` | *(required)* | Reference chromosome BED file |
| `-o` / `--out_file` | *(required)* | Output figure path (e.g., `panel.svg`) |
| `--width` | 20 | Panel image width |
| `--height` | 1 | Panel image height per row |
| `--show_clone_name` | False | Draw clone names on the CN profile |
| `--show_prop` | False | Draw clone proportions on the CN profile |
| `--show_ploidy` | False | Draw per-clone ploidy on the CN profile |
| `--min_prop` | 0.01 | Hide tumor clones below this proportion from the panel |
| `--dpi` | 300 | Image resolution |
| `--transparent` | False | Transparent background |
| `--title` | `panel` | Plot title |
| `--plot_1d2d` | False | Also run `plot-cn` per panel row (requires a `PATH_TO_BBC` column) |
| `--plot_summary` | False | Emit per-sample purity + ploidy barplots (one page per metric per `cancer_type`) |

## Output

```
output/my_sample/
  bbc/                             # cluster-bins output
    bulk.bbc                       # per-bin cluster assignments (optimal K)
    bulk.seg                       # per-cluster summary statistics (optimal K)
    labels/
      bulk<K>.bbc                  # per-bin assignments for each swept K
      bulk<K>.seg                  # per-cluster summary for each swept K
    cluster_infos/                 # per-K cluster-label TSVs
    plots/                         # ELBO traces, RDR-BAF scatter, model score
    model_scores.tsv               # BIC and ICL scores across K
  results/                         # compute-cn output
    best.bbc.ucn                   # model-selected CN solution (per-bin)
    best.seg.ucn                   # model-selected CN solution (per-cluster)
    chosen.<ploidy>.bbc.ucn        # best-n solution per ploidy (per-bin)
    chosen.<ploidy>.seg.ucn        # best-n solution per ploidy (per-cluster)
    results.<ploidy>.n*.bbc.ucn.tsv  # every (ploidy, n) solution (per-bin)
    results.<ploidy>.n*.seg.ucn.tsv  # every (ploidy, n) solution (per-cluster)
    gammas.tsv                     # RDR scaling factors per sample and ploidy
    summary.tsv                    # fit/regularization metrics per solution
    plots/
      scaling_2d.pdf               # RDR-vs-BAF scaling diagnostic
      model_selection.pdf          # Pareto front + elbow/BIC selection page
      <ploidy>_n*/                 # per-solution CN plots
  logs/                            # Snakemake log_dir
    cluster_bins.log
    compute_cn.log
    cluster_bins.benchmark.tsv     # wall time, CPU, peak memory
    compute_cn.benchmark.tsv
```

The final model-selected outputs are `results/best.bbc.ucn` and `results/best.seg.ucn`
(these are the Snakemake workflow targets). `compute-cn` produces its own plots under
`results/plots/`; the standalone `hatchet plot-cn` command can regenerate or customize them.
