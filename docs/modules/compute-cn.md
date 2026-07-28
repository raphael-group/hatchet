# `compute-cn`
`compute-cn` performs allele-specific integer copy number and clone proportion deconvolution with regularization using integer linear programming (ILP) or a coordinate descent algorithm.

## Input

The BBC and SEG tables produced by `cluster-bins` (`--bbc`, `--seg`; e.g. `bbc/bulk.bbc` and `bbc/bulk.seg`), plus the reference `--genome_size` and `--region_bed`. See [reference.md#output](../reference.md#output) for the `bbc/` file layout.

## Usage

```console
$ hatchet compute-cn --help
usage: hatchet compute-cn [-h] --result_dir RESULT_DIR --bbc BBC --seg SEG
                          [--mode {both,cd,ilp,cnt_cd}]
                          [--model_select {elbow,bic}] [--force]
                          [--solver {gurobi,cbc}] [--timelimit TIMELIMIT]
                          [--fcn_ci_alpha FCN_CI_ALPHA]
                          [--min_ci_margin MIN_CI_MARGIN]
                          [--obj_type {imf,ci}] [--minClone MINCLONE]
                          [--maxClone MAXCLONE] [--diploid] [--tetraploid]
                          [--reg_term {RAW,MAXCN,DBOX_L1,DBOX_L0,DROOT_SUM,DADJ_SUM}]
                          [--reg_steps REG_STEPS] [--reg_bound REG_BOUND]
                          [--fix_cn_dip FIX_CN_DIP] [--fix_cn_tet FIX_CN_TET]
                          [--zero_cn_thres ZERO_CN_THRES] [--cd_tol CD_TOL]
                          [--no_ampdel] [--num_cnstates NUM_CNSTATES]
                          [-eD DIPLOIDCMAX] [-eT TETRAPLOIDCMAX]
                          [--min_prop MIN_PROP] [--purities PURITIES]
                          [--cd_niters CD_NITERS]
                          [--cd_convergence_iters CD_CONVERGENCE_ITERS]
                          [--cd_nseeds CD_NSEEDS] [--cd_njobs CD_NJOBS]
                          [--cd_seed CD_SEED]
                          [--u_init {dirichlet,bubble,bin_dir}]
                          [--u_dir_alpha U_DIR_ALPHA] [--u_bin_p U_BIN_P]
                          [--solver_threads SOLVER_THREADS]
                          [--tree_file TREE_FILE] [--eps_fit EPS_FIT]
                          [--verbosity VERBOSITY] --genome_size GENOME_SIZE
                          --region_bed REGION_BED [--patient_id PATIENT_ID]
```

## Main parameters

Here we describe the main parameters. See [reference.md#compute-cn](../reference.md#compute-cn) for the full parameter description.

### Fractional copy-number scaling factor estimation

compute-cn first estimates a per-sample RDR scaling factor (gamma) that maps read-depth ratios to fractional copy numbers (FCN), together with tumor purity and clonal CN anchors (`get_scaling_factor`), run separately for the no-WGD (diploid) and WGD (tetraploid) hypotheses. Balanced clusters are taken from the cluster-bins `is_balanced` calls; the user can override the anchors with `--fix_cn_dip` / `--fix_cn_tet`. The decision tree:

- **Baseline (s0).** A user-pinned balanced cluster (`--fix_cn_dip` `1|1` / `--fix_cn_tet` `2|2`), otherwise the largest balanced cluster; sets `gamma = 2 / RD(s0)`.
- **User-pinned anchor.** If a `--fix_cn_dip`/`--fix_cn_tet` cluster is imbalanced, it sets the imbalanced anchor directly.
- **Grid-fit check.** A candidate (purity, gamma) is accepted only if every cluster's observed RDR and BAF are representable by integer copy numbers up to the ploidy's maximum CN.
- **Anchor search.** Otherwise, imbalanced clusters are searched over enumerated CN candidates (filtered by purity validity, an RDR/BAF concordance test, and grid-fit), keeping the one with minimum RDR MSE.
- **Fallback (no-WGD).** If no imbalanced anchor is found, keep only the balanced baseline (s0 as (1,1)) and leave purity to the downstream deconvolution.
- **Fallback (WGD).** Otherwise read a second balanced cluster s1 as (1,1) against s0 as (2,2), giving `gamma = 2 / RD(s1)`, `purity = RD(s0) / RD(s1) - 1`; skipped if no valid s1.

The estimated per-(sample, ploidy) gamma is written to `results/gammas.tsv`.

### Distance-based Constrained Allele-specific Copy-number Factorization (D-CACF)
The (D-CACF) problem is solved by either ILP (`ilp`) or coordinate-descent (`cd` default) algorithm set by `--mode` using an external ILP solver (`--solver`), see [Setup-ILP-Solver](../../README.md#setup-ilp-solver).

#### Main model parameters
- **Maximum copy number (`-eD`/`--diploidcmax`, `-eT`/`--tetraploidcmax`).** Caps the per-segment integer CN at 8 (diploid) / 12 (tetraploid) by default; set to 0 to infer the cap from the scaled fractional copy numbers.

- **Purity and proportions (`--purities`, `--min_prop`).** `--purities` fixes each sample's normal fraction via `sample:purity` pairs; `--min_prop` (default 0.01) is the smallest clone proportion retained.

#### Regularization terms

A regularization term avoids overfitting the CN solution to observation noise. `--reg_term` selects the penalty; compute-cn sweeps `--reg_steps` (default 15) penalty weights up to `--reg_bound` (default 0.15) and picks the elbow of the fit-vs-penalty (Pareto) curve. `--num_cnstates` can additionally cap the number of distinct CN states per clone (`-1` = unconstrained).

For cluster $m$ (weight $w_m$), tumor clones $n = 1,\dots,N$ (clone $0$ = normal), with A/B-allele copy numbers $a_{m,n}$ / $b_{m,n}$:

| `--reg_term` | Penalty | Description |
|---|---|---|
| `RAW` | $0$ | No regularization; fit term only. |
| `DBOX_L1` (default) | $\sum_m w_m\left[(\max_n a_{m,n} - \min_n a_{m,n}) + (\max_n b_{m,n} - \min_n b_{m,n})\right]$ | Per-cluster allelic CN range (max - min) across tumor clones. |
| `DBOX_L0` | $\sum_m w_m \cdot \mathbf{1}[\text{span}_m > 0]$ | Per-cluster indicator of subclonality (nonzero allelic spread); $\text{span}_m$ is the `DBOX_L1` term for cluster $m$. |
| `MAXCN` | $\sum_m w_m (\max_n a_{m,n} + \max_n b_{m,n})$ | Per-cluster maximum tumor-clone total copy number. |
| `DROOT_SUM` | $\sum_m w_m \sum_n \left(\lvert a_{m,n} - a_{m,0}\rvert + \lvert b_{m,n} - b_{m,0}\rvert\right)$ | L1 distance of each tumor clone from the normal clone. |
| `DADJ_SUM` | $\sum_m w_m \sum_{n_1<n_2} \left(\lvert a_{m,n_1} - a_{m,n_2}\rvert + \lvert b_{m,n_1} - b_{m,n_2}\rvert\right)$ | Total pairwise L1 distance between tumor-clone pairs. |

### Ploidy and clone number

- **Ploidy (`--diploid`, `--tetraploid`).** **diploid** assumes no WGD, and **tetraploid** assumes one WGD. If neither flag is set, compute-cn solves both instances.

- **Number of clones (`--minClone`, `--maxClone`).** Integer CN and clone proportions are inferred for every clone count n in [`--minClone` (default 2), `--maxClone` (default 4)]. We recommend to set higher `--maxClone` if more than 1 tumor samples are provided.

- **Model selection (`--model_select`).** For each ploidy case, the clone number is picked by either BIC (default) or an elbow along the log-likelihood curve using the `kneed` library. The log-likelihood scores each candidate by the observed RDR (Gaussian) and phased allele counts (Beta-Binomial) given the deconvolved integer states and clone proportion. with per-cluster RDR variance and BAF dispersion provided from the seg file. Then, the ploidy with lower number of clones is chosen as final model-selected solution.

> [!TIP]
> Model selection is a heuristic. We recommend users to inspect the selection curves (`results/plots/model_selection.pdf`) alongside `results/summary.tsv`, and if another (ploidy, n) fits better, use its solution instead - every candidate is written to `results/chosen.<ploidy>.*.ucn` and `results/results.<ploidy>.n*.ucn.tsv`.

## Output

Copy-number solutions written to `--result_dir`: the model-selected `best.bbc.ucn` / `best.seg.ucn`, per-(ploidy, n) solutions, `gammas.tsv`, `summary.tsv`, and plots under `plots/`. See [reference.md#output](../reference.md#output) for the full directory tree.
