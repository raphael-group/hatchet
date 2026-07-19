# Changelog

## [3.0.0] - 2026-07-19

**Changed**

- **Gaussian RDR + Beta-Binomial BAF factorial HMM with phase switch correction.** This replaces HATCHet2's [`combine_counts.py`](https://github.com/raphael-group/hatchet/blob/84ebfbac765a8329899a3d03c3791324bbb1fe3e/src/hatchet/utils/combine_counts.py) and [`cluster_bins.py`](https://github.com/raphael-group/hatchet/blob/84ebfbac765a8329899a3d03c3791324bbb1fe3e/src/hatchet/utils/cluster_bins.py). The new segmentation algorithm models phasing switch errors as prior information through genetic-map-derived switch probabilities, and is implemented with a C++/OpenMP backend (`run_baum_welch`, [hmm_model.py:143](src/hatchet/cluster_bins/hmm/hmm_model.py#L143)).
- **HMM model selection.** Added ICL scoring to account for cluster-assignment uncertainty when selecting the cluster number K (`model_select_K`, [hmm_utils.py:6](src/hatchet/cluster_bins/hmm/hmm_utils.py#L6)).
- **Balanced-cluster identification.** Each HMM-decoded cluster is tested for allelic balance (BAF≈0.5) by an interval likelihood-ratio test with parametric bootstrap on a symmetric Beta-Binomial mixture (`label_balanced_clusters`, [cluster_utils.py:300](src/hatchet/cluster_bins/cluster_utils.py#L300)). The `is_balanced` column in the `SEG.UCN` file marks the chosen balanced clusters. This replaces HATCHet2's user-defined hard cutoff `diploidbaf` ([cluster_bins.py:446](https://github.com/raphael-group/hatchet/blob/84ebfbac765a8329899a3d03c3791324bbb1fe3e/src/hatchet/utils/cluster_bins.py#L446)).
- **Scaling factor estimation.** A **decision-tree strategy** identifies a clonal imbalanced cluster to infer the RDR scaling factor (gamma), tumor purity, and clonal CN anchors (`get_scaling_factor`, [scaling.py:8](src/hatchet/compute_cn/scaling.py#L8)). This replaces HATCHet2's `scale_rdr`, which derived a single purity/scale in closed form from user-supplied `--clonal` cluster copy numbers ([solve/utils.py:63](https://github.com/raphael-group/hatchet/blob/84ebfbac765a8329899a3d03c3791324bbb1fe3e/src/hatchet/utils/solve/utils.py#L63)).
  - **Baseline.** A user-pinned balanced cluster (`--fix_cn_dip` `1|1` / `--fix_cn_tet` `2|2`), or otherwise the largest balanced cluster.
  - **User-pinned anchor.** If provided, a `--fix_cn_dip`/`--fix_cn_tet` cluster sets the imbalanced anchor directly.
  - **Grid-fit check.** Shared filter accepting a candidate (purity, gamma) only if every cluster's observed RDR and BAF fall within the ranges representable by integer copy numbers up to the ploidy's maximum CN under that scaling.
  - **Anchor search.** Otherwise, imbalanced clusters are searched over enumerated CN candidates, filtered by purity validity, a statistical RDR/BAF concordance test, and grid-fit, then pick the candidate imbalanced cluster with minimum MSE between observed and CN-expected RDR.
  - **Fallback (no-WGD).** If no imbalanced anchor is found, keep only the balanced baseline (s0 read as (1,1)) and leave purity to be inferred by the downstream deconvolution.
  - **Fallback (WGD).** If no imbalanced anchor is found, read a second balanced cluster s1 as (1,1) against the baseline s0 as (2,2), giving gamma = 2/RD(s1) and purity = RD(s0)/RD(s1) - 1. If no such s1 gives a valid bounded purity and grid-fit, WGD path is skipped.
- **Clone-number model selection.** For each ploidy case, the clone number is picked by an elbow on the log-likelihood (or BIC) vs. clone-number curve using the `kneed` library (`model_selection_ploidy`, [model_select.py:128](src/hatchet/compute_cn/model_select.py#L128)). This replaces HATCHet2's model selection against ILP objectives.
- **[Universal-Genotyping-Pipeline](https://github.com/raphael-group/Universal-Genotyping-Pipeline).** This replaces the error-prone preprocessing sub-modules in HATCHet2: [count_alleles.py](https://github.com/raphael-group/hatchet/blob/84ebfbac765a8329899a3d03c3791324bbb1fe3e/src/hatchet/utils/count_alleles.py), [count_reads.py](https://github.com/raphael-group/hatchet/blob/84ebfbac765a8329899a3d03c3791324bbb1fe3e/src/hatchet/utils/count_reads.py), and [phase_snps.py](https://github.com/raphael-group/hatchet/blob/84ebfbac765a8329899a3d03c3791324bbb1fe3e/src/hatchet/utils/phase_snps.py).
- **[Modular Snakemake workflow](./Snakefile).** One rule per stage replaces HATCHet2's monolithic [run.py:30-31, 211-268](https://github.com/raphael-group/hatchet/blob/84ebfbac765a8329899a3d03c3791324bbb1fe3e/src/hatchet/utils/run.py#L211-L268) driver.
- **1D/2D plotting library.** Genome-wide 1D RDR/BAF tracks and a 2D RDR-vs-BAF joint scatter
  with cluster coloring, marginal distribution, gap masking, and density-based transparency (`plot_1d`, `plot_2d`,
  [plot_1d2d.py:147](src/hatchet/plot/plot_1d2d.py#L147)).

**Added**

- **Integer CN deconvolution with regularization.** CN deconvolution solver is extended with multiple selectable regularizers (`MAXCN`, `DBOX_L1`/`DBOX_L0`, `DROOT_SUM`, `DADJ_SUM`, or unregularized `RAW`, default: `DBOX_L1`) to avoid overfitting the observation noise. The pareto-front solution is picked using the `kneed` library.
  - `MAXCN`: penalizes the per-cluster maximum total copy number.
  - `DBOX_L{0,1}`: penalizes the `L{0,1}` distance for maximum and minimum allele CN for tumor clones per cluster.
  - `DROOT_SUM`: penalizes the L1 distance of each tumor clone's copy numbers from the normal clone.
  - `DADJ_SUM`: penalizes the total pairwise L1 distance between all tumor-clone pairs.