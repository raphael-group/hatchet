# Changelog

HATCHet3 [commit xxx] is the next generation of HATCHet2 that infers allele-specific and clone-specific CNAs from multiple tumor samples of short-read WGS, WES, and long-read sequencing data.

* HATCHet2 commit: [`84ebfba`](https://github.com/raphael-group/hatchet/tree/84ebfbac765a8329899a3d03c3791324bbb1fe3e).

All notable changes to HATCHet3 are documented in this file.

### Pipeline

**Changed**

- **Modular Snakemake workflow.** One Snakemake rule per stage (`run_cluster_bins` →
  `run_compute_cn`) with a parser-introspecting config→CLI forwarder that renders only
  valid flags from nested `cluster-bins:`/`compute-cn:` config blocks (`render_cli_args`,
  [Snakefile:12](Snakefile#L12)).
  <br>_HATCHet2:_ monolithic `run.py` driver sequentially invoking each step, configured
  by a single INI file via `configparser`
  ([run.py:30-31, 211-268](https://github.com/raphael-group/hatchet/blob/84ebfbac765a8329899a3d03c3791324bbb1fe3e/src/hatchet/utils/run.py#L211-L268)).

### Preprocessing

**Changed**

- **Standalone Universal Genotyping pipeline.** Read counting, panel-based SNP
  genotyping, population/long-read-based phasing, bias correction, and segmentation are now a
  dedicated Snakemake pipeline [Universal-Genotyping-Pipeline](https://github.com/raphael-group/Universal-Genotyping-Pipeline/tree/dev) supporting multiple sequencing platforms,
  where each external tool is an explicit, reproducible, independently-resumable rule,
  emitting the `bb_dir` matrices HATCHet3 consumes
  (`universal-genotyping/workflow/Snakefile`, `universal-genotyping/workflow/rules/`).
  <br>_HATCHet2:_ the error-prone pattern of driving external binaries in-process (e.g.
  calling `bcftools`/`samtools`/`tabix` from within Python), with subprocess calls buried
  inside Python wrappers
  ([count_alleles.py](https://github.com/raphael-group/hatchet/blob/84ebfbac765a8329899a3d03c3791324bbb1fe3e/src/hatchet/utils/count_alleles.py),
  [count_reads.py](https://github.com/raphael-group/hatchet/blob/84ebfbac765a8329899a3d03c3791324bbb1fe3e/src/hatchet/utils/count_reads.py),
  [phase_snps.py](https://github.com/raphael-group/hatchet/blob/84ebfbac765a8329899a3d03c3791324bbb1fe3e/src/hatchet/utils/phase_snps.py)).

### cluster-bins

**Changed**

- **Emission model → phase-aware Beta-Binomial factorial HMM.** Joint hidden state over K
  clusters × 2 haplotype-phase orientations, with a Gaussian RDR emission and a
  **count-based Beta-Binomial** minor-haplotype BAF emission (dispersion `τ`) evaluated
  under both phase orientations, so allelic imbalance is inferred rather than folded a
  priori (`_fwd_bwd_seg`,
  [hmm_fwd_bwd.py:47](src/hatchet/cluster_bins/hmm/hmm_fwd_bwd.py#L47);
  `compute_loglik`,
  [hmm_likelihoods.py:12](src/hatchet/cluster_bins/hmm/hmm_likelihoods.py#L12)).
  <br>_HATCHet2:_ `hmmlearn.GaussianHMM` with a **multivariate Gaussian** over standardized
  [RDR, pre-folded BAF] — no count model, no phase state
  ([cluster_bins.py:176-193](https://github.com/raphael-group/hatchet/blob/84ebfbac765a8329899a3d03c3791324bbb1fe3e/src/hatchet/utils/cluster_bins.py#L176-L193)).
- **Genetic-map phase coupling.** Per-bin switch/stay probabilities from the upstream
  genetic map gate the phase axis of the transition kernel, an explicit factorial phase
  state
  ([hmm_fwd_bwd.py:85](src/hatchet/cluster_bins/hmm/hmm_fwd_bwd.py#L85)).
  <br>_HATCHet2:_ HMM transitions run over cluster states only; BAF is pre-folded to the
  minor-allele fraction, so no haplotype orientation is modeled
  ([cluster_bins.py:176-191](https://github.com/raphael-group/hatchet/blob/84ebfbac765a8329899a3d03c3791324bbb1fe3e/src/hatchet/utils/cluster_bins.py#L176-L191)).
- **Initialization → likelihood-based D2 (cna++).** Centroids seeded by squared
  negative-log-likelihood distance (k-means++/D2 in likelihood space) with a fixed diploid
  anchor at (RDR=1, BAF=0.5) (`init_hmm_cna_plus_plus`,
  [hmm_init.py:15](src/hatchet/cluster_bins/hmm/hmm_init.py#L15)).
  <br>_HATCHet2:_ `hmmlearn` default init (`init_params="mc"`) seeds means with sklearn
  **KMeans (k-means++)** in Euclidean [RDR, BAF] feature space, ignoring the emission
  likelihood
  ([cluster_bins.py:186-193](https://github.com/raphael-group/hatchet/blob/84ebfbac765a8329899a3d03c3791324bbb1fe3e/src/hatchet/utils/cluster_bins.py#L186-L193)).
- **Compute backend → C++/OpenMP Baum-Welch.** The full EM loop (emissions,
  forward-backward, closed-form RDR M-step, Brent τ/BAF M-steps) is compiled via pybind11
  with OpenMP parallelism, with a NumPy fallback (`run_hmm`,
  [_hmm_cpp/bindings.cpp:101](src/hatchet/cluster_bins/hmm/_hmm_cpp/bindings.cpp#L101)).
  <br>_HATCHet2:_ single-threaded pure-Python `hmmlearn` `model.fit`, no compiled kernel
  ([cluster_bins.py:193](https://github.com/raphael-group/hatchet/blob/84ebfbac765a8329899a3d03c3791324bbb1fe3e/src/hatchet/utils/cluster_bins.py#L193)).
- **Model selection → ICL/BIC/elbow over K.** K is swept and selected by ICL (default),
  BIC, or a kneedle-elbow criterion (`score_model`, `model_select_K`,
  [hmm_utils.py:54](src/hatchet/cluster_bins/hmm/hmm_utils.py#L54)).
  <br>_HATCHet2:_ silhouette score or BIC only — no ICL — selected via `state_selection`
  ([cluster_bins.py:196-199](https://github.com/raphael-group/hatchet/blob/84ebfbac765a8329899a3d03c3791324bbb1fe3e/src/hatchet/utils/cluster_bins.py#L196-L199)).
- **Balanced-cluster detection → statistical test.** Each cluster is tested for allelic
  balance (BAF≈0.5) by an interval likelihood-ratio test with parametric bootstrap on a
  symmetric Beta-Binomial mixture, labeling `is_balanced` across all samples
  (`label_balanced_clusters`,
  [cluster_utils.py:300](src/hatchet/cluster_bins/cluster_utils.py#L300));
  balanced clusters then seed the diploid RDR baseline and are pinned as solver seeds
  (`get_scaling_factor`,
  [scaling.py:252](src/hatchet/compute_cn/scaling.py#L252)).
  <br>_HATCHet2:_ a user-defined hard cutoff `diploidbaf` collapses BAF to exactly 0.5 when
  `|0.5 − BAF| < threshold` for all samples
  ([cluster_bins.py:285](https://github.com/raphael-group/hatchet/blob/84ebfbac765a8329899a3d03c3791324bbb1fe3e/src/hatchet/utils/cluster_bins.py#L285)).

### compute-cn

**Added**

- **CNT tree-constrained solver.** New `cnt_cd` mode alternates a per-chromosome
  copy-number-tree MILP (lexicographic: minimize fit loss, then tree parsimony under a
  near-optimal-fit constraint) with a global U-step LP over enumerated clone-tree
  topologies (`solve_c_step`,
  [cnt_model.py:229](src/hatchet/compute_cn/solve/cnt_model.py#L229);
  `run_coordinate_descent`,
  [inference.py:369](src/hatchet/compute_cn/solve/inference.py#L369)).
- **Selectable regularizers.** The objective interpolates the fit term with a chosen
  penalty — `MAXCN`, `DBOX_L1`/`DBOX_L0` (subclonal allelic spread), `DROOT_SUM`,
  `DADJ_SUM`, or unregularized `RAW` — swept over a regularization path
  (`build_regularization`,
  [regularization.py:169](src/hatchet/compute_cn/solve/regularization.py#L169);
  `build_final_objective`,
  [objectives.py:32](src/hatchet/compute_cn/solve/objectives.py#L32)).

**Changed**

- **Objective & model selection.** The regularized objective above is paired with automated
  clone-number/ploidy selection by a kneedle elbow on the log-likelihood (or BIC) and along
  the (fit vs. regularization) Pareto front (`model_selection_ploidy`,
  `model_select_elbow_from_regularization`,
  [model_select.py:128](src/hatchet/compute_cn/model_select.py#L128)).
  <br>_HATCHet2:_ a single fixed objective inside `ILPSubset`, with the number of clones
  `n` user-specified (no regularization path, no Pareto/elbow selection in the solver)
  ([solve/cd.py:27](https://github.com/raphael-group/hatchet/blob/84ebfbac765a8329899a3d03c3791324bbb1fe3e/src/hatchet/utils/solve/cd.py#L27),
  [solve/ilp_subset.py](https://github.com/raphael-group/hatchet/blob/84ebfbac765a8329899a3d03c3791324bbb1fe3e/src/hatchet/utils/solve/ilp_subset.py)).
- **Solver backend.** MILP/LP models are abstracted through Pyomo, running on Gurobi **or**
  open-source CBC (`create_solver`,
  [inference.py:48](src/hatchet/compute_cn/solve/inference.py#L48)).
  <br>_HATCHet2:_ Gurobi-only, via a compiled C++ ILP and a Python coordinate-descent wrapper
  ([solve/cd.py:67](https://github.com/raphael-group/hatchet/blob/84ebfbac765a8329899a3d03c3791324bbb1fe3e/src/hatchet/utils/solve/cd.py#L67),
  [src/solve.cpp](https://github.com/raphael-group/hatchet/blob/84ebfbac765a8329899a3d03c3791324bbb1fe3e/src/solve.cpp)).

### Plotting

**Added**

- **1D/2D library.** Genome-wide 1D RDR/BAF tracks and a 2D RDR-vs-BAF joint scatter with
  cluster coloring, gap masking, and density-based transparency (`plot_1d`, `plot_2d`,
  [plot_1d2d.py:147](src/hatchet/plot/plot_1d2d.py#L147)).

**Changed**

- **Output format.** Figures emit editable-font vector SVG (`svg.fonttype="none"`) with
  dpi/transparency control (`run`,
  [plot_cn.py:34](src/hatchet/plot/plot_cn.py#L34)).
  <br>_HATCHet2:_ raster **PNG-only** output, no vector/editable-font support
  ([plot_cn_1d2d.py:389](https://github.com/raphael-group/hatchet/blob/84ebfbac765a8329899a3d03c3791324bbb1fe3e/src/hatchet/utils/plot_cn_1d2d.py#L389)).

---

## TODO / Known limitations

### Universal-genotyping part
1. Het SNPs genotyping from high purity tumor sample without matched-normal.

### cluster-bins

- **Post-merge step after HMM decoding.** ICL/BIC under-penalize complexity at the very
  large observation counts here (the data log-likelihood scales with N while the penalty
  scales with log N), so K-selection tends to over-segment; a post-decoding merge of
  statistically indistinguishable clusters — extending the Beta-Binomial LRT in
  `label_balanced_clusters` to pairwise cluster equivalence — would yield more robust
  cluster counts than information-criterion selection alone.

### compute-cn

- **`cnt_cd` runtime.** The tree-constrained solver is the throughput bottleneck: it
  enumerates clone-tree topologies × Dirichlet seeds and solves a per-chromosome MILP
  C-step each iteration, scaling poorly with clone number and genome size; it needs
  tree-space pruning, MILP warm-starting/relaxation, or a cheaper C-step to be practical
  at scale.
