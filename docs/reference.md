# Reference

## Output Structure

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
    results.<ploidy>.n*.bbc.ucn.tsv   # per-n solutions
  summary/                # plot-cn output
    <ploidy>_n*/         # per-solution plots
      <sample>.1D.png   # 1D genome-wide CN profile
      <sample>.2D.png   # 2D RDR-vs-BAF scatter plot
```

## Parameter Reference

## cluster-bins

| Parameter | Default | Description |
|---|---|---|
| `--bb_dir` | *(required)* | Directory containing NPZ count matrices and `bb.tsv.gz` input files |
| `--bbc_dir` | *(required)* | Output directory for BBC/seg TSV files |
| `--genome_size` | *(optional)* | Reference chromosome sizes file |
| `--minK` | 3 | Minimum number of HMM cluster states |
| `--maxK` | 30 | Maximum number of HMM cluster states |
| `-t` | 1e-6 | Initial off-diagonal transition mass |
| `--restarts` | 30 | Number of random restarts per K |
| `--top_restarts` | 30 | Number of top-scoring inits to run full EM on |
| `--n_local_trials` | 3 | Candidate bins to evaluate per k-means++ seeding step (1 = original single-draw) |
| `--niters` | 50 | Number of EM iterations per restart |
| `--decode_method` | `map` | HMM decoding method: `viterbi` (most-likely path) or `map` (marginal per-bin posterior) |
| `--score_method` | `icl` | Model selection criterion: `bic` or `icl` |
| `--init_method` | `cna_plus_plus` | Initialization method: `cna_plus_plus` (HMM-aware seeding) or `kmeans_plus_plus` (sklearn KMeans++) |
| `--min_tau` | 50 | Minimum Beta-Binomial dispersion tau |
| `--max_tau` | 200 | Maximum Beta-Binomial dispersion tau |
| `--baf_eps` | 1e-3 | BAF mean Brent search bounds [baf_eps, 1-baf_eps]; related to sequencing error floor |
| `--min_covar` | 1e-3 | Minimum RDR variance floor applied after each M-step |
| `--ig_alpha` | 10.0 | Inverse-gamma prior shape parameter for RDR variance updates |
| `--tau_iters` | 1 | Number of EM iterations during which BAF dispersion tau is updated |
| `--seed` | 42 | Random seed for HMM init step |
| `--log_rdr` | False | Use log(RDR) instead of raw RDR in the Gaussian emission |
| `--verbosity` | 0 | Verbose level: 0, 1, or 2 |

## compute-cn

| Parameter | Default | Description |
|---|---|---|
| `--bbc` | *(required)* | Filename for BBC table (e.g., `bbc/bulk.bbc`) |
| `--seg` | *(required)* | Filename for SEG table (e.g., `bbc/bulk.seg`) |
| `--result_dir` | `results` | Output directory for computed CN results |
| `--genome_size` | *(required)* | Reference chromosome sizes file |
| `--region_bed` | *(required)* | Reference chromosome BED file |
| `--mode` | `ilp` | Solver mode: `ilp`, `cd`, or `both` |
| `--solver` | `gurobi` | ILP solver backend: `gurobi` or `cbc` |
| `--timelimit` | None | ILP solver time limit in seconds |
| `--filter_cluster` | False | Enable RDR/BAF variance outlier filtering of clusters before optimization |
| `--filter_std` | 2.0 | Filter clusters whose variance deviates from mean by filter_std × std(variances) |
| `--min_nbins` | 10 | Minimum number of bins a cluster must have to be retained |
| `--ub_nbins` | 50 | Variance-outlier filtering only applies to clusters with #bins ≤ ub_nbins |
| `--balanced_baf_tol` | 0.03 | BAF tolerance for locating balanced clusters |
| `-tR` / `--toleranceRDR` | 0.08 | RDR tolerance for locating clonal CN clusters |
| `-tB` / `--toleranceBAF` | 0.04 | BAF tolerance for locating clonal CN clusters |
| `--minClone` | 2 | Minimum number of tumor clones |
| `--maxClone` | 4 | Maximum number of tumor clones |
| `--diploid` | False | Solve under diploid assumption (cn_max=diploidcmax) |
| `--tetraploid` | False | Solve under tetraploid/WGD assumption (cn_max=tetraploidcmax) |
| `--reg_term` | `MAXCN` | Regularization term: `RAW`, `MAXCN`, `DROOT_SUM`, or `DADJ_SUM` |
| `--reg_steps` | 10 | Number of steps in the regularization path |
| `--reg_stepsize` | 0.01 | Multiplicative step size between regularization path values |
| `--no_ampdel` | False | Disable the amp/del symmetry constraint |
| `--num_cnstates` | -1 | Constrain the number of distinct CN states per clone (-1 = unconstrained) |
| `-eD` / `--diploidcmax` | 6 | Maximum copy-number value for diploid mode (0 = inferred from scaled FCN) |
| `-eT` / `--tetraploidcmax` | 12 | Maximum copy-number value for tetraploid mode (0 = inferred from scaled FCN) |
| `--min_prop` | 0.01 | Minimum clone proportion |
| `--purities` | None | Semicolon-separated `sample:purity` pairs (e.g., `s1:0.80;s2:0.70`); fixes normal-clone proportion to 1 − purity |
| `--cd_niters` | 10 | CD: max outer iterations per seed |
| `--cd_convergence_iters` | 2 | CD: consecutive convergence iterations required to stop |
| `--cd_nseeds` | 400 | CD: number of random restarts |
| `--cd_njobs` | 8 | CD: number of parallel worker processes |
| `--cd_seed` | 42 | CD: random seed for reproducibility |
| `--verbosity` | 0 | Verbose level: 0, 1, or 2 |
