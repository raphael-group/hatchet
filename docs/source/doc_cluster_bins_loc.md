# cluster-bins-loc

This step globally clusters genomic bins along the entire genome and jointly across tumor samples.
`cluster-bins-loc` clusters bins while also taking into account their locations on the genome to form clusters that tend to be contiguous genomic segments
The input/output files for `cluster-bins-loc` are exactly the same as those for `cluster-bins`.

## Input

`cluster-bins-loc` takes in input a tab-separated file with the following fields.

| Field | Description |
|-------|-------------|
| `CHR` | Name of a chromosome |
| `START` | Starting genomic position of a genomic bin in `CHR` |
| `END` | Ending genomic position of a genomic bin in `CHR` |
| `SAMPLE` | Name of a tumor sample |
| `RD` | RDR of the bin in `SAMPLE` |
| `#SNPS` | Number of SNPs present in the bin in `SAMPLE` |
| `COV` | Average coverage in the bin in `SAMPLE` |
| `ALPHA` | Alpha parameter related to the binomial model of BAF for the bin in `SAMPLE`, typically total number of reads from A allele |
| `BETA` | Beta parameter related to the binomial model of BAF for the bin in `SAMPLE`, typically total number of reads from B allele |
| `BAF` | BAF of the bin in `SAMPLE` |

The fields `#SNPS`, `COV`, `ALPHA`, and `BETA` are currently deprecated and their values are ignored.

## Output

`cluster-bins-loc` produces two tab-separated files:

1. A file of clustered genomic bins, specified by the flag `-O`, `--outbins`. The tab separated file has the same fields as the input plus a last field `CLUSTER` which specifies the name of the corresponding cluster.

2. A file of clustered genomic bins, specified by the flag `-o`, `--outsegments`. The tab separated file has the following fields.

| Field | Description |
|-------|-------------|
| `ID` | The name of a cluster |
| `SAMPLE` | The name of a sample |
| `#BINS` | The number of bins included in `ID` |
| `RD` | The RDR of the cluster `ID` in `SAMPLE` |
| `#SNPS` | The total number of SNPs in the cluster `ID` |
| `COV` | The average coverage in the cluster `ID` |
| `ALPHA` | The alpha parameter of the binomial model for the BAF of the cluster `ID` |
| `BETA` | The beta parameter of the binomial model for the BAF of the cluster `ID` |
| `BAF` | The BAF of the cluster `ID` in `SAMPLE` |


## Main parameters

1. `cluster-bins-loc` has a parameter `-d`, `--diploidbaf` that specifies the maximum expected shift from 0.5 for BAF for a diploid or tetraploid cluster (i.e. with copy-number states (1, 1) or (2, 2)). This threshold is used for two goals: (1) To identify the diploid or tetraploid cluster which is used to correct the estimated BAF of potentially biased clusters. (2) To identify potentially biased clusters.
The default value of this parameter (0.1) is typically sufficient for most of the datasets, but its value can be changed or tuned to accommodate the features of special datasets.
In particular, the value of this threshold depends on the variance in the data (related to noise and coverage); generally, higher variance requires a higher shift.
Information provided by `plot-bins` can be crucial to decide whether one needs to change this value in special datasets.

2. By default, `cluster-bins-loc` takes as input a minimum number of clusters (`--minK`, default `2`) and maximum number of clusters (`--maxK`, default `30`), and chooses the number `K` of clusters in this closed interval that maximizes the silhoutette score. Users can also specify an exact number of clusters (`--exactK`) to infer, which skips the model selection step.

3. Other options are available to change aspects of the GHMM model that is used by `cluster-bins-loc`:

| Name | Description | Usage | Default |
|------|-------------|-------|---------|
| `-t`, `--transmat` | Type of transition matrix to infer | `fixed` (to off-diagonal = tau), `diag` (all diagonal elements are equal, all off-diagonal elements are equal) or `full` (freely varying) | `diag` |
| `--tau` | Off-diagonal value for initializing transition matrix | must be `<= 1/(K-1)` | `1e-6` |
| `-c`, `--covar` | Type of covariance matrix to infer | options described in [hmmlearn documentation](https://hmmlearn.readthedocs.io/en/latest/api.html#hmmlearn.hmm.GaussianHMM) | `diag` |
| `-x`, `--decoding` | Decoding algorithm to use to infer final estimates of states | `map` for MAP infernece, `viterbi` for Viterbi algorithm | `map` |

## Optional parameters

| Name | Description | Usage | Default |
|------|-------------|-------|---------|
| `-e`, `--seed`  | Random number generator seed used in model fitting | 0 |
