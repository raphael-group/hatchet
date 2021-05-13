# hatchet

This step computes allele-specific fractional copy numbers, solves a constrained distance-based simultaneous factorization to compute allele and clone-specific copy numbers and clone proportions, and deploys a model-selection criterion select the number of clone by explicitly considering the trade-off between subclonal copy-number aberrations and whole-genome duplication.
The step offers some parameters to control each of these features.

## Input

This step requires in input two tab-separated files that must have the same prefix `${PRE}` and the only difference is the extension:

1. The name of the first file must be `${PRE}.bbc` and is a tab-separated file describing the RDR and BAF of clustered genomic bins with the following fields

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
| `CLUSTER` | The name of a cluster |

2. The name of the second file must be `${PRE}.seg` and is a tab-separated file describing the RDR and BAF of clusters of genomic regions with the following fields

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

## Output

This step produces 2 files, a BBC result file and a SEG result file, for every considered number of clones and for both the assumption of the presence and absence of a WGD.
For every considered number of clones `${N}` and assuming `${P}` is either equal to `diploid` or `tetraploid`, the steps generates a pair of files as follows:

1. The first BBC file has the following name `results.${P}.n${N}.bbc.ucn.tsv` and is a copy of the BBC input file with `2 (${N} + 1)` additional fields according to the following:

| Field | Description |
|-------|-------------|
| ... | ... |
| `cn_normal` | The copy number state of the normal diploid clone equal to <code>1&#124;1</code> |
| `u_normal` | The normal admixture of the normal diploid cells in the corresponding sample |
| `cn_clone${n}` | The copy number state of the `${n}` tumor clone in the format <code>A&#124;B</code> where `A` and `B` are the two allele-specific copy numbers of the corresponding genomic bin |
| `u_clone${n}` | The clone proportion of the `${n}` tumor clone in the corresponding sample |

2. The second SEG file has the following name `results.${P}.n${N}.seg.ucn.tsv` and is a file of genomic segments with the following fields:

| Field | Description |
|-------|-------------|
| `CHR` | The name of a chromosome |
| `START` | The genomic position that starts the corresponding genomic segment |
| `END` | The genomic position that ends the corresponding genomic segment |
| `SAMPLE` | The name of a sample |
| `cn_normal` | The copy number state of the normal diploid clone equal to <code>1&#124;1</code> |
| `u_normal` | The normal admixture of the normal diploid cells in the corresponding sample |
| `cn_clone${n}` | The copy number state of the `${n}` tumor clone in the format <code>A&#124;B</code> where `A` and `B` are the two allele-specific copy numbers of the corresponding genomic bin |
| `u_clone${n}` | The clone proportion of the `${n}` tumor clone in the corresponding sample |

In addition, the two files for the best solution assuming both the absence or presence of a WGD are copied in the files `chosen.diploid.bbc.ucn`, `chosen.diploid.seg.ucn` and `chosen.tetraploid.bbc.ucn`, `chosen.tetraploid.seg.ucn`, respectively.
Last, hatchet copies the best solution according to the prediction of WGD in the two files `best.bbc.ucn` and `best.seg.ucn`.

## Fractional copy numbers

The first feature of hatchet is the explicit estimation of allele-specific fractional copy numbers.
hatchet performs this estimation based on rigorous criterion stating that allele-specific fractional copy numbers can be obtained from RDR and BAF of clusters of genomic regions.
More specifically, the criterion states that:
- when there is no WGD, the identification of the diploid cluster with copy-number state (1, 1) is sufficient
- when there is a WGD, the identification of the tetraploid cluster with copy-number state (2, 2) and of a second clonal cluster with a copy-number state (A, B) in all tumor clones are sufficient.

As such, hatchet identifies a heuristic to identifies the required clusters and copy-number states.
This heuristic can be controlled by the following parameters:

| Name | Description | Usage | Default |
|------|-------------|-------|---------|
| `-c`, `--clonal` | The required clusters and corresponding copy number states | User can directly specifies the required clusters and the corresponding copy number states to compute the allele-specific fractional copy number states. These must be speficied in the format `IDX-1:A-1:B-1, ..., IDX-M:A-M:B-M` where `IDX-S` is the name of cluster `S` and `(A-S, B-S)` is the corresponding copy-number state. Moreover, user can specify addittional clusters and copy numbers beyond the required ones. The copy-numbers for these clusters will be fixed during the computation. This can be an usefule feature for especially noisy datasets. | None |
| `-ts`, `--minsize` | Threshold for size of clusters | The minimum size of the clusters to consider for the heuristic that identifies clonal clusters. The non-selected clusters will not be considered as potential tumor-clonal clusters. The threshold must be expressed as a fraction of the entire genome. | 0.02, e.g. `2%` of genome |
| `-tc`, `--minchrs` | Threshold for number of chromosomes | The minimum number of chromosomes covered by the clusters to consider for the heuristic that identifies clonal clusters. The non-selected clusters will not be considered as potential tumor-clonal clusters. | 1 |
| `-td`, `--maxneutralshift` | Maximum BAF shift allowed for diploid cluster | The maximum expected shift from 0.5 for BAF for a diploid or tetraploid cluster (i.e. with copy-number states (1, 1) or (2, 2)). This threshold is used for two goals to identify the diploid or tetraploid cluster. | 0.1 |
| `-tR`, `--toleranceRDR` | Maximum RDR tolerance | The maximum RDR tolerance used by the heuristic when estimating the position of all clonal clusters | 0.04 |
| `-tB`, `--toleranceBAF` | Maximum BAF tolerance | The maximum BAF tolerance used by the heuristic when estimating the position of all clonal clusters | 0.03 |
| `--merge` | Activate merging of clusters | When activated, the heuristic will merge together clusters that appear to have the same values of RDR and BAF, according to the values below. This procedure can help the heuristic by refining the clustering and merging clusters that are likely to have the same copy-number states and unlikely to be clonal. | False, not used |
| `-mR`, `--mergeRDR` | RDR merging threhsold | The maximum difference in RDR considered by the merging procedure | 0.1 |
| `-mB`, `--mergeBAF` | BAF merging threhsold | The maximum difference in BAF considered by the merging procedure | 0.03 |

## Simultaneous factorization

The second feature corresponds to the simultaneous factorization of allele-specific copt numbers and clone proportions.
hatchet solves a constrained and distance-based variant of the factorization that can be controlled by the following parameters.

| Name | Description | Usage | Default |
|------|-------------|-------|---------|
| `-n`, `--clones` | Interval for number of clones | User can specify the minimum and maximum number of clones to consider | `"2,8"` |
| `-eD`, `--diploidcmax` | Maximum copy number with no WGD | The value of maximum copy number that is considered when assuming no WGD. When `0`-value is specified the maximum copy number is directly inferred from the data by rounding the maximum fractional copy number | 8 |
| `-eT`, `--tetraploidcmax` | Maximum copy number with a WGD | The value of maximum copy number that is considered when assuming there is a WGD. When `0`-value is specified the maximum copy number is directly inferred from the data by rounding the maximum fractional copy number | 8 |
| `-u`, `--minprop` | Minimum clone proportion | In every sample, each clone either is non present (clone proportion equal to 0.0) or has a clone proportion higher than this threshold | 0.03 |
| `-f`, `--noampdel` | Activate clone evolutionary contraints | User can decide whether to enable or not constrained about the evolution of tumor clones. These constrained force each allele to be either amplified or deleted across all tumor clones | Activated |
| `-d`, `--cnstates` | Maximum number of distinct copy-number states per cluster | When enabled, the maximum number of distinct copy-number states per cluster is fixed. This option is deprecated | Not used |

HATCHet implements two methods to solve the constrained and distance-based simultaneous factorization: (1) a integer-linear programming (ILP) and (2) a coordinate-descent method (CD).
These methods can be combined in 3 different modes:
- (0) CD + ILP: the solution found by CD is used to start the ILP. As such, ILP attempts to improve the solution found by CD.
- (1) ILP only
- (2) CD only

The mode can be specified through the flag `--mode` whose default is mode 2.
In addition, the solving methods can be controlled by the following parameters.

| Name | Description | Usage | Default |
|------|-------------|-------|---------|
| `-j`, `--jobs` | Number of parallel jobs | Every job will run on parallel a different seed of the CD algorithm and all the jobs can be used for solving the ILP | Number of available processes on machine |
| `-p`, `--seeds` | Number of CD seeds | Number of random restarts executed by the CD algorithm. In general, the higher the number of restarts, the better the result | 400 |
| `-r`, `--randomseed` | Seed for random generator | Every execution started from the same seed have a deterministic behaviour and can be exactly replicated | None, random seed |
| `-s`, `--timelimit` | Time limit | The time limit, expressed in seconds, is imposed to every step of the CD algorithm or to the whole ILP | None |
| `-m`, `--timelimit` | Memory limit | The memory limit, expressed in megabytes, is imposed to every step of the CD algorithm or to the whole ILP. The execution will not be interrupted when reaching the threshold but disk is used | None |
| `--maxiterations` | Maximum number of iteration per seed | This number is imposed as the maximum number of iterations executed for every restart of the CD algorithm | 40 |
| `--diploid` | Assume no WGD | When enabled, HATCHet assumes the absence of a WGD | Not used |
| `--tetraloid` | Assume a WGD | When enabled, HATCHet assumes the occurrence of a WGD | Not used |

## Model selection

This steps have two main parameters to control the model-selection criterion:

| Name | Description | Usage | Default |
|------|-------------|-------|---------|
| `-l`, `--limitinc` | Sensitivity level | The sensitivity level is used to control the confidence in evaluating the presence of tumor clones characterized by small CNAs. By decreasing the value of the sensitivity, HATCHet is more sensible to the presence of small CNAs and small clusters or with small shifts in RDR/BAF are more likely considered as the signal of an additional tumor clone. The possible values of this parameter are between 1.0 and 0.0 and specifically corresponds to an upper bound for the left relative improvement of the objective function. | None, reasonable values to use can be 0.6, 0.5, 0.4, 0.3, 0.2, ... according to the values of the objective function. |
| `-g`, `--ghostprop` | Confidence in the presence of a single tumor clone | This value expresses the confidence of HATCHet when evaluating the presence of a single tumor clone. The higher the value the more likely the presence of a single clone is considered | 0.2 |

## Additional parameters
