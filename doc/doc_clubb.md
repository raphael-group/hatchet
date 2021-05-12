# cluBB

This step globally clusters genomic bins along the entire genome and jointly across tumor samples, and estimate the corresponding values of RDR and BAF for every cluster in every sample.
cluBB uses BNPY for clustering; the main parameters can be tuned for dealing with special datasets, especially those with high variance or low tumor purity (see Main parameters below).

## Input

cluBB takes in input a tab-separated file with the following fields.

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

cluBB produces two tab-separated files:

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

cluBB has 4 main features with some main parameters that allow  to improve the clustering.

1. cluBB has a parameter `-d`, `--diploidbaf` that specifies the maximum expected shift from 0.5 for BAF for a diploid or tetraploid cluster (i.e. with copy-number states (1, 1) or (2, 2)). This threshold is used for two goals: (1) To identify the diploid or tetraploid cluster which is used to correct the estimated BAF of potentially biased clusters. (2) To identify potentially biased clusters.
The default value of this parameter (0.08) is typically sufficient for most of the datasets, but its value can be changed or tuned to accommodate the features of special datasets.
In particular, the value of this threshold depends on the variance in the data (related to noise and coverage); generally, higher variance requires a higher shift.
Information provided by BBot can be crucial to decide whether one needs to change this value in special datasets.

2. cluBB has some main parameters to control the clustering; the default values for most of these parameters allow to deal with most of datasets, but their values can be changed or tuned to accommodate the features of special datasets.
BBot provides informative plots that can be used to assess the quality of the clustering and evaluate the need of changing some parameters for special datasets.
If your clusters do not appear to be cohesive, try lowering the maximum number of clusters (`-K`) which will force cluBB to infer fewer clusters.

| Name | Description | Usage | Default |
|------|-------------|-------|---------|
| `-K`, `--initclusters` | Maximum number of clusters | The parameter specifies the maximum number of clusters to infer, i.e., the maximum number of GMM components | 50 |
| `-c`, `--concentration` | Concentration parameter for clustering | This parameter determines how much confidence the GMM has in different types of clusterings. Higher values (e.g., 10 or 100)  favor fewer clusters, and smaller values (e.g., 0.01 or 0.001) favor more clusters. For experts, this is the alpha parameter for the Dirichlet process prior. | 1/K |

3. cluBB offers a bootstraping approach that allows a succesfull clustering even when there is a limited number genomic bins that are considred. The bootstraping approach generates sinthetic (i.e. used only for clustering) bins based on the data of the given bins. The bootstraping is controlled by the following parameters.

| Name | Description | Usage | Default |
|------|-------------|-------|---------|
| `-u`, `--bootclustering` | Number of sinthetic bins to generate | Sinthetic bins can be generated based on the RDR and BAF of given bins and are added only to the clustering to improve it when the total number of bins is low (e.g. when considering data from WES) | 0, not used |
| `-dR`,`--ratiodeviation` | Standard deviation for generate RDR of sinthetic bins | The parameter affects the variance of the generated data, this value can be estimated from given bins and BBot generates informative plots to do this | 0.02 |
| `-dB`,`--bafdeviation` | Standard deviation for generate BAF of sinthetic bins | The parameter affects the variance of the generated data, this value can be estimated from given bins and BBot generates informative plots to do this | 0.02 |
| `-s`, `--seed` | Random seed | The value is used to seed the random generation of RDR and BAF of synthetic bins | 0 |

4. cluBB offers a basic iterative process to merge clusters according to given tolerances. This feature can be used to refine the results of BNPY clustering and merge distinct clusters that are not sufficiently distinguished. This process can be controlled by the following parameters.

| Name | Description | Usage | Default |
|------|-------------|-------|---------|
| `-tR`, `--tolerancerdr` | Tolerance for RDR | The value is used to determine when two clusters should be merged in terms of RDR | 0.0, merging is not performed |
| `-tB`, `--tolerancebaf` | Tolerance for BAF | The value is used to determine when two clusters should be merged in terms of BAF | 0.0, merging is not performed |

## Optional parameters

| Name | Description | Usage | Default |
|------|-------------|-------|---------|
| `-v`, `--verbose`  | Verbose logging flag | When enabled, comBBo outputs a verbose log of the executiong | Not used |
| `-r`, `--disablebar` | Disabling progress-bar flag | When enabled, the output progress bar is disabled | Not used |
