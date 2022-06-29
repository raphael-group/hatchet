# plot-bins-1d2d

This step produces alternate plots that show bins in terms of their computed read-depth ratios (RDR), B-allele frequencies (BAF), and cluster assignments.
These plots show bins colored by cluster, where the color is consistent between the "2D" (RDR x BAF) view and the "1D" (genomic location x RDR/BAF) view.
These plots should be used to review the clustering and tune clustering parameters.

When `plot_bins = True` is indicated in `hatchet.ini`, both this command and the command [plot-bins](doc_plot_bins.html) will be run.

----------------------

## Input

plot-bins considers two different inputs which are tab-separated files, one is mandatory and the other is optional (required to show cluster centers):

1. A file of clustered genomic bins, specified with the flag `-b, --bbc` and with the following fields.

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
| `CLUSTER` | The name of the corresponding cluster of the bin in `SAMPLE` |

2. A file of clusters, which is required only by commands and it is specified through the flag `-s`, `--segfile`. The file has the following fields.

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

## Parameters

| Name | Description | Usage | Default |
|------|-------------|-------|---------|
| `-O`, `--outdir` | Output directory | Directory in which output plots will be stored | "" (current directory) |
| `--baflim` | Axis limits for mirrored BAF | Specify limits as comma-separated values, e.g., '0,0.51' | None (determined by extremes of data) |
| `--rdrlim` | Axis limits for read-depth ratio | Specify limits as comma-separated values, e.g., '0,3' | None (determined by extremes of data) |
| `--centers` | Show cluster centers | Requires SEG file provided via `-s, --seg`. Enable by including this flag. | False |
| `--centromeres` | Mark centromere locations with grey rectangles | Enable by including this flag. | False |
| `-a, --alpha` | Opacity value for scatterplots | Valid values are (0, 1]. Smaller values make points more transparent. | 1 |
