# plot-cn-1d2d

This step produces alternate plots that show bins in terms of their computed read-depth ratios (RDR), B-allele frequencies (BAF), and assigned copy-number states.
These plots show bins colored by cluster, where the color is consistent between the "2D" (RDR x BAF) view and the "1D" (genomic location x RDR/BAF) view.
Additionally, the labeled points in the 2D plots and the black bars in the 1D plots show the *expected* positions of the assigned copy-number states (determined by the mixture proportions and fractional copy number scaling). These indicators can be used to evaluate the consistency of the HATCHet solution.

When `plot_cn = True` is indicated in `hatchet.ini`, both this command and the command [plot-cn](doc_plot_cn.html) will be run.

These plots should be used to review the results of the copy-number assignment step [compute_cn](doc_compute_cn.html) and tune its parameters. They are most useful when there are few samples - when there are many samples, the visualizations from `plot-cn` may be more interpretable.

## Input

plot-bins requires a single tab-separated file as input. This file (e.g., `best.bbc.ucn`, `results.diploid.n3.bbc.ucn.tsv`, ...) is typically produced by [compute_cn](doc_compute_cn.html) and contains the following fields:

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
| `cn_normal` | The copy number state of the normal diploid clone equal to <code>1&#124;1</code> |
| `u_normal` | The normal admixture of the normal diploid cells in the corresponding sample |
| `cn_clone${n}` | The copy number state of the `${n}` tumor clone in the format <code>A&#124;B</code> where `A` and `B` are the two allele-specific copy numbers of the corresponding genomic bin |
| `u_clone${n}` | The clone proportion of the `${n}` tumor clone in the corresponding sample |

## Parameters

| Name | Description | Usage | Default |
|------|-------------|-------|---------|
| `-O`, `--outdir` | Output directory | Directory in which output plots will be stored | "" (current directory) |
| `--baflim` | Axis limits for mirrored BAF | Specify limits as comma-separated values, e.g., '0,0.51' | None (determined by extremes of data) |
| `--fcnlim` | Axis limits for fractional copy numbers | Specify limits as comma-separated values, e.g., '0,3' | None (determined by extremes of data) |
| `--centromeres` | Mark centromere locations with grey rectangles | Enable by including this flag with no argument. | False |
| `--bysample` | Split plots into distinct sample-specific files (instead of grouping all 1D plots together and all 2D plots together) | Enable by including this flag with no argument. | False |
