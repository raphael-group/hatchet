# combine-counts-fw

NOTE: This function (formerly called `comBBo`) uses the legacy fixed-width binning described in the HATCHet paper. We recommend using `count-reads` and `combine-counts` which apply an adaptive binning scheme to ensure that each genomic bin has comparable SNP signal.

This step combines the read counts and the allele counts for the identified germline SNPs to compute the read-depth ratio (RDR) and B-allele frequency (BAF) of every genomic bin.

## Input

combine-counts-fw takes in input three tab-separate files:

1. A file of read counts for genomic bins obtained from the matched-normal sample, specified by the flag `-c`, `--normalbins`. The tab separated file has the following fields.

| Field | Description |
|-------|-------------|
| `SAMPLE` | Name of the matched-normal sample |
| `CHR` | Name of a chromosome |
| `START` | Starting genomic position of a genomic bin in `CHR` |
| `END` | Ending genomic position of a genomic bin in `CHR` |
| `COUNT` | Count of sequencing reads in the corresponding bin  |

2. A file of read counts for genomic bins obtained from all the tumor samples, specified by the flag `-C`, `--tumorbins`. The tab separated file has the following fields.

| Field | Description |
|-------|-------------|
| `SAMPLE` | Name of a tumor sample |
| `CHR` | Name of a chromosome |
| `START` | Starting genomic position of a genomic bin in `CHR` |
| `END` | Ending genomic position of a genomic bin in `CHR` |
| `COUNT` | Count of sequencing reads in the corresponding bin |

3. A file of allele counts for heterozygous germline SNPs obtained from all the tumor samples, specified by the flag `-B`, `--tumorbafs`. The tab separated file has the following fields.

| Field | Description |
|-------|-------------|
| `SAMPLE` | Name of a tumor sample |
| `CHR` | Name of a chromosome |
| `POS` | Genomic position corresponding to a heterozygous germline in `CHR` |
| `REF_COUNT` | Count of reads covering `POS` with reference allele |
| `ALT_COUNT` | Count of reads covering `POS` with alternate allele |

## Output

combine-counts-fw produces a tab-separated file with the following fields.

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

## Main parameters

combine-counts has some main parameters; the main values of these parameters allow to deal with most of datasets, but their values can be changed or tuned to accommodate the features of special datasets.

| Name | Description | Usage | Default |
|------|-------------|-------|---------|
| `-d`, `--diploidbaf` | Maximum expected shift from 0.5 for BAF of diploid or tetraploid clusters | The maximum shift is used to identify all the potential clusters with base states (1, 1) or (2, 2). The value depends on the variance in the data (related to noise and coverage); generally, higher variance requires a higher shift. Information provided by plot-bins can help to decide this value in special datasets. | 0.08 (other typically suggested values are 0.1-0.11 for higher variance and 0.06 for low variance) |

## Optional parameters

| Name | Description | Usage | Default |
|------|-------------|-------|---------|
| `-v`, `--verbose`  | Verbose logging flag | When enabled, combine-counts outputs a verbose log of the executiong | Not used |
| `-r`, `--disablebar` | Disabling progress-bar flag | When enabled, the output progress bar is disabled | Not used |
| `-b`, `--normalbafs` | File of allele counts for SNPs in matched-normal sample | When provided, combine-counts attempts to correct the estimated BAF using the variance in matched-normal sample. | Not used (deprecated) |
| `-d`, `--diploidbaf` | Maximum expected shift from 0.5 for BAF of diploid or tetraploid clusters | The maximum shift is used to identify potential potential bins with base states (1, 1) or (2, 2) whose BAF needs to be corrected. The value depends on the variance in the data (related to noise and coverage); generally, higher variance requires a higher shift. Information provided by plot-bins can help to decide this value in special datasets. | 0.08 (other typically suggested values are 0.1-0.11 for higher variance and 0.06 for low variance) |
| `-t`, `--totalcounts` | File of total read counts | When provided, the total read counts are used to normalize the read counts from the corresponding sample | Not used (deprecated) |
| `-m`, `--mode` | Mode used to estimate BAFs | Different modes are provided to combine to combine the allele counts of SNPs | The counts from the allele with minor count are combined (deprecated) |
