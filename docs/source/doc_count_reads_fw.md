# count-reads-fw

NOTE: This function (formerly called `comBBo`) uses the legacy fixed-width binning described in the HATCHet paper. We recommend using [`count-reads`](doc_count_reads.md) and [`combine-counts`](doc_combine-counts.md) which apply an adaptive binning scheme to ensure that each genomic bin has comparable BAF signal.

This step of HATCHet splits the human reference genome into fixed-width bins (i.e., small genomic regions), and computes the number of sequencing reads aligned to each bin from every given tumor samples and from the matched normal sample.

## Input

count-reads-fw takes in input sorted and indexed BAM files for multiple tumor samples from the same patient, a sorted and index BAM file from a matched-normal sample, and a indexed human reference genome.

| Name | Description | Usage |
|------|-------------|-------|
| `-T`, `--tumors` | A white-space separated list of sorted-indexed BAM files for tumor samples | The tumor samples from the same patient that are jointly analyzed by HATCHet |
| `-N`, `--normal` | A sorted-indexed BAM file for matched-normal sample | The matched normal sample for the same patient |
| `-r`, `--reference` | A FASTA file | The human reference genome used for germline variant calling |

## Output

count-reads-fw produces three tab-separated files: the first contains the read counts for every genomic bin in every tumor sample (raw counts not normalized), the second contains the read counts for every genomic bin the matched-normal sample (raw counts not normalized), and the third contains a list of the genomic positions that have been identified as germline heterozygous SNPs in the matched-normal sample.

| Name | Description |
|------|-------------|
| `-O`, `--outputnormal` | The output file for the read counts from matched-normal sample |
| `-o`, `--outputtumors` | The output file for the read counts from the tumor samples |

## Main parameters

count-reads-fw

## Optional parameters
