# binBAM

This step of HATCHet splits the human reference genome into bins, i.e. fixed-size small genomic regions, and computes the number of sequencing reads aligned to each bin from every given tumor samples and from the matched normal sample.

## Input

binBAM takes in input sorted and indexed BAM files for multiple tumor samples from the same patient, a sorted and index BAM file from a matched-normal sample, and a indexed human reference genome.

| Name | Description | Usage |
|------|-------------|-------|
| `-T`, `--tumors` | A white-space separated list of sorted-indexed BAM files for tumor samples | The tumor samples from the same patient that are jointly analyzed by HATCHet |
| `-N`, `--normal` | A sorted-indexed BAM file for matched-normal sample | The matched normal sample for the same patient |
| `-r`, `--reference` | A FASTA file | The human reference genome used for germline variant calling |

## Output

binBAM produces three tab-separated files: the first contains the read counts for every genomic bin in every tumor sample, the second contains the read counts for every genomic bin the matched-normal sample, and the third contains a list of the genomic positions that have been identified as germline heterozygous SNPs in the matched-normal sample.

| Name | Description |
|------|-------------|
| `-O`, `--outputnormal` | The output file for the read counts from matched-normal sample |
| `-o`, `--outputtumors` | The output file for the read counts from the tumor samples |
| `-l`, `--outputsnps` | the output file for the list of identified heterozygous germline SNPs |

## Main parameters

binBAM

## Optional parameters
