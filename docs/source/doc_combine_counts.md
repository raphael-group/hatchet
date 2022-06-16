# combine-counts

This step constructs variable-length bins that ensure that each bin has at least some number (`--msr`) of SNP-covering reads and at least some number (`--mtr`) of total reads per bin. Then, it combines the read counts and the allele counts for the identified germline SNPs to compute the read-depth ratio (RDR) and B-allele frequency (BAF) of every genomic bin.

## Input

`combine-counts` takes in input the output from `count-reads` (i.e., for each chromosome `ch`, two gzipped filed `ch.total.gz` and `ch.thresholds.gz`). Use the `-A, --array` flag to specify a directory containing these input files.

It also requires (specified by the flag `-b`, `--baffile`) a tab-separated file specifying the allele counts for heterzygous germline SNPs from all tumor samples. The tab separated file would typically be produced by the `count-alleles` command and has the following fields:

| Field | Description |
|-------|-------------|
| `CHR` | Name of a chromosome |
| `POS` | Genomic position corresponding to a heterozygous germline in `CHR` |
| `SAMPLE` | Name of a tumor sample |
| `REF_COUNT` | Count of reads covering `POS` with reference allele |
| `ALT_COUNT` | Count of reads covering `POS` with alternate allele |

Finally, `combine-counts` requires a TSV file (`-t, --totalcounts`) specifying the total number of reads aligned in each sample (also typically produced by `count-reads`).

In summary, **the following arguments are required to specify input**:

| Name | Description | Usage | Default |
|------|-------------|-------|---------|
| `-A`, `--array`  | Directory containing intermediate files | Typically populated by `count-reads`. For each chromosome `ch`, this directory should contain files `ch.total.gz` and `ch.thresholds.gz` (as well as `samples.txt` indicating sample names) |  |
| `-b, --baffile`  | Tab-separated file with allele counts | Typically produced by `count-alleles`. See description above. |  |
| `-t, --totalcounts`  | Tab-separated file with total aligned reads for each sample | Typically produced by `count-alleles`. | |
| `-V, --refversion` | Reference genome version | Either "hg19" or "hg38". This argument is used to select which centromere locations to use. |


## Output

combine-counts produces a tab-separated file (`-o, --outfile`) with the following fields.

| Field | Description |
|-------|-------------|
| `CHR` | Name of a chromosome |
| `START` | Starting genomic position of a genomic bin in `CHR` |
| `END` | Ending genomic position of a genomic bin in `CHR` |
| `SAMPLE` | Name of a tumor sample |
| `RD` | RDR of the bin in `SAMPLE` (corrected by the total reads in `SAMPLE` vs. the total reads in the matched normal sample) |
| `#SNPS` | Number of SNPs present in the bin in `SAMPLE` |
| `COV` | Average coverage in the bin in `SAMPLE` |
| `ALPHA` | Alpha parameter related to the binomial model of BAF for the bin in `SAMPLE`, typically total number of reads from A allele |
| `BETA` | Beta parameter related to the binomial model of BAF for the bin in `SAMPLE`, typically total number of reads from B allele |
| `BAF` | BAF of the bin in `SAMPLE` |
| `TOTAL_READS` | Total number of reads in the bin in `SAMPLE` |
| `NORMAL_READS` |  Total number of reads in the bin in the matched normal sample |
| `CORRECTED_READS` |  Total number of reads in the bin in `SAMPLE`, corrected by the total reads in `SAMPLE` vs. the total reads in matched normal. |

Currently, it produces one such file that excludes sex chromosomes (for use in HATCHet), and one that includes sex chromosomes (for future use).

## Main parameters

combine-counts has some main parameters; the main values of these parameters allow to deal with most of datasets, but their values can be changed or tuned to accommodate the features of special datasets.

| Name | Description | Usage | Default |
|------|-------------|-------|---------|
| `--msr`  | Minimum SNP-covering reads for each bin | Each bin constructed by this command must have at least this many reads covering heterozygous SNPs in each sample | 5000 |
| `--mtr`  | Minimum total reads for each bin | Each bin constructed by this command must have at least this many total reads in each sample | 5000 |

## Optional parameters

| Name | Description | Usage | Default |
|------|-------------|-------|---------|
| `-j`, `--processes` | Number of parallel processes to use (default 1) |  | 1 |
| `--use_mm`  | Use MM BAF inference | If this flag is set, an exhaustive "maximize-maximize" approach inference is used to infer BAF and phasing for each bin (instead of EM). MM results are higher-likelihood and may be closer to the inference from earlier versions of HATCHet, but often produces poorer clusters. |  |
| `-z, --not_compressed`  | Indicates that intermediate files are not compressed | For compatibility with legacy versions of previous step -- set this flag if your `.total` and `.thresholds` files are plaintext rather than gzipped. |  |
| `-p`, `--phase`  | vcf.gz with phasing for all het. SNPs | File containing phasing data for germline SNPs, typically `phased.vcf.gz` if using the HATCHet pipeline. |  |
| `-s`, `--blocksize`  | Maximum phasing block size | Maximum distance (in bp) between a pair of SNPs included in the same phasing block (ignored if `-p, --phase` is not used) | 25000 |
| `-m`, `--max_spb`  | Maximum number of SNPs per phased block | No more than this many SNPs can be included in the same phasing block (included to minimize phasing errors in high-LD regions) | 10 |
| `-a`, `--alpha`  | Significance threshold to allow adjacent SNPs to be merged | If adjacent SNPs have significantly different BAFs (at this significance level) after taking the phasing into account, they are not merged a priori. Higher means less trust in phasing. | 0.1 |

## Example usage
`hatchet combine-counts -b baf/bulk.1bed -o abin/bulk.bb -j 24 -V hg19 -A array -t array/total.tsv -V hg19`
