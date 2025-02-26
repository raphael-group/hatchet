# count-alleles

Given one or more BAM files and lists of heterozygous SNP positions, this step of HATCHet counts the number of reads covering both the alleles of each identified heterozgyous SNP in every tumor sample.

## Input

count-alleles takes in input sorted and indexed BAM files for multiple tumor samples from the same patient, a sorted and index BAM file from a matched-normal sample, and a indexed human reference genome.

| Name | Description | Usage |
|------|-------------|-------|
| `-T`, `--tumors` | A white-space separated list of sorted-indexed BAM files | The tumor samples from the same patient that are jointly analyzed by HATCHet |
| `-N`, `--normal` | A sorted-indexed BAM file | The matched normal sample for the same patient |
| `-L`, `--snps` | VCF files | One or more files listing heterozygous SNP positions |
| `-r`, `--reference` | A FASTA file | The human reference genome used for germline variant calling |

## Output

count-alleles produces three tab-separated files: the first contains the read counts for every genomic bin in every tumor sample, the second contains the read counts for every genomic bin the matched-normal sample, and the third contains a list of the genomic positions that have been identified as germline heterozygous SNPs in the matched-normal sample.

| Name | Description | Format |
|------|-------------|--------|
| `-O`, `--outputnormal` | The output file for the read counts from matched-normal sample | `#SAMPLE  CHR  POS  REF_COUNT  ALT_COUNT` |
| `-o`, `--outputtumors` | The output file for the read counts from the tumor samples | `#SAMPLE  CHR  POS  REF_COUNT  ALT_COUNT` |
| `-l`, `--outputsnps` | the output directory for the list of identified heterozygous germline SNPs | `#CHR POS` |

The format fields are described in the following.

| Field | Description |
|-------|-------------|
| `SAMPLE` | Name of a sample |
| `CHR` | Name of the chromosome |
| `POS` | Genomic position in `CHR` |
| `REF_COUNT` | Number of reads harboring reference allele in `POS` |
| `ALT_COUNT` | Number of reads harboring alternate allele in `POS` |

## Main parameters

count-alleles has some main parameters; the main values of these parameters allow to deal with most of datasets, but their values can be changed or tuned to accommodate the features of special datasets.

| Name | Description | Usage | Default |
|------|-------------|-------|---------|
| `-S`, `--samples` | White-space separater list of a names | The first name is used for the matched-normal sample, while the others are for the tumor samples and they match the same order of the corresponding BAM files | File names are used |
| `-st`, `--samtools` | Path to `bin` directory of SAMtools | The path to this direcoty needs to be specified when it is not included in `$PATH` | Path is expected in the enviroment variable `$PATH` |
| `-bt`, `--bcftools` | Path to `bin` directory of BCFtools | The path to this direcoty needs to be specified when it is not included in `$PATH` | Path is expected in the enviroment variable `$PATH` |
| `-c`, `--mincov` | Minimum coverage | Minimum number of reads that have to cover a variant to be called, the value can be increased when considering a dataset with high depth (>60x) | 8 |
| `-C`, `--maxcov` | Maximum coverage | Maximum number of reads that have to cover a variant to be called, the typically suggested value should be twice higher than expected coverage to avoid sequencing and mapping artifacts | 300 |
| `-j`, `--processes` | Number of parallele jobs | Parallel jobs are used to consider the chromosomes in different samples on parallel. The higher the number the better the running time | 22 |


## Optional parameters

count-alleles has some optional parameters; changes in the default values of these parameters are not expected to have a significant impact but they can be tuned to better fit the given data.

| Name | Description | Usage | Default |
|------|-------------|-------|---------|
| `-v`, `--verbose`  | Verbose logging flag | When enabled, count-alleles outputs a verbose log of the executiong | Not used |
| `-g`, `--gamma` | Level of confidence for selecting germline heterozygous SNPs | This value is the level of confidence used for the binomial model used to assess whether a called SNPs is in fact germline heterozygous | 0.05 |
| `-q`, `--readquality` | Threshold for phred-score quality of sequencing reads | The value can be either decreased (e.g. 10) or increased (e.g. 30) to adjust the filtering of sequencing reads | 20 |
| `-Q`, `--basequality` | Threshold for phred-score quality of sequenced nucleotide bases | The value can be either decreased (e.g. 10) or increased (e.g. 30) to adjust the filtering of sequenced nucleotide bases | 20 |
| `-U`, `--snpquality` | Threshold for phred-score quality of called variants | The value can be either decreased (e.g. 10) or increased (e.g. 30) to adjust the filtering of called variants | 20 |
| `-L`, `--snps` | Path to file of SNPs in the format `#CHR POS` | When provided, only the included genomic positions will be considered for calling germline SNPs. Using well-known lists (e.g. dbSNP) help to significantly speed up this step | Not used, SNPs are called across all genome |
| `-E`,`--newbaq` | Flag to enable `newbaq` veafute of SAMtools | When selected, the user asks SAMtools to recompute alignment of reads on the fly during SNP calling | Not used |
| `-b`, `--maxshift` | Maximum BAF difference from 0.5 | When used, only SNPs with an absolute difference between the BAF and 0.5 below the maximum are selected | Not used |
