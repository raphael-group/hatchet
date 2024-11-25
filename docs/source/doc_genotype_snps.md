# genotype-snps

Given the normal BAM file, this step of HATCHet2 identifies heterozygous germline SNP positions. The user can restrict candidate positions to a given list (e.g., dbSNP) using the `-R, --snps` argument.

## Input

genotype-snps takes in input a sorted and index BAM file from a matched-normal sample and an indexed human reference genome (preferrably the same version as was used for alignment).

| Name | Description | Usage |
|------|-------------|-------|
| `-N`, `--normal` | A sorted-indexed BAM file | The matched normal sample |
| `-r`, `--reference` | A FASTA file | The human reference genome used for germline variant calling |

## Output

genotype-snps produces a tab-separated VCF file for each chromosome which contains a list of the genomic positions that have been identified as germline heterozygous SNPs in the matched-normal sample.

| Name | Description | Format |
|------|-------------|--------|
| `-l`, `--outputsnps` | the output file for the list of identified heterozygous germline SNPs | `#CHR POS` |

## Main parameters

| Name | Description | Usage | Default |
|------|-------------|-------|---------|
| `-R`, `--snps` | VCF files | Optional list of candidate SNP positions to consider | None* |
| `-st`, `--samtools` | Path to `bin` directory of SAMtools | The path to this direcoty needs to be specified when it is not included in `$PATH` | Path is expected in the enviroment variable `$PATH` |
| `-bt`, `--bcftools` | Path to `bin` directory of BCFtools | The path to this direcoty needs to be specified when it is not included in `$PATH` | Path is expected in the enviroment variable `$PATH` |
| `-c`, `--mincov` | Minimum coverage | Minimum number of reads that have to cover a variant to be called, the value can be increased when considering a dataset with high depth (>60x) | 0 |
| `-C`, `--maxcov` | Maximum coverage | Maximum number of reads that have to cover a variant to be called, the typically suggested value should be twice higher than expected coverage to avoid sequencing and mapping artifacts | 1000 |
| `-j`, `--processes` | Number of parallel jobs | Parallel jobs are used to consider the chromosomes in different samples on parallel. The higher the number the better the running time | 2 |

*When run as a standalone module, the no SNP list is used by default. When run as part of `hatchet run` and no `--snps` argument is supplied, the correct version of dbSNP is automatically downloaded and used to call SNPs.

## Optional parameters

| Name | Description | Usage | Default |
|------|-------------|-------|---------|
| `-v`, `--verbose`  | Verbose logging flag | When enabled, count-alleles outputs a verbose log of the executiong | Not used |
| `-q`, `--readquality` | Threshold for phred-score quality of sequencing reads | The value can be either decreased (e.g. 10) or increased (e.g. 30) to adjust the filtering of sequencing reads | 0 |
| `-Q`, `--basequality` | Threshold for phred-score quality of sequenced nucleotide bases | The value can be either decreased (e.g. 10) or increased (e.g. 30) to adjust the filtering of sequenced nucleotide bases | 11 |
| `-E`,`--newbaq` | Flag to enable `newbaq` veafute of SAMtools | When selected, the user asks SAMtools to recompute alignment of reads on the fly during SNP calling | Not used |
