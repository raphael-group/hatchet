# count-reads

This step of HATCHet uses the locations of heterozygous SNPs (called by `count-alleles`) to identify candidate bin thresholds between SNPs. Then, it counts the total number of reads in each sample between each set of candidate thresholds for use in constructing variable-length bins.

## Input

`count-reads` takes in input sorted and indexed BAM files for multiple tumor samples from the same patient, a sorted and index BAM file from a matched-normal sample, an indexed human reference genome, and a 1bed file containing SNP information for this individual (output from `count-alleles` command, normally `baf/bulk.1bed`).

| Name | Description | Usage |
|------|-------------|-------|
| `-T`, `--tumors` | A white-space separated list of sorted-indexed BAM files for tumor samples | The tumor samples from the same patient that are jointly analyzed by HATCHet |
| `-N`, `--normal` | A sorted-indexed BAM file for matched-normal sample | The matched normal sample for the same patient |
| `-b`, `--baffile` | A 1bed file containing locations of heterozygous germline SNPs | Typically, a user would run `count-alleles` to obtain this file. |
| `-V`, `--refversion` | Reference genome version (hg19 or hg38 supported) | |

## Output

`count-reads` writes all output files to a given output directory. For each chromosome `chr`, count-reads produces two gzipped files needed to construct adaptive bins: `chr.threhsolds.gz` and `chr.total.gz`. `count-reads` also produces a tab-separated file `total.tsv` containing the total number of reads in each sample, and a text file `samples.txt` containing the list of sample names.

| Name | Description |
|------|-------------|
| `-O`, `--outdir` | Output directory | Directory in which output will be written to (must already exist before running `count-reads`)

## Main parameters

| Name | Description | Usage |
|------|-------------|-------|
| `-V`, `--refversion` | Reference genome version ("hg19" or "hg38" supported) | |


## Optional parameters

| Name | Description | Usage | Default |
|------|-------------|-------|---------|
| `-S`, `--samples` | White-space separator list of a names | The first name is used for the matched-normal sample, while the others are for the tumor samples and they match the same order of the corresponding BAM files | File names are used |
| `-j`, `--processes` | Number of parallel jobs | Parallel jobs are used to consider the chromosomes in different samples in parallel. The higher the number the better the running time (up to 22 * n_samples) | 22 |
| `-st`, `--samtools` | Path to `bin` directory of SAMtools | The path to this directory needs to be specified when it is not included in `$PATH` | Path is expected in the environment variable `$PATH` |
| `-md`, `--mosdepth` | Path to the `mosdepth` executable | The path to this executable needs to be specified when it is not included in `$PATH` | None (expected on `PATH`) |
| `-tx`, `--tabix` | Path to the `tabix` executable | The path to this executable needs to be specified when it is not included in `$PATH` | None (expected on `PATH`) |
| `-q`, `--readquality` | Minimum mapping quality (MAPQ) for a read to be counted | Values range from 0 to 41, see alignment documentation for details  | 11 |
| `-i`, `--intermediates` | Keep intermediate files | Retain intermediate files (read starts and per-position coverage) that are used to compute arrays for binning | False  (these files are deleted, subsequent arrays are retained instead) |
## Example usage

Given samtools, mosdepth, and tabix are on the PATH, the referenced files are in the current directory, and the intended output directory `array` is present:

`hatchet count-reads -T first_sample.bam second_sample.bam -N normal_sample.bam -S normal tumor1 tumor2 -V hg19 -j 24 -O array -b baf/bulk.1bed`
