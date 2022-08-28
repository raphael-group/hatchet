# Command for running the HATCHet workflow

The entire end-end `HATCHet` pipeline can be run by using the `hatchet run` command. This command requires an *ini-file*
from which it gets its configuration values. A sample
[hatchet.ini](https://raw.githubusercontent.com/raphael-group/hatchet/master/script/hatchet.ini) file is provided in
this folder for you to get started. You can name this file anything you want and specify it during `hatchet run`, but we
will use `hatchet.ini` in the writeup below.

## Set variables

Set all variables in `hatchet.ini` with appropriate values. You will likely not need to modify anything at all other
than the paths to the reference genome, paths to the normal and tumor bam files, and unique names for the tumor
samples, all in the `run` section of `hatchet.ini`:

```
reference = "/path/to/reference.fa"
normal = "/path/to/normal.bam"
bams = "/path/to/tumor1.bam /path/to/tumor2.bam"
samples = "Primary Met"
```

Optionally, if you wish to run the HATCHet pipeline only on select chromosome(s), specify their name(s) under the
'chromosomes' key, separated by whitespace. For example:

```
chromosomes = chr21 chr22
```

This can be very useful when trying to validate your pipeline relatively quickly before running it on all chromosomes.
As an example, this should be set to `chr22` for [HATCHet Demo data](https://zenodo.org/record/4046906).
To run the pipeline on all chromosomes, leave the key blank.

```
chromosomes = 
```

## Run HATCHet without phasing

Use the following command To run HATCHet without phasing:

```
hatchet run hatchet.ini
```

As explained above, you can leave all values to their defaults, but you will want to override the `reference`, `normal`,
`bams` and `samples` values in the ini file.

## Run HATCHet with phasing

Running HATCHet with phasing is currently a two part process. It's a little more labor intensive but may produce cleaner
results.

First run `hatchet run hatchet.ini`, but **enable only the first 3 steps** of the HATCHet pipeline in `hatchet.ini`:

```
genotype_snps = True
count_alleles = True
count_reads = True
combine_counts = False
cluster_bins = False
plot_bins = False
compute_cn = False
plot_cn = False
```

After the run finishes, go to the `snps` subdirectory within the output directory (`output/` by default) specified in
`hatchet.ini`. Here you will find a collection of VCF files, one for each chromosome. These must then be phased (e.g.
using the [Michigan Imputation Server](https://imputationserver.sph.umich.edu/index.html#!)), and the location of the
phased VCF file specified in `hatchet.ini` as the `phase` variable under the `combine_counts` section. If you use the
Michigan imputation server:

1. You may have to use `bcftools annotate` to convert between chromosome names (e.g. chr20 -> 20)
2. Results are always returned in hg19 coordinates, so you may need to convert coordinates back to hg38 using e.g.
   Picard's [LiftoverVcf](https://broadinstitute.github.io/picard/command-line-overview.html#LiftoverVcf)
3. The by-chromosome phased VCF files you receive must be combined with the `bcftools concat` command to give HATCHet a
   single phased VCF file.

Also in `hatchet.ini`, under the `combine_counts` section is a `blocklength` parameter, which is the haplotype block
size used for combining SNPs when estimating B-Allele frequencies. This should ideally be 20kb to 50kb. While larger
haplotype block sizes allow you to combine more SNPs, the accuracy of phasing declines with the block size used (e.g.
see this [paper](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1007308) comparing various
phasing methods).

Then, run the HATCHet workflow again using `hatchet run hatchet.ini`, after enabling only the remaining steps of
the HATCHet pipeline. This should have a shorter runtime than when you ran the first 3 steps:

```
genotype_snps = False
count_alleles = False
count_reads = False
combine_counts = True
cluster_bins = True
plot_bins = True
compute_cn = True
plot_cn = True
```
