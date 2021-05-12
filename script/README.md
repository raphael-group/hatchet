# Scripts for running the HATCHet workflow

## 1) Set environmental variables

Set all variables in `hatchet_config` with appropriate values. You likely will not need to change anything in the `*.sh` scripts.

Note that if the `LIST` variable is blank, HATCHet automatically downloads a list of known germline SNP positions based on your reference genome. If your internet connection is unstable (or if the SNPcaller step is failing), please specify the path to a locally stored list in VCF format.

## 2a) HATCHet without phasing

Use the following command To run HATCHet without phasing:

```
bash hatchet_unphased > out.txt 2> err.txt
```

Feel free to rename the standard out (out.txt) and standard error (err.txt) files to whatever you wish.

## 2b) HATCHet with phasing

Running HATCHet with phasing is currently a two part process. It's a little more labor intensive on the user end but may produce cleaner results.

First run `hatchet_phased_1`, which executes the first three steps of HATCHet:
```
bash hatchet_phased_1 > out1.txt 2> err1.txt
```

After this script finishes, go to the `snps` subdirectory within the working directory given to HATCHet in `hatchet_config`. Here you will find a collection of VCF files, one for each chromosome. These must then be phased (e.g. [Michigan Imputation Server](https://imputationserver.sph.umich.edu/index.html#!)), and the location of the phased VCF file is specified in `hatchet_config` under the `PHASE` variable. If you use the Michigan imputation server: 

1. you may have to use `bcftools annotate` to convert between chromosome names (e.g. chr20 -> 20)
2. results are always returned in hg19 coordinates, so you may need to convert coordinates back to hg38 using e.g. Picard's [LiftoverVcf](https://broadinstitute.github.io/picard/command-line-overview.html#LiftoverVcf)
3. the by-chromosome phased VCF files you receive must be combined with the `bcftools concat` command to give HATCHet a single phased VCF file.

Also in the `hatchet_config` file is the `BLOCK` parameter, which is the haplotype block size used for combining SNPs when estimating B-Allele frequencies. This should ideally be 20kb to 50kb. While larger haplotype block sizes allow you to combine more SNPs, the accuracy of phasing declines with the block size used (e.g. see this [paper](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1007308) comparing various phasing methods).

Then, run the second half of the HATCHet workflow, which should have a shorter runtime than the first part:

```
bash hatchet_phased_2 > out2.txt 2> err2.txt
```
