# Scripts for running the HATCHet workflow

## 1) Set environmental variables

Set all variables in `config.sh` with appropriate values. You likely will not need to change anything in the `*.sh` scripts.


## 2a) HATCHet without phasing

Use the following commend To run HATCHet without phasing:

```
bash runUnphased.sh > out.txt 2> err.txt
```

Feel free to rename the standard out (out.txt) and standard error (err.txt) files to whatever you wish.

## 2b) HATCHet with phasing

Running HATCHet with phasing is currently a two part process. It's a little more labor intensive on the user end but may produce cleaner results.

First run `runPhased_01.sh`, which executes the first three steps of HATCHet:
```
bash runPhased_01.sh > out1.txt 2> err1.txt
```

After this script finishes, go to the `snps` subdirectory within the working directory given to HATCHet in `config.sh`. Here you will find a collection of VCF files, one for each chromosome. These must then be phased (e.g. [Michigan Imputation Server](https://imputationserver.sph.umich.edu/index.html#!)), and the location of the phased VCF file is specified in `config.sh` under the `PHASE` variable. If you use the Michigan imputation server: 

1. you may have to use `bcftools annotate` to convert between chromosome names (e.g. chr20 -> 20)
2. results are always returned in hg19 coordinates, so you may need to convert coordinates back to hg38 using e.g. Picard's [LiftoverVcf](https://broadinstitute.github.io/picard/command-line-overview.html#LiftoverVcf)
3. the by-chromosome phased VCF files you receive must be combined with the `bcftools concat` command to give HATCHet a single phased VCF file.

Then, run the second half of the HATCHet workflow, which should have a shorter runtime than the first part:

```
bash runPhased_02.sh > out2.txt 2> err2.txt
```
