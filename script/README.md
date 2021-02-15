# Scripts for running the HATCHet workflow

## 1) Set environmental variables

Set all variables in `config.txt` with appropriate values. You likely will not need to change anything in the `*.sh` scripts.


## 2a) HATCHet without phasing

Use the following commend To run HATCHet without phasing:

```
bash runUnphased.sh > out.txt 2> err.txt
```

Feel free to rename the standard out (out.txt) and standard error (err.txt) files to whatever you wish.

## 2b) HATCHet with phasing

Running HATCHet with phasing is currently a two part process. It's a little more labor intensive on the user end but may produce cleaner results.

First run `01_runPhased.sh`, which executes the first three steps of HATCHet:
```
bash 01_runPhased.sh > out1.txt 2> err1.txt
```

After this script finishes, go to the `snps` subdirectory within the working directory given to HATCHet in `config.txt`. Here you will find a collection of VCF files, one for each chromosome. Give these to an imputation server (e.g. [Michigan Imputation Server](https://imputationserver.sph.umich.edu/index.html#!)), download the phased VCF files, and combine these with the `bcftools concat` command. Give HATCHet the location of this combined phased VCF file in `config.txt` under the `PHASE` variable, meaning you have to run `01_runPhased.sh` and phase the output in the `snps` directory before doing so.

Then, run the second half of the HATCHet workflow, which should have a shorter runtime than the first part:

```
bash 02_runPhased.sh > out2.txt 2> err2.txt
```
