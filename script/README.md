# Scripts for running the HATCHet workflow

First set all environmental variables in `config.txt` with appropriate values. You do not need to change anything in the `*.sh` scripts.

To run HATCHet without phasing, simply run the following command:

```
bash runUnphased.sh > out.txt 2> err.txt
```

Feel free to rename the standard out (out.txt) and standard error (err.txt) files to whatever you wish.


Running HATCHet **with** phasing is currently a two part process. First run `01_runPhased.sh`, which executes the first three steps of HATCHet:
```
bash 01_runPhased.sh > out1.txt 2> err1.txt
```

After this script finishes, go to the `snps` subdirectory within the working directory given to HATCHet in `config.txt`. Here you will find a collection of VCF files, one for each chromosome. Give these to an imputation server (e.g. [Michigan Imputation Server](https://imputationserver.sph.umich.edu/index.html#!)), download the phased VCF files, and combine these with the `bcftools concat` command. Give HATCHet the location of this combined phased VCF file in `config.txt` under the `PHASE` variable.

Then, run the second half of the HATCHet workflow, which should have a shorter runtime than the first part:

```
bash 02_runPhased.sh > out2.txt 2> err2.txt
```
