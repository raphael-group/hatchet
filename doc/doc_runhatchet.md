# Tutorial

This tutorial illustrates how to use the complete pipeline which is encoded in the main script [runHATCHet.sh](../script/runHATCHet.sh) in `script` folder.
The tutorial is subdivided into some subsections and each of these describes sequential parts of the full pipeline:
1. [Preliminaries](#preliminaries)
2. [Setting up running directory](#rundir)
3. [binBAM](#binbam)
4. [deBAF](#debaf)
5. [comBBo](#combbo)
6. [cluBB](#clubb)
7. [BBot](#bbot)
8. [solve](#solve)
9. [BBeval](#bbeval)

We suggest to make a copy of the script, place the script into the designated running directory, and follow the tutorial.

## Preliminaries
<a name="preliminaries"></a>

```shell
REF="/path/to/reference.fa"
SAM="/path/to/samtools-home/bin/"
BCF="/path/to/bcftools-home/bin/"

XDIR="/path/to/running-dir/"
NORMAL="/path/to/matched-normal.bam"
BAMS="/path/to/tumor-sample1.bam /path/to/tumor-sample2.bam"
ALLNAMES="Normal Primary Met"
NAMES="Primary Met"
J=22

set -e
set -o xtrace
PS4='\''[\t]'\'
export PATH=$PATH:${SAM}
export PATH=$PATH:${BCF}
#source /path/to/virtualenv-python2.7/bin/activate
```

This preliminary part of the script contains all the preliminary information that are required to execute the full pipeline.

***

```shell
REF="/path/to/reference.fa"
SAM="/path/to/samtools-home/bin/"
BCF="/path/to/bcftools-home/bin/"
```

First, one needs to specify the full path to the reference genome with the variable `${REF}`, to the home directory of SAMtools with the variable `${SAM}`, and to the home directory of BCFtools with the variable `${BCF}`. Simply, substitute the value of each variable with the corresponding full path between double apices `"..."`.

***

```shell
XDIR="/path/to/running-dir/"
NORMAL="/path/to/matched-normal.bam"
BAMS="/path/to/tumor-sample1.bam /path/to/tumor-sample2.bam"
ALLNAMES="Normal Primary Met"
```

Second, the user needs to specify the full paths to the running directory where all the files/directories produced by this run of HATCHet will be made with variable `${XDIR}` and the full paths to the required data:
1. `${NORMAL}` is the full path to the BAM file of matched-normal samples
2. `${BAMS}` is a white-space separated list of the BAM files for the multiple tumor samples from the considered patient.

Moreover, `${ALLNAMES}` is a white-space separated list of sample's names where the first is the name of the matched-normal sample and this is follows by the tumor-sample names in the same order as specified in `${BAMS}`.
The variable `${NAMES}` is also a white-space separated list of tumor-sample names and is generally equal to `${ALLNAMES}` but without the matched-normal sample name.
Last, `${J}` is the maximum number of threads that the execution can use. Typically, we suggest to use at least 22 threads if possible in order to consider each chromosome in parallele in the pre-processing steps.

***

```shell
set -e
set -o xtrace
PS4='\''[\t]'\'
export PATH=$PATH:${SAM}
export PATH=$PATH:${BCF}
#source /path/to/virtualenv-python2.7/bin/activate
```

Third, three commands activate the log trace for the script which terminates in case of error and add time stamps to this.
The next two commands add the paths to SAMtools and BCFtools directly in `${PATH}`; this is optional and the paths could be explicitly specified in the corresponding steps.
Last, there is an example command which activates a corresponding python environment.
The usage of virtual environment or anaconda's environments (which command would replace this) is recommended.

***

## Setting up running directory
<a name="rundir"></a>

```shell
BIN=${XDIR}bin/
mkdir -p ${BIN}
BAF=${XDIR}baf/
mkdir -p ${BAF}
BB=${XDIR}bb/
mkdir -p ${BB}
BBC=${XDIR}bbc/
mkdir -p ${BBC}
ANA=${XDIR}analysis/
mkdir -p ${ANA}
RES=${XDIR}results/
mkdir -p ${RES}
EVA=${XDIR}evaluation/
mkdir -p ${EVA}

cd ${XDIR}
```

Next, the script prepare the directories that are used to organize all the outputs from this run of HATCHet.
To avoid condlicts, user should make sure the running directory is an empty directory, expect when re-executing the same run.

## binBAM
<a name="binbam"></a>

```shell
\time -v python2 -m hatchet binBAM -N ${NORMAL} -T ${BAMS} -S ${ALLNAMES} \
                                   -b 50kb -g ${REF} -j ${J} \
                                   -q 20 -O ${BIN}normal.bin -o ${BIN}bulk.bin -v &> ${BIN}bins.log
```

binBAM computes the read counts in genomic bins of the reference genome from the matched-normal sample `${NORMAL}` and from all tumor samples `${BAMS}` with names `${ALLNAMES}`.
The size of genomic bins depend on the size of the CNAs that user aims to infer and on the noise in the data.
More specifically, shorter genomic bins allow to infer more refined CNAs, while larger genomic bins allow to better estimate the values of read-depth ratio (RDR) and B-allele frequency (BAF) for every bin.
The standard size 50kb typically represents a good compromise in whole-genome sequencing (WGS) data, while larger size (e.g. `200kb` or `250kb`) are better suited for whole-exome sequencing (WES) where we expect a sparser distribution of germline SNPs.
The reference genome `${REF}` is provided and `${REF}` must be properly indeced such that in the same folder there is a corresponding dictionary with the same name but `.dict` extenion; this is used to know the length of the sequenced chromosomes.
A standard read quality of 20 is considered and simple parameters are specified including: number of parallel threads, output filenames, verbosity of log, and the log filename of this step.

## deBAF
<a name="debaf"></a>

```shell
\time -v python2 -m hatchet deBAF -N ${NORMAL} -T ${BAMS} -S ${ALLNAMES} \
                                  -r ${REF} -j ${J} -q 20 -Q 20 -U 20 -c 4 \
                                  -C 300 -O ${BAF}normal.baf -o ${BAF}bulk.baf -v \
                                  &> ${BAF}bafs.log
```

deBAF call germline SNPs from the matched-normal sample `${NORMAL}` and computes the corresponding allele counts from all tumor samples `${BAMS}` with names `${ALLNAMES}`.
The same reference genome used for the BAM files has to be specified by `${REF}`.
Minimum and maximum thresholds for the read counts of germline SNPs to consider are fixed between 4 and 300.
GATK best practices suggest that the maximum should be at least twice the expected average coverage to avoid mapping artifacts.
Observe that WES generally requires higher thresholds (e.g. 100 and 3000).
Also, an increasing in the minimum threshold allow to improve the quality of the estimated BAF, while a decreasing allow to consider more SNPs.
Several simple parameters are specified including: number of parallel threads, read/allele/variant standard-quality values at 20, output filenames, verbosity of log, and the log filename of this step.

## comBBo
<a name="combbo"></a>

```shell
\time -v python2 -m hatchet comBBo -c ${BIN}normal.bin -C ${BIN}bulk.bin -B ${BAF}bulk.baf -m MIRROR -e 12 > ${BB}bulk.bb
```

comBBo estimates the RDR and BAF from the read counts of all genomic bins from matched-normal sample in `${BIN}normal.bin`, from the read counts of all genomic bins from all tumor samples in `${BIN}bulk.bin`, and germline SNPs allele counts in `${BAF}bulk.baf`.
The standard method `MIRROR` for estimation is used and a random seed of `12` which allows to replicate the analysis is specified.
The output from standard output is correspondingly written in a BB file `${BB}bulk.bb`.

## cluBB
<a name="clubb"></a>

```shell
\time -v python2 -m hatchet cluBB ${BB}bulk.bb -by ${BNPY} -o ${BBC}bulk.seg -O ${BBC}bulk.bbc \
                                               -e 12 -tB 0.04 -tR 0.15 -d 0.08
#\time -v python2 -m hatchet cluBB ${BB}bulk.bb -by ${BNPY} -o ${BBC}bulk.seg -O ${BBC}bulk.bbc \
#                                               -e 12 -tB 0.04 -tR 0.15 -d 0.08 \
#                                               -u 20 -e 12 -dR 0.002 -dB 0.002
```

cluBB globally clusters genomic bins based on RDR and BAF jointly along the genome and across all tumor samples, specified in a BB file `${BB}bulk.bb`.
The home directory of BNPY is specified through `${BNPY}` to perform a Dirichelt-process clustering.
There are 2 main kind of parameters:
- The maximum expected BAF shift `-d` for diploid segments equal to `0.08`.
- Thresholds for clustering refinement `-tB` and `-tR` are used to merge clusters whose difference is no more than these values in all samples.

For all these parameters, these standard values are good for most of the datasets, however very noisy datasets may require higher thresholds.
In these cases, the values can be estimated and tuned by using the informative plots from BBot.
Other simple parameters allow to specify the output files containing the data of clusters in `${BBC}bulk.seg`, the one containing the clusterd genomic bins in `${BBC}bulk.bbc`, and a random seed of `12`.

This set up is generally effective when considering whole-genome sequencing (WGS) data.
However, a bootstraping approach allows to improve and empower the clustering when considering more sparse data as the ones from whole-exome sequeing (WES).
The bootsraping can be controlled by adding the following parameters (as in the similar commented command):

```shell
-u 20 -dR 0.002 -dB 0.002
```

A total of `20` synthetic genomic bins are added only for the clustering by bootstraping each bin and generating RDR and BAF following normal distributions with `0.002` and `0.002`  variances, respectively.

## BBot
<a name="bbot"></a>

```shell
cd ${ANA}
\time -v python2 -m hatchet BBot -c RD --figsize 6,3 ${BBC}bulk.bbc &
\time -v python2 -m hatchet BBot -c CRD --figsize 6,3 ${BBC}bulk.bbc &
\time -v python2 -m hatchet BBot -c BAF --figsize 6,3 ${BBC}bulk.bbc &
\time -v python2 -m hatchet BBot -c BB ${BBC}bulk.bbc &
\time -v python2 -m hatchet BBot -c CBB ${BBC}bulk.bbc &
wait
```

BBot produces informative plots which are described [here](doc_bbot.md).
Many of these plots can be very useful to assess the performance of the various steps of HATCHet, especially in the case of noisy datasets.
The different plots are generated in parallel; this feature can be disabled by removing `&` from the end of each command and removing `wait`.

## solve
<a name="solve"></a>

```shell
cd ${RES}
\time -v python2 -m hatchet solve -i ${BBC}bulk -n2,8 -p 400 -v 3 \
                                  -u 0.03 -r 12 -j ${J} -eD 6 -eT 12
                                  -g 0.35 -l 0.6 |& tee hatchet.log
```

solve computes the fractional copy numbers, factorize these into allele and clone-specific copy numbers and clone proportions, and applies a model-selection step to jointly infer the number of clones and predict a whole-genome duplication (WGD).
This step requires the common prefix `${BBC}bulk` of the cluster and clustered bins files.
Some basic parameter include:
- the interval `2,8`  specifying the minimum and maximum number of clones which can be increased as much as wanted. User can start to use this interval and increasing if the number of clones appear to be higher or close to the maximum.
- the number of restarts `400` for the coordinate-descent method. This value is typically enough to obtain a good solution close to the optimum, however an increase in this value allows to improve the search of a solution. This can be usefule when observing a number of clusters significantly higher than usual.
- Verbose log with a level of `3`, a more limited log can be obtained with a level of `2`.
- Random initializing seed of `12`.
- Number of parallel threads `${J}`

There are also some important parameters which need to be tuned for special or noisy datasets.
The details of these parameters are the following:
- Maximum copy number without a WGD `-eD` and with a WGD `-eT`. The specified values are typically enough to obtain the copy numbers of large CNAs. While higher values increase the number of reasonable solutions, the values can be increased to investigate the presence of large genomic regions with higher copy numbers. In particular, when the value `0` is specified, then the maximum copy numbers are automatically chosen from the fractional copynumbers.
- Minimum clone proportion `-u`. The value depends on the power to identify tumor clones present in samples at low tumor proportions. The default value is reasonable for most of the datasets, but datasets with high noise may need higher values to avoid overfitting. We thus suggest, especially in these special cases, to repeat the hatchet step by considering increasing values of this parameter up to a maximum or until clone proportions are different from this threshold.
- Two parameters are used to control the elbow method during model selection to choose the best solution. The user should repeat this step with different values of these parameters, especially when one wants to investigate the presence of a single tumor clone or the presence of more tumor clones with small different CNAs and small clone proportions. These two parameters are the following:
  - Confidence in the presence of a single tumor clone `g`. This value determines the confidence in having a single tumor clone in several of the given samples. The user can investigate the hypothesis that a single tumor clone is present by increasing the confidence, e.g. using values of `0.4, 0.5, ...`. Vice versa, the confidence is lowered by decreasing these values. For example, one can do that to investigate the presence of tumor clones with small CNAs or small proportions.
  - Sensitivity `l`. This value determines the sensitivity of the method. Lower values correspond to higher sensitivities which allow to consider solutions with a higher number of clones. The user can decrease this value for investigating the presence of more tumor clones, especially tumor clones with small different CNAs or present in small proportions.

## BBeval
<a name="bbeval"></a>

```shell
cd ${EVA}
\time -v python -m hatchet BBeval ${RES}/best.bbc.ucn -rC 10 -rG 1
```

BBeval produces informative and summary plots for the inferred results.
Examples and details of these plots are reported [here](doc_bbeval.md).
