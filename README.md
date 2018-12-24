# HATCHet <br/> <sub><u>H</u>olistic <u>A</u>llele-specific <u>T</u>umor <u>C</u>opy-number <u>Het</u>erogeneity</sub> #

HATCHet is an algorithm to infer allele and clone-specific copy-number aberrations (CNAs), clone proportions, and whole-genome duplications (WGD) for several tumor clones jointly from multiple bulk-tumor samples of the same patient or from a single bulk-tumor sample. HATCHet has been designed and developped by Simone Zaccaria in the group of prof. Ben Raphael at Princeton University. The full description of the algorithm, the comparison with previous state-of-the-art methods, and its application on published cancer datasets are described in

[Simone Zaccaria and Ben Raphael, 2018](https://www.biorxiv.org/content/early/2018/12/17/496174)

The novel simulationg framework, MASCoTE, is available at

[MASCoTE](https://github.com/raphael-group/mascote)

The simulated data, the results of all the methods considered in the comparison, the results of HATCHet on the published whole-genome multi-sample tumor sequencing datasets ([(Gundem et al., *Nature*, 2015)](...) and [(Makohon-Moore et al., *Nature genetics*, 2017)](...)) and all the related analyses are available at

[manuscript's data](https://github.com/raphael-group/hatchet-paper)

## Contents ##

1. [Overview](#overview)
    - [Algorithm](#algorithm)
    - [Software](#software)
2. [Setup](#setup)
    - [Dependencies](#dependencies)
    - [Compilation](#compilation)
    - [Required data](#requireddata)
3. [Usage](#usage)
    - [Full pipeline and tutorial](#fullpipelineandtutorial)
    - [Detailed steps](#detailedsteps)
    - [Tips and reccomendations](#tipsandreccomendations)
4. [Current issues](#currentissues)
5. [Contacts](#contacts)

## Overview
<a name="overview"></a>

### Algorithm
<a name="algorithm"></a>

![](doc/hatchet-cartoon.png "HATCHet algorithm")

**Overview of \hatchetFull (\hatchet) algorithm.**
**(A)** HATCHet analyzes the read-depth ratio (RDR) and the B-allele frequency (BAF) in bins of the reference genome (black squares) jointly from multiple tumor samples. Here, we show two tumor samples *p* and *q*. **(B)** HATCHet globally clusters the bins based on RDR and BAF along the entire genome and jointly across samples *p* and *q*. Each cluster (color) includes bins with the same copy-number state within each clone present in *p* or *q*. **(C)** HATCHet estimates the fractional copy number of each cluster. If there is no WGD, the identification of the cluster (magenta) with copy-number state *(1, 1)* is sufficient and RDRs are scaled correspondingly. If a WGD occurs, HATCHet finds the cluster with copy-number state _(2, 2)_ (same magenta cluster) and a second cluster having an identical copy-number state in all tumor clones. **(D)** HATCHet factorizes the allele-specific fractional copy numbers *F^A, F^B* into the allele-specific copy numbers *A, B*, respectively, and the clone proportions *U*. Here there is a normal clone and 3 tumor clones. **(E)** HATCHet's model selection criterion identifies the matrices *A*, *B* and *U* in the factorization while evaluating the fit according to both the inferred number of clones and presence/absence of a WGD. **(F)** Clusters are classified by their inferred copy-number states in each sample. *Sample-clonal clusters* have a unique copy-number state in the sample and correspond to evenly-spaced positions in the scaled RDR-BAF plot (vertical grid lines in each plot). *Sample-subclonal clusters* (e.g. cyan in *p*) have different copy-number states in a sample and thus correspond to intermediate positions in the scaled RDR-BAF plot. *Tumor-clonal clusters* have identical copy-number states in all tumor clones -- thus they are sample-clonal clusters in every sample and preserve their relative positions in scaled-RDR-BAF plots. In contrast, *tumor-subclonal clusters* have different copy-number states in different tumor clones and their relative positions in the scaled RDR-BAF plot varies across samples (e.g. purple cluster).


### Software
<a name="software"></a>

The current implementation of HATCHet is composed of two sets of modules:

(1) The *core* modules of HATCHet are designed to efficiently solve a challenging constrained and distance-based simultaneous matrix factorization which aim to infer allele and clone-specific copy numbers and clone proportins from fractional copy numbers. The module are implemented in C++11 and are included in `src` folder.

(2) The *utility* modules of HATCHet perform several different tasks that are needed to process the raw data, perform steps of the HATCHet's algorithm needed for the factorization, and process the results. These task include reading/calling germinal single-point mutations, counting reads from a BAM file, combining the read counts and other information, segmenting through HATCHet's global approach, plotting very useful information, etc. THese modules are implemented in python2.7 and are included in `util` and `bin` folders.

## Setup
<a name="setup"></a>

### Dependencies
<a name="dependencies"></a>

#### > Core

The core module of HATCHet is written in C++11 and thus require a modern C++ compiler that is supporting this (GCC >= 4.8.1, or Clang).
In addition, the core module has the following dependencies:
- [CMake](https://cmake.org/) (>= 3.0), this is needed to drive the compilation. OPTIONAL: we suggest the use of CCMake to facilitate the specifications of other dependencies.
- [Gurobi](http://www.gurobi.com/) (>= 6.0), this is used becase the coordinate-method applied by HATCHet is based on several integer linear programming (ILP) formulations. Gurobi is a commercial ILP solver with two licensing options: (1) a single-host license where the license is tied to a single computer and (2) a network license for use in a compute cluster (using a license server in the cluster). Both options are freely and easily available for users in academia [here](http://www.gurobi.com/academia/academia-center).

#### > Utilities

The utility modules are written in python2.
These modules have the following dependencies:

| Name and link | REQUIRED | Usage | Comments |
|---------------|----------|-------|----------|
| subprocess and multiprocess |  REQUIRED | Multiprocessing and threading  | Standard python modules  |
| [SAMtools and BCFtools](http://www.htslib.org/doc/)  | REQUIRED | reading BAM files, read counting, allele counting, and SNP calling | **CURRENT ISSUE**: Unfortunately, HATCHet currently support only the following versions of these softwares: 1.5, 1.6, and 1.7 |
| [BNPY-dev](https://bitbucket.org/michaelchughes/bnpy-dev) | REQUIRED  | Non-parametric Bayesian clustering algorithm | We highly reccomend to install BNPY by simply cloning the repository in a specific directory and verify the related dependencies (Section `Installing *bnpy*` in the [installation instructions](https://bitbucket.org/michaelchughes/bnpy-dev/wiki/Installation.md)). We also suggest to install the additional C++ libraries to improve the performance (Section `Optional: Fast C++ libraries for HMM inference` in the [installation instructions](https://bitbucket.org/michaelchughes/bnpy-dev/wiki/Installation.md)) |
| [pandas](https://pandas.pydata.org/), [matplotlib](https://matplotlib.org/), [seaborn](https://seaborn.pydata.org/) | OPTIONAL | Plotting | Required only for plotting modules. Seaborn has been tested with version 0.8.1 |

While all these packages can be easily installed using standard package managment tools (e.g. `pip`) for the default python version on your OS, we suggest the usage of python virtual enviroments created directly through python or (reccomended) using [Anaconda](https://www.anaconda.com/).
In addition to packages with better perfomance, these enviroments guarantee to keep the same version of all pacjkages without affecting the standard installation of python on your OS.
An enviroment specific for HATCHet can be created with the suggested versions of all required packages.

### Compilation
<a name="compilation"></a>

Compilation is only required when starting from source code.
Only the core modules of HATCHet requires compilation, while all the other modules are ready to use.
To perform the compilation, execute the following commands from the root of HATCHet's repository.

```shell
$ mkdir build
$ cd build/
$ ccmake ..
$ make
```

HATCHet's compilation process attempts to automatically find the following Gurobi's paths.

| Name | Path | Comment |
|------|------|---------|
| `GUROBI_CPP_LIB` | `/to/gurobiXXX/YY/lib/libgurobi_c++.a`  | <ul><li>`/to/gurobi` is the path to Gurobi's home, typically `/opt/gurobiXXX` for linux and `/Library/gurobiXXX` for mac</li><li>`XXX` is the Gurobi full version, e.g. 702 or 751</li><li>`YY` depends on os, typically `linux64` for linux or `mac64` for mac</li></ul> |
| `GUROBI_INCLUDE_DIR` | `/to/gurobiXXX/YY/include`  | <ul><li>`/to/gurobi` is the path to Gurobi's home</li><li>`XXX` is the Gurobi full version</li><li>`YY` depends on os</li></ul> |
| `GUROBI_LIB` | `/to/gurobiXXX/YY/lib/libgurobiZZ.so`  | <ul><li>`/to/gurobiXXX` is the path to Gurobi's home</li><li>`XXX` is the Gurobi full version</li><li>`YY` depends on os</li><li>`ZZ` are typically the first 2 numbers of `XXX`</li></ul> |

If the compilation fails to find the Gurobi's paths, these need to be specified directly by either using

```shell
$ ccmake ..
```

or by directly running `CMake` with proper flags as following

```shell
$ cmake .. \
        -DGUROBI_CPP_LIB=/to/gurobiXXX/YY/lib/libgurobi_c++.a \
        -DGUROBI_INCLUDE_DIR=/to/gurobiXXX/YY/include \
        -DGUROBI_LIB=/to/gurobiXXX/YY/lib/libgurobiZZ.so
```

### Required data
<a name="requiredata"></a>

HATCHet requires 3 input data:
1. One or more BAM files containing DNA sequencing reads obtained from tumor samples of a single patient. Every BAM file contains the sequencing reads from a sample ad needs to be indexed and sorted. For example, each BAM file can be easily indexed and sorted using [SAMtools](http://www.htslib.org/workflow/#mapping_to_variant). In addition, one can improve the quality of the data by processing the BAM files according to the [GATK Best Practices](https://software.broadinstitute.org/gatk/best-practices/).
2. A BAM file containg DNA sequecing reads obtained from a matched-normal sample of the same patient of the considered tumor samples. The BAM file needs to be indexed and sorted as the tumor BAM files. Also, the BAM files can be processed as the tumor BAM files.
3. A human reference genome. Ideally, one should consider the same human reference genome used to align the sequencing reads in the given BAM files. The most-used human reference genome are available at [GRC](https://www.ncbi.nlm.nih.gov/grc/human) or [UCSC](http://hgdownload.cse.ucsc.edu/downloads.html#human). Observe that human reference genomes use two different notations for chromosomes: either `1, 2, 3, 4, 5 ...` or `chr1, chr2, chr3, chr4, chr5 ...`. One needs to make sure all BAM files and reference genome share that same chromosome notation. When this is not the case, one needs to change the reference to guarantee consistency and needs to re-index the new reference (e.g. using [SAMtools](http://www.htslib.org/workflow/#mapping_to_variant)). Also, HATCHet currently requires the name of the chromosomes in the reference genome is exact and exactly equal to the one in the BAM files, no other labels should be present having `>1 ... \n>2 ... \n>3 ...` or `>chr1 ... \n>chr2 ... \n>chr3`.


## Usage
<a name="usage"></a>

The repository includes all the components that are required to cover every step of the entire pipeline, starting from the processing of raw data reported in a BAM file through the analysis of the final results.
We  provide a script representing the full pipeline of HATCHet and we describe in details the whole script through a tutorial with instructions for usage.
However, the implementation of HATCHet is highly modular and one can replace any component with any other method to obtain the required results.
As such, we also provide here an overview of the entire pipeline and we describe the details of each step in a dedicted section of the manual.


### Full pipeline and tutorial
<a name="fullpipelineandtutorial"></a>

We provide an exemplary [BASH script](script/runHATCHet.sh) that implements the entire pipeline of HATCHet.
This script and its usage are described in detailed in a guided [tutorial](doc/doc_runhatchet.md).

### Detailed steps
<a name="detailedsteps"></a>

The full pipeline of HATCHet is composed of 7 sequential steps, starting from the required input data.
The descrition of each step also includes the details of the corresponding input/output that are especially useful when one wants to replace or change some of the steps in the pipeline while guaranteeing the correct functioning of HATCHet.

| Order | Step | Description |
|--|--|--|
| (1) | [*binBAM*](doc/doc_binbam.md) | This step splits the human reference genome into bins, i.e. fixed-size small genomic regions, and computes the number of sequencing reads aligned to each bin from every given tumor samples and from the matched normal sample. |
| (2) | [*deBAF*](doc/doc_debaf.md) | This step calls heterozygous germline SNPs from the matched-normal sample and counts the number of reads covering both the alleles of each identified heterozgyous SNP in every tumor sample. |
| (3) | [*comBBo*](doc/doc_combbo.md) | This step combines the read counts and the allele counts for the identified germline SNPs to compute the read-depth ratio (RDR) and B-allele frequency (BAF) of every genomic bin. |
| (4) | [*cluBB*](doc/doc_clubb.md) | This step globally clusters genomic bins along the entire genome and jointly across tumor samples, and estimate the corresponding values of RDR and BAF for every cluster in every sample. |
| (5) | [*BBot*](doc/doc_bbot.md) | This step produces informative plots concerning the computed RDRs, BAFs, and clusters. The information produced by this step are important to validate the compute clusters of genomic regions. |
| (6) | [*hatchet*](doc/doc_hatchet.md) | This step computes allele-specific fractional copy numbers, solves a constrained distance-based simultaneous factorization to compute allele and clone-specific copy numbers and clone proportions, and deploys a model-selection criterion select the number of clone by explicitly considering the trade-off between subclonal copy-number aberrations and whole-genome duplication. |
| (7) | [*BBeval*](doc/doc_bbeval.md) | This step analyzes the inferred copy-number states and clone proportions and produces informative plots jointly considering all samples from the same patient. In addition, this step can also combine results obtained for different patients and perform integrative analysis. |

### Tips and reccomendations
<a name="tipsandreccomendations"></a>

- All the components of HATCHet's pipeline depend on some basic parameters that allow to deal with data characterized by different features. The default values of these parameters allow one to succesfully apply HATCHet on most datasets. However, special or noisy datasets may require some tuning of the parameters. The user can deal with these cases by following the suggestions reported here, reading the descriptions of the various steps, and using the informative plots to verify the results. In particular:
  - *Overfitting of minimum clone proportion*. If the results contain several clone proportions equal to the given minimum clone proportion, try to consider increasing values of this threshold as this mya indicate potential overfitting.
  - *Increase specificity*. HATCHet is based on a criterion of parsimony and tumor clones characterized by small CNAs or very low clone proportions may be missed. If one wants to investigate the presence of such clones (especially with small CNAs), the user can vary the corresponding parameter in [*hatchet*](doc/doc_hatchet.md) (see also [tutorial](doc/doc_runhatchet.md)). Observe that the parameter shall be decreased to increase the sensitivity.
  - *Non-decreasing objective function*. There are two possible explanations for observing values of the objective function for increasing number of clones which do not decrease or decrease only when considering a single tumor clone: (1) there is a single tumor clone, (2) the heuristic of HATCHet identified wrong tumor-clonal clusters. While (1) can be assessed by varying the corresponding parameter in [*hatchet*](doc/doc_hatchet.md) (see also [tutorial](doc/doc_runhatchet.md)), (2) has to be excluded. In particular, the heuristic of HATCHet may fail to identify a second tumor-clonal cluster required with a WGD in very noisy datasets or the ones only comprising low-purity samples. One can verifies the identification by using the informative plots from BBot (e.g. BB command) and can correct potential errors by either (1) changing the parameters of the related heuristic of HATCHet (see [*hatchet*](doc/doc_hatchet.md)) or by manually identifying this second cluster by using the informative plots (the second and additional clusters can be manually specified).
- *Control clustering*. The global clustering is a crucial feature of HATCHet and the quality of the final results is affected by the quality of the clustering. The default parameters allow to deal with most of the datasets, however the user can validate the results and improve it. In particular there are 2 parameters to control the clustering and 2 parameters to refine the clusters (see [*cluBB*](doc/doc_clubb.md)). These parameters can be used to obtain the best result. The user can repeat cluBB with different settings and choose the best results considering the plots, especially BB plots, produced by BBot.
  - *WGS/WES data*. The default values used by HATCHet are for analyzing whole-genome sequencing (WGS) data. However, when considering whole-exome sequencing (WES) data some of the parameters need to be adjusted due to the different features of this kind of data. More specifically, there are 3 main points to consider when analyzing WES data:
    - *Larger bin sizes*. While a size of 50kb is standard for CNA analysis when considering whole-genome sequencing (WGS) data, data from whole-exome sequencing (WES) generally require to use large bin sizes in order to guarantee that each bin contains a sufficient number of heterozygous germline SNPs. Indeed, having a sufficient number of germline SNPs is needed to have good estimations for the B-allele frequency (BAF) of each bin. As such, more appropriate bin sizes to consider may be 200kb or 250k when analyzing WES data. However, one can use the informative plots to test different bin sizes and obtain the smallest size that allows low-variance estimates of BAF.
    - *Read-count thresholds*. As suggested in the GATK best practices, deBAF requires two parameters -c (the minimum coverage for SNPs) and -C (the maximum coverage for SNPs) to reliably call SNPs and exclude those in regions with artifacts. GATK suggests to consider a value of -C that is at least twice larger than the average coverage and -c should be large enough to exclude non-sequenced regions. For example, `-c of 50 and -C 3000` are values previously used but the user should ideally pick values according to the considered data.
    - *Bootstrapping for clustering*. WES has much less data points than WGS. As such, the global clustering of cluBB may generally benefit from the integrated bootstrapping approach. This approach allow to generate a certain number of synthetic bins from the real ones to increase the power of the clustering. For example, the fllowing cluBB parameters `-u 20 -dR 0.002 -dB 0.002` allow to activate the bootstraping which introduces 20 synthetic bins for each real bin with low variances.
- *SNP calling from scratch*. HATCHet allows to provide to deBAF a list of known germline SNPs. This allows to significantly improve the performance. However, running deBAF without this list results in deBAF calling germline SNPs along the genome and allowing to identify private germline SNPs and increase the total number. The user can consider this trade-off.

## Current issues
<a name="currentissues"></a>

HATCHet is in active development, please report any issue or question as this could help the devolment and imporvement of HATCHet. Current known issues with current version are reported here below:

- HATCHet currently only supports versions 1.5 -- 1.6 -- 1.7 of SAMtools and BCFtools
- HATCHet currently only supports dev-bitbucket version of BNPY
- The allele-swapping feature of comBBo has been temporarily disabled due to conflicts with recent SAMtools versions
- Binaries will be available soon

## Contacts
<a name="contacts"></a>

HATCHet's repository is actively mantained by Simone Zaccaria, currently a postodoctoral research associate at Princeton University in the research group of prof. Ben Raphael.
