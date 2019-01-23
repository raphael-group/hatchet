# HATCHet <br/> <sub><u>H</u>olistic <u>A</u>llele-specific <u>T</u>umor <u>C</u>opy-number <u>Het</u>erogeneity</sub> #

HATCHet is an algorithm to infer allele and clone-specific copy-number aberrations (CNAs), clone proportions, and whole-genome duplications (WGD) for several tumor clones jointly from multiple bulk-tumor samples of the same patient or from a single bulk-tumor sample. HATCHet has been designed and developped by Simone Zaccaria in the group of prof. Ben Raphael at Princeton University. The full description of the algorithm, the comparison with previous state-of-the-art methods, and its application on published cancer datasets are described in

[Simone Zaccaria and Ben Raphael, 2018](https://www.biorxiv.org/content/early/2018/12/17/496174)

The novel simulationg framework, MASCoTE, is available at

[MASCoTE](https://github.com/raphael-group/mascote)

The simulated data, the results of all the methods considered in the comparison, the results of HATCHet on the published whole-genome multi-sample tumor sequencing datasets ([(Gundem et al., *Nature*, 2015)](...) and [(Makohon-Moore et al., *Nature genetics*, 2017)](...)) and all the related analyses are available at

[manuscript's data](https://github.com/raphael-group/hatchet-paper)

This repository includes a quick overview of the HATCHet's algorithm and software, detailed instructions for installation and requirements, a script to execute the full pipeline of HATCHet, demos to illustrate how to apply HATCHet on different datasets with different features, a list of current issues, and contacts.

## Contents ##

1. [Overview](#overview)
    - [Algorithm](#algorithm)
    - [Software](#software)
2. [Setup](#setup)
    - [Dependencies](#dependencies)
    - [Compilation](#compilation)
    - [Using Gurobi](#usinggurobi)
    - [Required data](#requireddata)
3. [Usage](#usage)
    - [Full pipeline and tutorial](#fullpipelineandtutorial)
    - [Demos](#demos)
    - [Custom pipelines](#custompipelines)
    - [Detailed steps](#detailedsteps)
    - [Tips and reccomendations](#tipsandreccomendations)
4. [Current issues](#currentissues)
5. [Contacts](#contacts)

## Overview
<a name="overview"></a>

### Algorithm
<a name="algorithm"></a>

![](doc/hatchet-cartoon.png "HATCHet algorithm")

**Overview of HATCHet algorithm.**
**(A)** HATCHet analyzes the read-depth ratio (RDR) and the B-allele frequency (BAF) in bins of the reference genome (black squares) jointly from multiple tumor samples. Here, we show two tumor samples *p* and *q*. **(B)** HATCHet globally clusters the bins based on RDR and BAF along the entire genome and jointly across samples *p* and *q*. Each cluster (color) includes bins with the same copy-number state within each clone present in *p* or *q*. **(C)** HATCHet estimates the fractional copy number of each cluster. If there is no WGD, the identification of the cluster (magenta) with copy-number state _(1, 1)_ is sufficient and RDRs are scaled correspondingly. If a WGD occurs, HATCHet finds the cluster with copy-number state _(2, 2)_ (same magenta cluster) and a second cluster having an identical copy-number state in all tumor clones. **(D)** HATCHet factorizes the allele-specific fractional copy numbers *F^A, F^B* into the allele-specific copy numbers *A, B*, respectively, and the clone proportions *U*. Here there is a normal clone and 3 tumor clones. **(E)** HATCHet's model selection criterion identifies the matrices *A*, *B* and *U* in the factorization while evaluating the fit according to both the inferred number of clones and presence/absence of a WGD. **(F)** Clusters are classified by their inferred copy-number states in each sample. *Sample-clonal clusters* have a unique copy-number state in the sample and correspond to evenly-spaced positions in the scaled RDR-BAF plot (vertical grid lines in each plot). *Sample-subclonal clusters* (e.g. cyan in *p*) have different copy-number states in a sample and thus correspond to intermediate positions in the scaled RDR-BAF plot. *Tumor-clonal clusters* have identical copy-number states in all tumor clones -- thus they are sample-clonal clusters in every sample and preserve their relative positions in scaled-RDR-BAF plots. In contrast, *tumor-subclonal clusters* have different copy-number states in different tumor clones and their relative positions in the scaled RDR-BAF plot varies across samples (e.g. purple cluster).


### Software
<a name="software"></a>

The current implementation of HATCHet is composed of two sets of modules:

(1) The *core* modules of HATCHet are designed to efficiently solve a challenging constrained and distance-based simultaneous matrix factorization which aim to infer allele and clone-specific copy numbers and clone proportins from fractional copy numbers. The module are implemented in C++11 and are included in `src` folder.

(2) The *utility* modules of HATCHet perform several different tasks that are needed to process the raw data, perform steps of the HATCHet's algorithm needed for the factorization, and process the results. These task include reading/calling germinal single-point mutations, counting reads from a BAM file, combining the read counts and other information, segmenting through HATCHet's global approach, plotting very useful information, etc. THese modules are implemented in python2.7 and are included in `util` and `bin` folders.

## Setup
<a name="setup"></a>

The setup process is composed of 4 steps:
1. [Dependencies](#dependencies): the installation of all the required dependecies; this is a one-time process.
2. [Compilation](#compilation): the compilation of HATCHet; this is a one-time process.
3. [Using Gurobi](#usinggurobi): the requirements to run Gurobi; these need to be satisfied before every run.
4. [Required data](#requireddata): the requirements for considered data; these need to be satisfied whenever using new data.

### Dependencies
<a name="dependencies"></a>

#### > Core

The core module of HATCHet is written in C++11 and thus require a modern C++ compiler that is supporting this (GCC >= 4.8.1, or Clang).
In addition, the core module has the following dependencies:
- [CMake](https://cmake.org/) (>= 3.0), this is needed to drive the compilation. OPTIONAL: we suggest the use of CCMake to facilitate the specifications of other dependencies.
- [Gurobi](http://www.gurobi.com/) (>= 6.0), this is used becase the coordinate-method applied by HATCHet is based on several integer linear programming (ILP) formulations. Gurobi is a commercial ILP solver with two licensing options: (1) a single-host license where the license is tied to a single computer and (2) a network license for use in a compute cluster (using a license server in the cluster). Both options are freely and easily available for users in academia [here](http://www.gurobi.com/academia/academia-center).

#### > Utilities

The utility modules are written in python2.
These modules have the following dependencies, which need to be installed once:

| Name and link | REQUIRED | Usage | Comments |
|---------------|----------|-------|----------|
| subprocess and multiprocess |  REQUIRED | Multiprocessing and threading  | Standard python modules  |
| [SAMtools and BCFtools](http://www.htslib.org/doc/)  | REQUIRED | reading BAM files, read counting, allele counting, and SNP calling | **CURRENT ISSUE**: Unfortunately, HATCHet currently support only the following versions of these softwares: 1.5, 1.6, and 1.7 |
| [BNPY-dev](https://bitbucket.org/michaelchughes/bnpy-dev) | REQUIRED  | Non-parametric Bayesian clustering algorithm | We highly reccomend to install BNPY by simply cloning the repository in a specific directory and verify the related dependencies (Section `Installing *bnpy*` in the [installation instructions](https://bitbucket.org/michaelchughes/bnpy-dev/wiki/Installation.md)). We also suggest to install the additional C++ libraries to improve the performance (Section `Optional: Fast C++ libraries for HMM inference` in the [installation instructions](https://bitbucket.org/michaelchughes/bnpy-dev/wiki/Installation.md)) |
| [pandas](https://pandas.pydata.org/), [matplotlib](https://matplotlib.org/), [seaborn](https://seaborn.pydata.org/) | OPTIONAL | Plotting | Required only for plotting modules. Seaborn has been tested with version 0.8.1. While all these packages can be easily installed using standard package managment tools (e.g. `pip`) for the default python version on your OS, we suggest the usage of python virtual enviroments created directly through python or (reccomended) using [Anaconda](https://www.anaconda.com/). In addition to packages with better perfomance, these enviroments guarantee to keep the same version of all pacjkages without affecting the standard installation of python on your OS. An enviroment specific for HATCHet can be created with the suggested versions of all required packages; the enviroment needs to be activated before every run or explicitly included in the pipeline script. |


### Compilation
<a name="compilation"></a>

Compilation is required to be executed only once for the core modules of HATCHet, while all the other modules are ready to use after the installation of the required packages and dependecies. This repository includes an automatic compilation process which requires only 4 simple steps.

1. **Get [Gurobi](http://www.gurobi.com/)**: user simply needs to unpack Gurobi (Linux version, no other step is required) or to install it (Mac version, standard installation is in `/Library/`).
2. **Build Gurobi**: user only needs to build Gurobi when using newer GXX compilers (>= 5.0) on Linux or having linking issues by following these steps (assuming Gurobi's home is `/path/to/gurobiXXX` where `XXX` is the 3-digit version and `YYYYY64` is equal to either `linux64` or `mac64` when a subfolder is present present):
    ```shell
    $ cd /path/to/gurobiXXX/YYYYY64/src/build/
    $ make
    $ cp libgurobi_c++.a ../../lib/
    ```
3. **Link HATCHet to Gurobi**: substitute the following line at the beginning of `FindGUROBI.cmake`
    ```shell
    set(GUROBI_HOME "" )
    ```
    with
    ```shell
    set(GUROBI_HOME "/path/to/gurobiXXX" )
    ```
    to link HATCHet with Gurobi.
4. **Build HATCHet**: execute the following commands from the root of HATCHet's repository.
    ```shell
    $ mkdir build
    $ cd build/
    $ ccmake ..
    $ make
    ```

When the compilation process fail or when the enviroment has special requirements, the user can manually specify the required paths to Gurobi by following the [detailed intructions](doc/doc_compilation.md).

### Using Gurobi
<a name="usinggurobi"></a>

Every run of HATCHet (especially, the `hatchet` step) needs to use Gurobi which requires a valid license pointed by the enviromental variable `GRB_LICENSE_FILE`. This can be easily obtained depending on the type of free academic license available:

1. **Individual license**. This license can be obtained [easily](http://www.gurobi.com/academia/academia-center) by any academic user with an istitutional email. This license is user and machine-specific, meaning that the user needs to require a different license for every used machine. Assuming the license is stored at `/path/to/gurobi.lic`, the user can easily use it by the following command:
    ```shell
    export GRB_LICENSE_FILE="/path/to/gurobi.lic"
    ```
2. **Multi-use license**. This license can be used by multiple users on any machine in a cluster. This license can be obtained [easily](http://www.gurobi.com/academia/academia-center) but needs to be requested by the IT staff of the user's institution. This license is tipically used in a machine cluster and typically only requires the following command (where the name of gurobi module can slightly change):
    ```shell
    module load gurobi
    ```

As a valid license needs to be available for every run, user can simply add the required command to the pipeline script used to run HATCHet, instead of executing it before every new access to a machine.

### Required data
<a name="requireddata"></a>

HATCHet requires 3 input data:
1. One or more BAM files containing DNA sequencing reads obtained from tumor samples of a single patient. Every BAM file contains the sequencing reads from a sample ad needs to be indexed and sorted. For example, each BAM file can be easily indexed and sorted using [SAMtools](http://www.htslib.org/workflow/#mapping_to_variant). In addition, one can improve the quality of the data by processing the BAM files according to the [GATK Best Practices](https://software.broadinstitute.org/gatk/best-practices/).
2. A BAM file containg DNA sequecing reads obtained from a matched-normal sample of the same patient of the considered tumor samples. The BAM file needs to be indexed and sorted as the tumor BAM files. Also, the BAM files can be processed as the tumor BAM files.
3. A human reference genome. Ideally, one should consider the same human reference genome used to align the sequencing reads in the given BAM files. The most-used human reference genome are available at [GRC](https://www.ncbi.nlm.nih.gov/grc/human) or [UCSC](http://hgdownload.cse.ucsc.edu/downloads.html#human). Observe that human reference genomes use two different notations for chromosomes: either `1, 2, 3, 4, 5 ...` or `chr1, chr2, chr3, chr4, chr5 ...`. One needs to make sure all BAM files and reference genome share that same chromosome notation. When this is not the case, one needs to change the reference to guarantee consistency and needs to re-index the new reference (e.g. using [SAMtools](http://www.htslib.org/workflow/#mapping_to_variant)). Also, HATCHet currently requires the name of the chromosomes in the reference genome is exact and exactly equal to the one in the BAM files, no other labels should be present having `>1 ... \n>2 ... \n>3 ...` or `>chr1 ... \n>chr2 ... \n>chr3`.


## Usage
<a name="usage"></a>

The repository includes all the components that are required to cover every step of the entire HATCHet's pipeline, starting from the processing of raw data reported in a BAM file through the analysis of the final results.
We  provide a script representing the [full pipeline](#fullpipelineandtutorial) of HATCHet and we describe in details the whole script through a tutorial with instructions for usage.
In addition we provide some [demos](#demos) which correspond to guided executions of HATCHet os some small examples and explain in detail the usage of HATCHet when considering standard datasets, real datasets with high noise, and different kind of data.
Moreover, the implementation of HATCHet is highly modular and one can replace any HATCHet's module with any other method to obtain the required results (especially for the pre-processing modules).
As such, we also provide here an overview of the entire pipeline and we describe the [details of each step](#detailedsteps) in a dedicated section of the manual.
Last, we provide some tips and suggestions which allow users to apply HATCHet on datasets with different features.


### Full pipeline and tutorial
<a name="fullpipelineandtutorial"></a>

We provide an exemplary [BASH script](script/runHATCHet.sh) that implements the entire pipeline of HATCHet.
This script and its usage are described in detailed in a guided [tutorial](doc/doc_runhatchet.md).
The user can simply use the script for every execution of HATCHet on different data by copying the script inside the running directory and changing the corresponding paths of the required data and dependecies at the beginning of the script, as described in the guided [tutorial](doc/doc_runhatchet.md).

### Demos
<a name="demos"></a>

Each demo is an exemplary and guided execution of HATCHet on a dataset included in the corresponding demo's folder of this repository (inside `examples`). The demos are meant to illustrate how the user should apply HATCHet on different datasets characterized by different features, noise, and kind of data. In fact, the default parameters of HATCHet allow to succesfully analyze most of the datasets but some of these may be characterized by special features or higher-than-expected variance of the data. Understanding the functioning of HATCHet, assessing the quality of the results, and tuning the few parameters needed to fit the unique features of the considered data thus become crucial to guarantee to always obtain the best-quality results. These are the goals of these demos. More specifically, each demo is simultaneously a guided description of the entire example and a BASH script, which can be directly executed to run the complete demo after setting the few required paths at the beginning of the file. As such, the user can both read the guided description as a web page and run the same script to execute the demo. At this time the following demos are available (more demos will be added in the near future):

| Name | Demo | Folder | Description |
|------|------|--------|-------------|
| `demo-WGS-sim` | [demo-wgs-sim](examples/demo-WGS-sim/demo-wgs-sim.sh) | [demo-WGS-sim](examples/demo-WGS-sim/) | A demo on a typical WGS (whole-genome sequencing) multi-sample dataset with standard noise and variance of the data |
| `demo-WGS-cancer` | [demo-wgs-cancer](examples/demo-WGS-cancer/demo-wgs-cancer.sh) | [demo-WGS-cancer](examples/demo-WGS-cancer/) | A demo on a cancer WGS (whole-genome sequencing) multi-sample dataset with high noise and variance of the data |
| `demo-WES` | [demo-wes](examples/demo-WES/demo-wes.sh) | [demo-WES](examples/demo-WES/) | A demo on a cancer WES (whole-exome sequencing) multi-sample dataset, which is typycally characterized by very high variance of RDR |

### Custom pipelines
<a name="custompipelines"></a>

The repository includes custom pipelines which have been designed to adapt the complete pipeline of HATCHet to special conditions or to integrate the processed data produced by other pipelines. Each custom pipeline is a variation of the main HATCHet's pipeline, we thus recommend the user to always first carefully understand the main [BASH script](script/runHATCHet.sh) through the corresponding guided [tutorial](doc/doc_runhatchet.md) and to carefully understand the provided [demos](#demos) to properly apply HATCHet for best-quality results. Each custom pipeline also includes a specific demo which represent a guided and executable example on example data.

| Name | Description | Script | Demo | Variations |
|------|-------------|--------|------|------------|
| [GATK4-CNV](custom/GATK4-CNV) | Custom pipeline for segmented files from GATK4 CNV pipeline | [custom-gatk4-cnv.sh](custom/GATK4-CNV/custom-gatk4-cnv.sh) | [demo-gatk4-cnv.sh](custom/GATK4-CNV/demo-gatk4-cnv.sh) | This custom pipeline takes the input the segmented files which already contain the estimated RDR and BAF. As such, the first pre-processing steps of HATCHet (`binBAM`, `deBAF`, and `comBBo`) are not needed; for this reason, the following depdencies SAMtools and BCFtools and the following required data, human reference genome, matched-normal sample, and BAM files, are not needed in this case. |

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
  - *Read-count thresholds*. As suggested in the GATK best practices, deBAF requires two parameters -c (the minimum coverage for SNPs) and -C (the maximum coverage for SNPs) to reliably call SNPs and exclude those in regions with artifacts. GATK suggests to consider a value of -C that is at least twice larger than the average coverage and -c should be large enough to exclude non-sequenced regions. For example, `-c of 50 and -C 400` are values previously used but the user should ideally pick values according to the considered data.
  - *Bootstrapping for clustering*. Occasionally, WES may have very few points and much less data points than WGS. In these special cases with very few data points, the global clustering of cluBB may generally benefit from the integrated bootstrapping approach. This approach allow to generate a certain number of synthetic bins from the real ones to increase the power of the clustering. For example, the fllowing cluBB parameters `-u 20 -dR 0.002 -dB 0.002` allow to activate the bootstraping which introduces 20 synthetic bins for each real bin with low variances.
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
