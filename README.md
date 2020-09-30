![CI](https://github.com/raphael-group/hatchet/workflows/CI/badge.svg)
[![codecov](https://codecov.io/gh/raphael-group/hatchet/branch/master/graph/badge.svg)](https://codecov.io/gh/raphael-group/hatchet)


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
    - [Installation](#installation)
    - [Using Gurobi](#usinggurobi)
    - [Required data](#requireddata)
3. [Usage](#usage)
    - [Full pipeline and tutorial](#fullpipelineandtutorial)
    - [Demos](#demos)
    - [Custom pipelines](#custompipelines)
    - [Detailed steps](#detailedsteps)
    - [Recommendations and quality control](#recommendations)
4. [Current issues](#currentissues)
5. [Contacts](#contacts)

## Overview
<a name="overview"></a>

### Algorithm
<a name="algorithm"></a>

![](doc/hatchet-cartoon.png "HATCHet algorithm")

**Overview of HATCHet algorithm.**
1. HATCHet analyzes the read-depth ratio (RDR) and the B-allele frequency (BAF) in bins of the reference genome (black squares) jointly from multiple tumor samples. Here, we show two tumor samples *p* and *q*.
2. HATCHet globally clusters the bins based on RDR and BAF along the entire genome and jointly across samples *p* and *q*. Each cluster (color) includes bins with the same copy-number state within each clone present in *p* or *q*.
3. HATCHet estimates the fractional copy number of each cluster. If there is no WGD, the identification of the cluster (magenta) with copy-number state _(1, 1)_ is sufficient and RDRs are scaled correspondingly. If a WGD occurs, HATCHet finds the cluster with copy-number state _(2, 2)_ (same magenta cluster) and a second cluster having an identical copy-number state in all tumor clones.
4. HATCHet factorizes the allele-specific fractional copy numbers *F^A, F^B* into the allele-specific copy numbers *A, B*, respectively, and the clone proportions *U*. Here there is a normal clone and 3 tumor clones.
5. HATCHet's model selection criterion identifies the matrices *A*, *B* and *U* in the factorization while evaluating the fit according to both the inferred number of clones and presence/absence of a WGD.
6. Clusters are classified by their inferred copy-number states in each sample. *Sample-clonal clusters* have a unique copy-number state in the sample and correspond to evenly-spaced positions in the scaled RDR-BAF plot (vertical grid lines in each plot). *Sample-subclonal clusters* (e.g. cyan in *p*) have different copy-number states in a sample and thus correspond to intermediate positions in the scaled RDR-BAF plot. *Tumor-clonal clusters* have identical copy-number states in all tumor clones -- thus they are sample-clonal clusters in every sample and preserve their relative positions in scaled-RDR-BAF plots. In contrast, *tumor-subclonal clusters* have different copy-number states in different tumor clones and their relative positions in the scaled RDR-BAF plot varies across samples (e.g. purple cluster).


### Software
<a name="software"></a>

The current implementation of HATCHet is composed of two sets of modules:

(1) The *core* modules of HATCHet are designed to efficiently solve a challenging constrained and distance-based simultaneous matrix factorization which aim to infer allele and clone-specific copy numbers and clone proportins from fractional copy numbers. The module is implemented in C++11 and are included in `src` folder.

(2) The *utility* modules of HATCHet perform several different tasks that are needed to process the raw data, perform steps of the HATCHet's algorithm needed for the factorization, and process the results. These task include reading/calling germinal single-point mutations, counting reads from a BAM file, combining the read counts and other information, segmenting through HATCHet's global approach, plotting very useful information, etc. These modules are implemented in python2.7 and are available as the `util` and `bin` submodules.

## Setup
<a name="setup"></a>

The setup process is composed of 3 steps:
1. [Installation](#installation): the compilation and installation of Hatchet and its dependencies; this is a one-time process.
2. [Using Gurobi](#usinggurobi): the requirements to run Gurobi; these need to be satisfied before every run.
3. [Required data](#requireddata): the requirements for considered data; these need to be satisfied whenever using new data.

### Installation
<a name="installation"></a>

The core module of HATCHet is written in C++11 and thus requires a modern C++ compiler (GCC >= 4.8.1, or Clang).
As long as you have a recent version of GCC or Clang installed, `setuptools` should automatically be able to download a recent version of `cmake` and compile the Hatchet code into a working package.

The installation process can be broken down into the following steps:

1. **Get [Gurobi](http://www.gurobi.com/)** (>= 6.0)

    The coordinate-method applied by HATCHet is based on several integer linear programming (ILP) formulations. Gurobi is a commercial ILP solver with two licensing options: (1) a single-host license where the license is tied to a single computer and (2) a network license for use in a compute cluster (using a license server in the cluster). Both options are freely and [easily available](http://www.gurobi.com/academia/academia-center) for users in academia.
[Download](https://www.gurobi.com/downloads/gurobi-optimizer-eula) Gurobi for your specific platform.

2. **Set GUROBI_HOME environment variable**
    ```shell
    $ export GUROBI_HOME /path/to/gurobiXXX
    ```
    Set `GUROBI_HOME` to where you download Gurobi. Here `XXX` is the 3-digit version of gurobi.

3. **Build Gurobi**
    ```shell
    $ cd "${GUROBI_HOME}"
    $ cd linux64/src/build/
    $ make
    $ cp libgurobi_c++.a ../../lib
    ```
    Substitute `mac64` for `linux64` if using the Mac OSX platform.

4. **Create a new venv/conda environment for Hatchet**

    `Hatchet` is a Python 2.7 package. Unless you want to compile/install it in your default Python 2 environment, you will
want to create either a new Conda environment for Python 2.7 and activate it:
    ```
    conda create --name hatchet python=2.7
    conda activate hatchet
    ```
    or use `virtualenv` through `pip`:
    ```
    python2 -m pip virtualenv env
    source env/bin/activate
    ```

5. **Install basic packages**

    It is **highly recommended** that you upgrade your `pip` and `setuptools` versions to the latest, using:
    ```shell
    pip install -U pip
    pip install -U setuptools
    ```
    
6. **Build and install HATCHet**

    Execute the following commands from the root of HATCHet's repository.
    ```shell
    $ pip install .
    ```

    **NOTE**: If you experience a failure of compilation with an error message like:
    ```
    _undefined reference to symbol 'pthread_create@@GLIBC_2.2.5'_.
    ```

    you may need to set `CXXFLAGS` to `-pthread` before invoking the command:
    ```shell
    $ CXXFLAGS=-pthread pip install .
    ```

    When the compilation process fails or when the environment has special requirements, you may have to manually specify the required paths to Gurobi by following the [detailed intructions](doc/doc_compilation.md).

7. **Install required utilities**

    For reading BAM files, read counting, allele counting, and SNP calling, you need to install [SAMtools and BCFtools](http://www.htslib.org/doc/).
    *Currently, HATCHet support only the following versions of these software: 1.5, 1.6, 1.7*

### Using Gurobi
<a name="usinggurobi"></a>

Every run of HATCHet (especially, the `hatchet` step) needs to use Gurobi which requires a valid license pointed by the enviromental variable `GRB_LICENSE_FILE`. This can be easily obtained depending on the type of free academic license available:

1. **Individual license**. This license can be obtained [easily](http://www.gurobi.com/academia/academia-center) by any academic user with an institutional email. This license is user and machine-specific, meaning that the user needs to require a different license for every used machine. Assuming the license is stored at `/path/to/gurobi.lic`, set the environment variable `GRB_LICENSE_FILE` to point to it:

    ```shell
    export GRB_LICENSE_FILE="/path/to/gurobi.lic"
    ```
2. **Multi-use license**. This license can be used by multiple users on any machine in a cluster. This [license](http://www.gurobi.com/academia/academia-center) can be obtained but needs to be requested by the IT staff of the user's institution. This license is typically used in a machine cluster and only requires the following command (the exact form of this command may vary):
    ```shell
    module load gurobi
    ```

### Required data
<a name="requireddata"></a>

HATCHet requires 3 input data files:
1. One or more BAM files containing DNA sequencing reads obtained from tumor samples of a single patient. Every BAM file contains the sequencing reads from a sample ad needs to be indexed and sorted. For example, each BAM file can be easily indexed and sorted using [SAMtools](http://www.htslib.org/workflow/#mapping_to_variant). In addition, one can improve the quality of the data by processing the BAM files according to the [GATK Best Practices](https://software.broadinstitute.org/gatk/best-practices/).

2. A BAM file containg DNA sequencing reads obtained from a matched-normal sample of the same patient of the considered tumor samples. The BAM file needs to be indexed and sorted as the tumor BAM files. Also, the BAM files can be processed as the tumor BAM files.

3. A human reference genome. Ideally, one should consider the same human reference genome used to align the sequencing reads in the given BAM files. The most-used human reference genome are available at [GRC](https://www.ncbi.nlm.nih.gov/grc/human) or [UCSC](http://hgdownload.cse.ucsc.edu/downloads.html#human). Observe that human reference genomes use two different notations for chromosomes: either `1, 2, 3, 4, 5 ...` or `chr1, chr2, chr3, chr4, chr5 ...`. One needs to make sure all BAM files and reference genome share that same chromosome notation. When this is not the case, one needs to change the reference to guarantee consistency and needs to re-index the new reference (e.g. using [SAMtools](http://www.htslib.org/workflow/#mapping_to_variant)). Also, HATCHet requires that the name of each chromosome is the first word in each ID such that `>1 [ANYTHING] ... \n>2 [ANYTHING] ... \n>3 [ANYTHING] ...` or `>chr1 [ANYTHING] ... \n>chr2 [ANYTHING] ... \n>chr3 [ANYTHING]`.

   For the reference genome, HATCHet requires the existence of a a sequence dictionary (`.dict`), which is part of all standard pipelines for sequencing data, see for example [GATK](https://gatk.broadinstitute.org/hc/en-us/articles/360035531652-FASTA-Reference-genome-format) or [Galaxy](https://galaxyproject.org/admin/data-preparation/). Please note that the sequence dictionary is **NOT** the reference index `.fai`, which is a different structure, has a different function, and it is also recommended.
   
   The dictionary of a reference genome is often included in the available bundles for the reference genomes, see the [example for hg19](ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19) from Broad Institute. However, the dictionary can also be generated in seconds using either [SAMtools](http://www.htslib.org/doc/samtools-dict.html) or [Picard tools](https://gatk.broadinstitute.org/hc/en-us/articles/360036729911-CreateSequenceDictionary-Picard-).
   
   In the folder where you want to download and index the human genome, the steps would typically be:
   
   ```script
   curl -L https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz | gzip -d > hg19.fa
   samtools faidx hg19.fa
   samtools dict hg19.fa > hg19.dict
   ```

## Usage
<a name="usage"></a>

The repository includes all the components that are required to cover every step of the entire HATCHet's pipeline, starting from the processing of raw data reported in a BAM file through the analysis of the final results.
We provide:

- <a name="fullpipelineandtutorial"></a> A script representing the [full pipeline](doc/doc_fullpipeline.md#fullpipelineandtutorial) of HATCHet, and we describe in details the whole script through a tutorial with instructions for usage.
- <a name="demos"></a> [Demos](doc/doc_fullpipeline.md#demos) that correspond to guided executions of HATCHet on some examples, and explain in detail the usage of HATCHet when considering standard datasets, real datasets with high noise, and different kind of data.
- <a name="custompipelines"></a> [Custom pipelines](doc/doc_fullpipeline.md#custompipelines) which adapt the full HATCHet's pipeline to special conditions or integrates pre-processed data belonging to different pipelines.

  The implementation of HATCHet is highly modular and one can replace any HATCHet's module with any other method to obtain the required results (especially for the pre-processing modules).
As such, we also provide here an overview of the entire pipeline and we describe the <a name="detailedsteps"></a> [details of each step](doc/doc_fullpipeline.md#detailedsteps) in a dedicated section of the manual.

- <a name="recommendations"></a> [Recommendations](doc/doc_fullpipeline.md#recommendations), especially for noisy datasets or with different features, to guide the user in the interpretation of HATCHet's inference. We explain how to perform quality control to guarantee the best-quality results, and describe how the user can control and tune some of the parameters to obtain the best-fitting results.

## Current issues
<a name="currentissues"></a>

HATCHet is in active development, please report any issue or question as this could help the devolment and imporvement of HATCHet. Current known issues with current version are reported here below:

- HATCHet currently only supports versions 1.5, 1.6, 1.7 of SAMtools and BCFtools.
- The allele-swapping feature of comBBo has been temporarily disabled due to conflicts with recent SAMtools versions.
- HATCHet has not been tested on Windows yet.

A list of the major recent updates:
- HATCHet accepts now any reference genome (including non humans).
- HATCHet now handles genomic regions properly, recommended for WES data.

## Contacts
<a name="contacts"></a>

HATCHet is actively maintained by Simone Zaccaria, currently a post-doctoral research associate at Princeton University in the research group of prof. Ben Raphael.
