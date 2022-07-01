## Full pipeline and tutorial
<a name="fullpipelineandtutorial"></a>

We provide example [BASH scripts](script/README.md) that implement the entire pipeline of HATCHet.
This script and its usage are described in detailed in a guided [tutorial](doc_runhatchet.md).
The user can simply use the script for every execution of HATCHet on different data by copying the script inside the running directory and changing the corresponding paths of the required data and dependencies at the beginning of the script, as described in the guided [tutorial](doc_runhatchet.md).

## Demos
<a name="demos"></a>

Each demo is an example and guided execution of HATCHet on a dataset included in the corresponding demo's folder of this
repository (inside `examples`). The demos are meant to illustrate how the user should apply HATCHet on different
datasets characterized by different features, noise, and kind of data. In fact, the default parameters of HATCHet allow
to successfully analyze most of the datasets but some of these may be characterized by special features or
higher-than-expected variance of the data. Understanding the functioning of HATCHet, assessing the quality of the
results, and tuning the few parameters needed to fit the unique features of the considered data thus become crucial to
guarantee to always obtain the best-quality results. These are the goals of these demos.

More specifically, each demo is simultaneously a guided description of the entire example and a BASH script, which can
be directly executed to run the complete demo after setting the few required paths at the beginning of the file. As
such, the user can both read the guided description as a web page and run the same script to execute the demo. At this
time the following demos are available (more demos will be added in the near future):

| Demo | Description |
|------|-------------|
| [demo-complete](examples/demo-complete/demo-complete.html) | A demo of the complete HATCHet pipeline starting from an example dataset of tumour and matched normal BAM files |
| [demo-wgs-sim](examples/demo-WGS-sim/demo-wgs-sim.html) | A demo on a typical WGS (whole-genome sequencing) multi-sample dataset with standard noise and variance of the data |
| [demo-wgs-cancer](examples/demo-WGS-cancer/demo-wgs-cancer.html) | A demo on a cancer WGS (whole-genome sequencing) multi-sample dataset with high noise and variance of the data |
| [demo-wes](examples/demo-WES/demo-wes.html) | A demo on a cancer WES (whole-exome sequencing) multi-sample dataset, which is typycally characterized by very high variance of RDR |

## Custom pipelines
<a name="custompipelines"></a>

The repository includes custom pipelines which have been designed to adapt the complete pipeline of HATCHet to special
conditions or to integrate the processed data produced by other pipelines. Each custom pipeline is a variation of the
main HATCHet's pipeline, we thus recommend the user to always first carefully understand the main
[BASH script](script/README.md) through the corresponding guided [tutorial](doc_runhatchet.md) and to carefully
understand the provided [demos](#demos) to properly apply HATCHet for best-quality results. Each custom pipeline also
includes a specific demo which represent a guided and executable example on example data.

| Name | Description | Script | Demo | Variations |
|------|-------------|--------|------|------------|
| GATK4-CNV | Custom pipeline for segmented files from GATK4 CNV pipeline | [custom-gatk4-cnv.sh](custom-gatk4-cnv.sh) | [demo-gatk4-cnv.sh](custom/GATK4-CNV/demo-gatk4-cnv.html) | This custom pipeline takes the input the segmented files which already contain the estimated RDR and BAF. As such, the first pre-processing steps of HATCHet (`count-reads`, `count-alleles`, and `combine-counts`) are not needed; for this reason, the following depdencies SAMtools and BCFtools and the following required data, human reference genome, matched-normal sample, and BAM files, are not needed in this case. |

## Detailed steps
<a name="detailedsteps"></a>

The full pipeline of HATCHet is composed of 9 sequential steps (and an additional dependency-checking command), starting from the required input data.
The description of each step also includes the details of the corresponding input/output that are especially useful when
one wants to replace or change some of the steps in the pipeline while guaranteeing the correct functioning of HATCHet.
Each step `<step>` of HATCHet can be run with the following command within a HATCHet conda environment:
```shell
hatchet <step>
```

**Note**: This version of HATCHet uses variable-width bins to ensure that each bin has comparable B-allele frequency (BAF) signal from heterogeneous germline SNPs and to account for sequencing coverage. To run the older versions of HATCHet with fixed-width bins, use [*count-reads-fw*](doc_count_reads_fw.html) (formerly binBAM) instead of *count-reads* and [*combine-counts-fw*](doc_combine_counts_fw.html) (formerly comBBo) instead of *combine-counts*. If you are using the `run` command, set `fixed_width = True` in your `.ini` file.

*Older versions of HATCHet used different names for these steps. The `Old Name` column lists those names.*

| Order | Step | Old Name | Description |
|-------|------|----------|-------------|
| (1)   | [*genotype-snps*](doc_genotype_snps.html)   | deBAF    | This step calls heterozygous germline SNPs from the matched-normal sample. |
| (2)   | [*count-alleles*](doc_count_alleles.html)   | deBAF    | This step counts the number of reads covering both the alleles of each identified heterozgyous SNP in every tumor sample. |
| (3)   | [*count-reads*](doc_count_reads.html)       | N/A   | This step splits the human reference genome into very small regions (representing candidate bins) according to germline SNP positions, and counts the number of total reads aligned to each region in each of the tumor samples and in the matched normal sample. |
| (4)   | [*phase-snps*](doc_phase_snps.html)       | N/A   | This **optional** step uses reference-based phasing to group together germline SNPs that are likely on the same haplotype. This improves the quality of BAF estimates for balanced regions of the genome. |
| (5)   | [*combine-counts*](doc_combine_counts.html) | N/A   | This step constructs genomic bins with variable sizes, and combines the read counts and the allele counts for the identified germline SNPs to compute the read-depth ratio (RDR) and B-allele frequency (BAF) of every genomic bin. |
| (6)   | [*cluster-bins-loc*](doc_cluster_bins_loc.html)*  | N/A    | This step globally clusters genomic bins along the entire genome and jointly across tumor samples while incorporating local information. |
| (7)   | [*plot-bins*](doc_plot_bins.html)           | BBot     | This **optional**  step produces informative plots concerning the computed RDRs, BAFs, and clusters. The information produced by this step are important to validate the compute clusters of genomic regions. When this step is enabled in the `.ini` file, [plot-bins-1d2d](doc_plot_bins_1d2d.html) is run as well to generate additional plots. |
| (8)   | [*compute-cn*](doc_compute_cn.html)         | hatchet  | This step computes allele-specific fractional copy numbers, solves a constrained distance-based simultaneous factorization to compute allele and clone-specific copy numbers and clone proportions, and deploys a model-selection criterion select the number of clone by explicitly considering the trade-off between subclonal copy-number aberrations and whole-genome duplication. |
| (9)   | [*plot-cn*](doc_plot_cn.html)               | BBeval   | This **optional** step analyzes the inferred copy-number states and clone proportions and produces informative plots jointly considering all samples from the same patient. In addition, this step can also combine results obtained for different patients and perform integrative analysis. When this step is enabled in the `.ini` file, [plot-cn-1d2d](doc_plot_cn_1d2d.html) is run as well to generate additional plots. |
| (10)   | [*check*](doc_check.html)     | check-solver           | This **optional** step is a quick way to verify if your solver and other 3rd party dependencies are working correctly. It checks for the presence of several dependencies and tests `compute_cn` on a set of small data files pre-packaged.  |

*This is a new clusting module that replaces the previous GMM-based clustering (command named cluBB). Users can still use GMM-based clustering if they wish by setting `loc-clust` to `False` under the `[run]` header in `hatchet.ini`, or by running [*cluster-bins*](doc_cluster_bins.html) directly.

## Recommendations and quality control
<a name="recommendations"></a>

All the components of HATCHet's pipeline use some basic parameters that allow to deal with data characterized by different features. The default values of these parameters allow one to succesfully apply HATCHet on most datasets. However, special or noisy datasets may require to tune some parameters. The user can deal with these cases by following the recommendations reported here, reading the descriptions of the various steps, and using the informative plots to verify the results. In the following guides and recommentations, we guide the user in the interpration of HATCHet's inference, we explain how to perform quality control to guarantee the best-quality results, and we describe how the user can control and tune some of the parameters to obtain the best-fitting results. We thus split the recommendations into distinct topics with dedicated descriptions.

| Recommendation | Description |
|----------------|-------------|
| [Analyze HATCHet inference](recommendation_inference.html) | Interpret HATCHet's inference, quality and error control, and investigate alternative solutions. |
| [Analyze global clustering](recommendation_clustering.html) | Interprent global clustering, quality and error control, and parameter tuning |
| [Analyze different type of data](recommendation_datatype.html) | Tuning parameters to better analyzing different type of data as those from WES |
| [Improve running time](recommendation_runtime.html)| Tips for improving running time of the whole pipeline |
