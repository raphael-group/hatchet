# Tutorial

This tutorial illustrates how to use the complete pipeline through the `hatchet run` command, along with the configuration [hatchet.ini](hatchet.ini) file that stores all the variables the user needs to specify.

The tutorial is subdivided into some subsections and each of these describes sequential parts of the full pipeline:
1. [Preliminaries](#preliminaries)
2. [genotype-snps](#genotype-snps)
3. [count-alleles](#count-alleles)
4. [count-reads](#count-reads)
5. [combine-counts](#combine-counts)
6. [cluster-bins](#cluster-bins)
7. [plot-bins](#plot-bins)
8. [compute-cn](#compute-cn)
9. [plot-cn](#plot-cn)

We suggest you make a copy of the configuration file, place it into the designated running directory, and follow the tutorial.

## Preliminaries
<a name="preliminaries"></a>

The following variables, specified in [hatchet.ini](hatchet.ini), should be changed according to the user's data:

```shell
# Path to reference genome - make sure you have also generated the reference dictionary as /path/to/reference.dict
reference = "/path/to/reference.fa"
normal = "/path/to/normal.bam"
# Space-delimited list of tumor BAM locations
bams = "/path/to/tumor1.bam /path/to/tumor2.bam"
# Space-delimited list of tumor names
samples = "Primary Met"

```

Each of these are explained further below.

***

```shell
reference = "/path/to/reference.fa"
reference_version = ""
chr_notation = True
```

First, you need to specify the full path to the human reference genome with the variable `reference`, along with the specific version of this reference with the `reference_version` variable, which may be "hg19" or "hg38". This reference version is used to select a list of known germline SNPs to genotype samples. Lastly, the `chr_notation` must be set to true/false depending on whether or not the chromosomes in your reference are prefixed with "chr".

***

```shell
normal = "/path/to/normal.bam"
bams = "/path/to/tumor1.bam /path/to/tumor2.bam"
samples = "Primary Met"
```

Next, you need to specify the full paths to the required input data:
1. `normal` is the full path to the BAM file of matched-normal samples
2. `bams` is a white-space separated list of the BAM files for the multiple tumor samples from the considered patient.

The variable `samples` is also a white-space separated list of tumor sample names (specified in the same order as the BAM files in `bams`), and these names are used in the plots produced by HATCHet.

***

```shell
mincov = 8
maxcov = 300
```

Both `mincov` and `maxcov` specify the minimum and maximum sequencing depth for considering germline snps. For samples that have undergone whole-genome sequencing with depths >30X, values around 8 and 300 should be reasonable for the minimum and maximum, respectively. For samples in which only the exome was sequenced to ~100X, then values around 20 and 1000 should be fine for the minimum and maximum, respectively.

***

## genotype-snps
<a name="genotype-snps"></a>

```shell
python3 -m hatchet genotype-snps -N ${NORMAL} -r ${REF} -j ${J} -c ${MINREADS} -C ${MAXREADS} \
                            -R ${LIST} -o ${SNP} |& tee ${BAF}bafs.log
```

genotype-snps genotypes germline SNPs from the matched-normal sample `${NORMAL}`, using positions of known germline variation specified by the file in `${LIST}`. If this parameter is omitted, heterozygous SNPs are inferred from the normal sample. As mentioned above, SNPs are only considered for downstream analysis if they have a minimum and maximum sequencing depth of `${MINREADS}` and `${MAXREADs}`, respectively.

GATK best practices suggest that the maximum should be at least twice the expected average coverage to avoid mapping artifacts.
Observe that WES generally requires higher thresholds (e.g. 100 and 3000).
Also, an increasing in the minimum threshold allow to improve the quality of the estimated BAF (computed in the next step), while a decreasing  it allows us to consider more SNPs.

## [count-alleles](doc_count_alleles.md)
<a name="count-alleles"></a>

```shell
python3 -m hatchet count-alleles -N ${NORMAL} -T ${BAMS} -S ${ALLNAMES} -r ${REF} \
                          -j ${J} -c ${MINREADS} -C ${MAXREADS} -L ${SNP}*.vcf.gz \
                          -O ${BAF}normal.1bed -o ${BAF}tumor.1bed |& tee ${BAF}bafs.log
```

count-alleles uses these germline SNPs from genotype-snps (specified in the list `${SNP}*.vcf.gz`) and computes the corresponding allele counts from all tumor samples `${BAMS}`.
As for all steps in this pipeline, the same reference genome used for the BAM files has to be specified by `${REF}`.
Minimum and maximum thresholds for the read counts of germline SNPs are again specified with `${MINREADS}` and `${MAXREADS}`.

Several simple parameters are also specified including: number of parallel threads, output filenames, and the log filename of this step.


## [count-reads](doc_count_reads.md)
<a name="count-reads"></a>

```shell
python3 -m hatchet count-reads -N ${NORMAL} -T ${BAMS} -S ${ALLNAMES} -b ${BAF}tumor.1bed \
                              -j ${J} -O ${RDR}normal.1bed -O {RDR} -V {REF_VERSION} |& tee ${RDR}bins.log
```

`count-reads` identifies candidate bin thresholds using the SNPs identified by `count_alleles` (in `${BAF}tumor.1bed `), and counts the number of reads in each tumor sample `${BAMS}` as well as in the normal sample `{$NORMAL}` that fall within each interval as well as the number of reads overlapping each threshold. It uses the centromere locations corresponding the the reference genome version `${REF_VERSION}` (hg19 or hg38) to avoid constructing bins that overlap centromeres. This step produces a gzip-compressed file for each chromosome and each sample in the `${RDR}` directory, as well as an ordered list `${RDR}samples.txt` of sample labels. It also counts the total reads in each of the BAM files and writes it to `${RDR}total.tsv`.
Other simple parameters are also specified including number of parallel threads, output filenames, and the log filename of this step.


## [combine-counts](doc_combine_counts.md)
<a name="combine-counts"></a>

```shell
python3 -m hatchet combine-counts -A ${RDR} -t ${RDR}total.tsv -b {BAF}tumor.1bed -o #{BB}bulk.bb \
                          --msr ${msr} --mtr ${mtr} -j ${J} \
                          -p ${PHASE} -s ${max_blocksize} -m {max_spb} -a {alpha} # optional phasing args
```

combine-counts constructs genomic bins such that in all samples, each bin has at least `${msr}` SNP-covering reads and at least `${mtr}` total reads. Bins will not have the same width, but using this rule each bin will have comparable RDR and BAF signals for the following clustering steps. The BAF for each bin and the relative phase of all SNPs in the bin are inferred via EM, and the RDR for each bin is computed fom the read count files in the `${RDR}` folder. After this computation, RDR values are normalized using the total reads in each sample (from `${RDR}total.tsv`). As with other HATCHet commands, `-j ${J}` controls the number of parallel processes.

See the [script](script/README.md) directory for a guide on how to run HATCHet with phasing. If a phased VCF file is supplied via `-p, --phase ${PHASE}` (e.g., `-p phase/phased.vcf.gz`), SNPs are merged into blocks before BAF inference. Each block contains at most `${max_spb}` such that no two SNPs in the same block are further apart than `${max_blocksize}`, and such that no two adjacent SNPs have significantly different marginal BAF estimates (at significance level `${alpha}` -- higher `${alpha}` corresponds to less trust in the phasing results). Then, blocks are passed to the EM which determines the relative phase of each block.

## [cluster-bins](doc_cluster_bins.md)
<a name="cluster-bins"></a>

```shell
python3 -m hatchet cluster-bins ${BB}bulk.bb -o ${BBC}bulk.seg -O ${BBC}bulk.bbc \
                                      -minK 5 -maxK 10 -tau 0.0001 -d 0.08
```

cluster-bins globally clusters genomic bins based on RDR and BAF (specified in a BB file `${BB}bulk.bb`) jointly along the genome and across all tumor samples. This *locality-aware* clustering takes into account the contiguity of bins along the genome to inform the clustering.

The main parameters are as follows:
- The maximum expected BAF shift `-d` for a balanced cluster, here set to `0.08`. This parameter should be as small as possible while still including the clusters with the closest BAF to 0.5.
- Parameters `minK` and `maxK` which specify the minimum and maximum (respectively) of clusters to consider. Alternatively, one can specify `exactK` to fix the number of clusters at a given value.
- Parameter `tau` which controls the trade-off between local and global information: lower `tau` puts more weight on *local* information (i.e., adjacent bins on each chromosome arm) while higher `tau` puts more weight on *global* information (i.e., RDR and BAF values).

The plots from `plot-bins` are helpful to evaluate clusterings and tune parameters.

Other parameters are used to specify the output files containing the cluster information (`${BBC}bulk.seg`) and the clustered genomic bins (`${BBC}bulk.bbc`).

## [plot-bins](doc_plot_bins.md)
<a name="plot-bins"></a>

```shell
cd ${PLO}
python3 -m hatchet plot-bins -c RD --figsize 6,3 ../${BBC}bulk.bbc
python3 -m hatchet plot-bins -c CRD --figsize 6,3 ../${BBC}bulk.bbc
python3 -m hatchet plot-bins -c BAF --figsize 6,3 ../${BBC}bulk.bbc
python3 -m hatchet plot-bins -c BB ../${BBC}bulk.bbc
python3 -m hatchet plot-bins -c CBB ../${BBC}bulk.bbc -tS 0.01
```

plot-bins produces informative plots which are described [here](doc_plot_bins.md).
Many of these plots can be very useful to assess the performance of the various steps of HATCHet, especially in the case of noisy datasets.

## [compute-cn](doc_compute_cn.md)
<a name="compute-cn"></a>

```shell
cd ../${RES}
python3 -m hatchet compute-cn -i ../${BBC}bulk -n2,6 -p 400 -u 0.03 \
                          -r ${RANDOM} -j ${J} -eD 6 -eT 12 -g 0.35 \
                          -l 0.6 &> >(tee >(grep -v Progress > hatchet.log))
```

compute-cn computes the fractional copy numbers, factorizes these into allele and clone-specific copy numbers and clone proportions, and applies a model-selection step to jointly infer the number of clones and predict a whole-genome duplication (WGD).
This step requires the common prefix `${BBC}bulk` of the cluster and clustered bins files.
Some basic parameter include:
- the interval `2,6`  specifying the minimum and maximum number of clones which can be increased as much as wanted. User can start to use this interval and increasing if the number of clones appear to be higher or close to the maximum.
- the number of restarts `400` for the coordinate-descent method. This value is typically enough to obtain a good solution close to the optimum, however an increase in this value allows to improve the search of a solution. This can be useful when observing a number of clusters significantly higher than usual.
- Random initializing seed with `${RANDOM}`.
- Number of parallel threads `${J}`

There are also some important parameters which need to be tuned for special or noisy datasets.
The details of these parameters are the following:
- Maximum copy number without a WGD `-eD` and with a WGD `-eT`. The specified values are typically enough to obtain the copy numbers of large CNAs. While higher values increase the number of reasonable solutions, the values can be increased to investigate the presence of large genomic regions with higher copy numbers. In particular, when the value `0` is specified, then the maximum copy numbers are automatically chosen from the fractional copynumbers.
- Minimum clone proportion `-u`. The value depends on the power to identify tumor clones present in samples at low tumor proportions. The default value is reasonable for most of the datasets, but datasets with high noise may need higher values to avoid overfitting. We thus suggest, especially in these special cases, to repeat the hatchet step by considering increasing values of this parameter up to a maximum or until clone proportions are different from this threshold.
- Two parameters are used to control the elbow method during model selection to choose the best solution. The user should repeat this step with different values of these parameters, especially when one wants to investigate the presence of a single tumor clone or the presence of more tumor clones with small different CNAs and small clone proportions. These two parameters are the following:
  - Confidence in the presence of a single tumor clone `g`. This value determines the confidence in having a single tumor clone in several of the given samples. The user can investigate the hypothesis that a single tumor clone is present by increasing the confidence, e.g. using values of `0.4, 0.5, ...`. Vice versa, the confidence is lowered by decreasing these values. For example, one can do that to investigate the presence of tumor clones with small CNAs or small proportions.
  - Sensitivity `l`. This value determines the sensitivity of the method. Lower values correspond to higher sensitivities which allow to consider solutions with a higher number of clones. The user can decrease this value for investigating the presence of more tumor clones, especially tumor clones with small different CNAs or present in small proportions.

## [plot-cn](doc_plot_cn.md)
<a name="plot-cn"></a>

```shell
cd ../${SUM}
python3 -m hatchet plot-cn ../${RES}/best.bbc.ucn
```

plot-cn produces informative and summary plots for the inferred results.
Examples and details of these plots are reported [here](doc_plot_cn.md).
