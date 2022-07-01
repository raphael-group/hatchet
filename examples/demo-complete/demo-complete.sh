# Demo complete for the entire HATCHet pipeline

: ex: set ft=markdown ;:<<'```shell' #

The following HATCHet demo represents a guided example of the complete HATCHet pipeline starting from an exemplary dataset of tumour and matched normal
[BAM files](https://doi.org/10.5281/zenodo.4046906) publicly available. From this directory, simply run this file through BASH as a standard script to run
the complete demo. The demo can also be considered as a guided example of a complete execution and is correspondingly commented.

## Requirements and set up

The demo requires that HATCHet has been succesfully installed in the current python environment.
Please make sure that you can succesfully run the required dependencies `samtools`, `bcftools`, `tabix`, and `mosdepth`.
The demo includes the downloading of all the required files and will terminate in <20 minutes on machine with minimum requirements satisfied.

We guarantee that the running directory in the same directory of the demo and we remove previous results.

```shell
cd $( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )
rm -rf rdr/ baf/ snps/ bb/ bbc/ plots/ results/ summary/
:<<'```shell' # Ignore this line
```

We also ask the demo to terminate in case of errors and to print a trace of the execution by the following commands
```shell
set -e
set -o xtrace
PS4='[\t]'
:<<'```shell' # Ignore this line
```

## Downloading of data

The demo automatically downloads the required tumor and matched-normal BAM files in `data` folder.

```shell
# Creating data folder
mkdir -p data

# Downloading matched-normal BAM file
echo "Downloading matched-normal BAM file from Zenodo, please be patient as downloading time may vary."
curl -L 'https://zenodo.org/record/4046906/files/normal.bam?download=1' > data/normal.bam
curl -L 'https://zenodo.org/record/4046906/files/normal.bam.bai?download=1' > data/normal.bam.bai

# Downloading BAM file of tumor sample 1
echo "Downloading BAM file of tumor sample 1 from Zenodo, please be patient as downloading time may vary."
curl -L 'https://zenodo.org/record/4046906/files/bulk_03clone1_06clone0_01normal.sorted.bam?download=1' > data/bulk_03clone1_06clone0_01normal.sorted.bam
curl -L 'https://zenodo.org/record/4046906/files/bulk_03clone1_06clone0_01normal.sorted.bam.bai?download=1' > data/bulk_03clone1_06clone0_01normal.sorted.bam.bai

# Downloading BAM file of tumor sample 2
echo "Downloading BAM file of tumor sample 2 from Zenodo, please be patient as downloading time may vary."
curl -L 'https://zenodo.org/record/4046906/files/bulk_08clone1_Noneclone0_02normal.sorted.bam?download=1' > data/bulk_08clone1_Noneclone0_02normal.sorted.bam
curl -L 'https://zenodo.org/record/4046906/files/bulk_08clone1_Noneclone0_02normal.sorted.bam.bai?download=1' > data/bulk_08clone1_Noneclone0_02normal.sorted.bam.bai

# Downloading BAM file of tumor sample 3
echo "Downloading BAM file of tumor sample 3 from Zenodo, please be patient as downloading time may vary."
curl -L 'https://zenodo.org/record/4046906/files/bulk_Noneclone1_09clone0_01normal.sorted.bam?download=1' > data/bulk_Noneclone1_09clone0_01normal.sorted.bam
curl -L 'https://zenodo.org/record/4046906/files/bulk_Noneclone1_09clone0_01normal.sorted.bam.bai?download=1' > data/bulk_Noneclone1_09clone0_01normal.sorted.bam.bai
:<<'```shell' # Ignore this line
```

Next the corresponding reference genome is downloaded and unpacked

```shell
echo "Downloading human reference genome, please be patient as downloading time may vary."
curl -L https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz | gzip -d > data/hg19.fa
samtools faidx data/hg19.fa
samtools dict data/hg19.fa > data/hg19.dict
:<<'```shell' # Ignore this line
```

## Configuring the HATCHet's execution

We follow the template of the HATCHet's [script](../../doc/doc_fullpipeline.md#fullpipelineandtutorial).

1. We specify the correct path to the reference genome and the output folder, and other required flags
```shell
echo '[run]' > hatchet.ini
echo 'genotype_snps=True' >> hatchet.ini
echo 'count_alleles=True' >> hatchet.ini
echo 'count_reads=True' >> hatchet.ini
echo 'combine_counts=True' >> hatchet.ini
echo 'cluster_bins=True' >> hatchet.ini
echo 'plot_bins=True' >> hatchet.ini
echo 'compute_cn=True' >> hatchet.ini
echo 'plot_cn=True' >> hatchet.ini
echo 'reference=data/hg19.fa' >> hatchet.ini
echo 'output=output/' >> hatchet.ini
j=$(grep -c ^processor /proc/cpuinfo)
processes="processes=${j}"
echo $processes >> hatchet.ini

:<<'```shell' # Ignore this line
```

2. We specify the path to the matched-normal BAM files
```shell
echo 'normal=data/normal.bam' >> hatchet.ini
:<<'```shell' # Ignore this line
```

3. We specify the list of paths to the tumor BAM files and corresponding names
```shell
echo 'bams=data/bulk_03clone1_06clone0_01normal.sorted.bam data/bulk_08clone1_Noneclone0_02normal.sorted.bam data/bulk_Noneclone1_09clone0_01normal.sorted.bam' >> hatchet.ini
echo 'samples=TumorSample1 TumorSample2 TumorSample3' >> hatchet.ini
:<<'```shell' # Ignore this line
```

4. We specify the min/max coverage for the genotpe_snps step
```shell
echo '[genotype_snps]' >> hatchet.ini
echo 'mincov=8' >> hatchet.ini
echo 'maxcov=300' >> hatchet.ini
:<<'```shell' # Ignore this line
```

5. We specify the reference genome and chr notation
```shell
echo 'reference_version=hg19' >> hatchet.ini
echo 'chr_notation=True' >> hatchet.ini
:<<'```shell' # Ignore this line
```

6. We specify mincov/maxcov for the count_alleles step
```shell
echo '[count_alleles]' >> hatchet.ini
echo 'mincov=8' >> hatchet.ini
echo 'maxcov=300' >> hatchet.ini
:<<'```shell' # Ignore this line
```

7. We specify the minimum number of total and SNP-covering reads in each bin for the combine_counts step
```shell
echo '[combine_counts]' >> hatchet.ini
echo 'msr=3000' >> hatchet.ini
echo 'mtr=5000' >> hatchet.ini
:<<'```shell' # Ignore this line
```

## Running HATCHet

```shell
python -m hatchet run hatchet.ini
exit $?
```
