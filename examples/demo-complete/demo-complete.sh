# Demo complete for the entire HATCHet pipeline
: ex: set ft=markdown ;:<<'```shell' #

The following HATCHet demo represents a guided example of the complete HATCHet pipeline starting from an exemplary dataset of tumour and matched normal 
[BAM files](https://doi.org/10.5281/zenodo.4046906) publicly available. From this directory, simply run this file through BASH as a standard script to run 
the complete demo. The demo can also be considered as a guided example of a complete execution and is correspondingly commented.

## Requirements and set up

The demo requires that HATCHet has been succesfully installed in the current python environment.
Please make sure that you can succesfully run the required `samtools` and `bcftools`.
The demo includes the downloading of all the required files and will terminate in <20 minutes on machine with minimum requirements satisfied.

We gurantee that the running directory in the same directory of the demo and we remove previous results.

```shell
cd $( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )
rm -rf rdr/ baf/ snps/ bb/ bbc/ analysis/ results/ evaluation/
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

The demo auomatically downloads the required tumor and matched-normal BAM files in `data` folder.

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

1. We copy over the default config.sh and make custom changes to it at the end
```shell
cp ../../script/config.sh config.sh
:<<'```shell' # Ignore this line
```
 
2. We specify the correct path to the reference genome
```shell
echo 'REF="data/hg19.fa"' >> config.sh
:<<'```shell' # Ignore this line
```

3. We specify the current folder as the running one
```shell
echo 'XDIR="./"' >> config.sh
:<<'```shell' # Ignore this line
```

4. We specify the path to the matched-normal BAM files
```shell
echo 'NORMAL="data/normal.bam"' >> config.sh
:<<'```shell' # Ignore this line
```

5. We specify the list of paths to the tumor BAM files and corresponding names
```shell
echo 'BAMS="data/bulk_03clone1_06clone0_01normal.sorted.bam data/bulk_08clone1_Noneclone0_02normal.sorted.bam data/bulk_Noneclone1_09clone0_01normal.sorted.bam"' >> config.sh
echo 'NAMES="TumorSample1 TumorSample2 TumorSample3"' >> config.sh
:<<'```shell' # Ignore this line
```

6. We keep the default number of reads and number of parallel processes
```shell
echo 'J=$(python -c "import multiprocessing as mp; print(mp.cpu_count())")' >> config.sh
echo "MINREADS=8" >> config.sh
echo "MAXREADS=300" >> config.sh
:<<'```shell' # Ignore this line
```

7. We specify the reference genome and chr notation
```shell
echo 'REF_VERS="hg19"' >> config.sh
echo 'CHR_NOTATION=true' >> config.sh
:<<'```shell' # Ignore this line
```

8. We add the unphased version of HATCHet's script
```shell
cp ../../script/runUnphased.sh runUnphased.sh
:<<'```shell' # Ignore this line
```

## Running HATCHet

```shell
bash runUnphased.sh |& tee hatchet.log
exit $?
```

