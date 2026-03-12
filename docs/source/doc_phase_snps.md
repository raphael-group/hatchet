# phase-snps

**Note:** To run this step, you must first run [*download-panel*](doc_download_panel.html) to download the reference-based phasing panel, and specify its location via the argument `-D, --refpaneldir`. The `download-panel` command only needs to be run once per system.

This step of HATCHet phases genotypes found in VCF files. It automatically takes care of differences in coordinates if the user has aligned their reads to a version of the reference genome (e.g. hg38) that is different from the version used in the reference panel (e.g. hg19 for the 1000 genomes project), using the liftover utility within Picard. Once genotypes are lifted over and phased, we again perform one last liftover to the coordinates of the genome used for alignment. Liftover is skipped if the version of the reference genome used for alignments corresponds to the same version used in the reference panel. Lastly, in order to account for differences in naming conventions for chromosomes, with or without the "chr" prefix, we also add or remove these prefixes so that chromosome names correspond to those used in the reference panel (without "chr" for the 1000 genomes project hg19 panel).

## Input

`phase-snps` takes one or more VCF files containing heterozygous SNP positions as input, specified using `-L, --snps`. These are typically produced by [`genotype_snps`](doc_genotype_snps.md).

The following parameters are required to specify and describe the input data:

| Name | Description | Usage |
|------|-------------|-------|
| `-L`, `--snps` | A list of VCF files to phase, one per chromosome | Specify a list using a path along with a wildcard, e.g. /path/to/snps/*.vcf.gz |
| `-D`, `--refpaneldir` | Path to the Reference Panel | This is the location where the 1000 genome project reference panel will be downloaded |
| `-g`, `--refgenome` | Path to the reference genome used to align reads | Path should include the filename |
| `-V`, `--refversion` | Version of reference genome used to align reads | Specify the human reference genome used for aligning reads; hg19 or hg38 |
| `-N`, `--chrnotation` | Chromosome names contain "chr" (e.g., "chr1" instead of "1") | Set this flag if and only if your BAM files/reference genome include "chr" in chromosome names |
| `-o`, `--outdir` | Output folder for phased VCFs | Specify a path or relative path |

## Output

The following files will be placed in the directory indicated by `-o, --outdir`:

| Name | Description |
|------|-------------|
| phased.vcf.gz | VCF file containing phased genotypes for all chromosomes |
| phased.log | Table showing how many SNPs were present before and after phasing; SNPs may be lost in the process due to various reasons, e.g. they were not present in the reference panel so there is no information to phase them |
| *_alignments.log | Log files from the shapeit phasing program, one per chromosome |

## Main Parameters

If HATCHet is installed via `conda`, the dependencies (`shapeit`, `picard`, `bcftools`, and `bgzip`) should be installed automatically and available on the PATH.
You can verify that these dependencies are available by running [the check command](doc_check.md), i.e., `hatchet check`.

 If HATCHet is installed from source, you may need to install them yourself (i.e., via `conda` or from source) and/or specify their locations using the following arguments:

| Name | Description | Usage |
|------|-------------|-------|
| `-j`, `--processes` | Number of parallel jobs | Parallel jobs independently phase VCF files, which are split up by chromosome |
| `-si`, `--shapeit` | Path to the `shapeit` executable | `shapeit` is required to run this command |
| `-pc`, `--picard` | Path to the `picard` executable or JAR file | `picard` is required to run this command |
| `-bt`, `--bcftools` | Path to the `bcftools` executable | `bcftools` is required to run this command |
| `-bg`, `--bgzip` | Path to the `bgzip` executable | `bgzip` is required to run this command |
