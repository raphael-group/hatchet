# phase-snps

This step of HATCHet phases genotypes found in VCF files. It automatically takes care of differences in coordinates if the user has aligned their reads to a version of the reference genome (e.g. hg38) that is different from the version used in the reference panel (e.g. hg19 for the 1000 genomes project), using the liftover utility within Picard. Once genotypes are lifted over and phased, we again perform one last liftover to the coordinates of the genome used for alignment. Liftover is skipped if the version of the reference genome used for alignments corresponds to the same version used in the reference panel. Lastly, in order to account for differences in naming conventions for chromosomes, with or without the "chr" prefix, we also add or remove these prefixes so that chromosome names correspond to those used in the reference panel (without "chr" for the 1000 genomes project hg19 panel).

## Input

| Name | Description | Usage |
|------|-------------|-------|
| `-D`, `--refpaneldir` | Path to the Reference Panel | This is the location where the 1000 genome project reference panel will be downloaded |
| `-g`, `--refgenome` | Path to the reference genome used to align reads | Path should include the filename |
| `-V`, `--refversion` | Version of reference genome used to align reads | Specify the human reference genome used for aligning reads; hg19 or hg38 |
| `-N`, `--chrnotation` | Notation of chromosomes, with or without "chr" | Specify "true" or "false" according to whether chromosomes are prefixed by "chr" |
| `-o`, `--outdir` | Output folder for phased VCFs | Specify a path or relative path |
| `-L`, `--snps` | A list of VCF files to phase, one per chromosome | Specify a list using a path along with a wildcard, e.g. /path/to/snps/*.vcf.gz |
| `-j`, `--processes` | Number of parallel jobs | Parallel jobs independently phase VCF files, which are split up by chromosome |

## Output

| Name | Description |
|------|-------------|
| phased.vcf.gz | VCF file containing phased genotypes for all chromosomes |
| phased.log | Table showing how many SNPs were present before and after phasing; SNPs may be lost in the process due to various reasons, e.g. they were not present in the reference panel so there is no information to phase them |
| *_alignments.log | Log files from the shapeit phasing program, one per chromosome |
