# download-panel

This step of HATCHet downloads the 1000 genomes reference panel to phase germline mutations. It also downloads or creates other files necessary for phasing when the user has aligned their reads to a version of the human reference genome that is not the same as the version used in the 1000 genomes project (hg19, with no 'chr' prefix preceding chromosome names). 

This step only needs to be run once for each combination of "refversion" and "chrnotation" you wish to use -- for example, you must run `download_panel -D <my_panel> -V hg19 -N false` before running on samples aligned to hg19 without "chr" notation, but then download_panel does not need to be run until you would like to phase a sample aligned to hg38 and/or with "chr" notation (at which point you would run e.g. `download_panel -D <my_panel> -V hg38 -N true` to add the extra files needed to the directory `my_panel`).

NOTE: This step requires access to the internet in order to download files.

## Input

download-panel takes as input variables specified within the hatchet.ini file.

| Name | Description | Usage | Default |
|------|-------------|-------|---------|
| `-D`, `--refpaneldir` | Path to the Reference Panel | This is the location where the 1000 genome project reference panel will be downloaded | "" (current location)
| `-R`, `--refpanel` | Reference Panel to use; only 1000 genomes project currently supported | specify "1000GP_Phase3" to automatically download and use the panel from the 1000 genomes project | "1000GP_Phase3"
| `-V`, `--refversion` | Version of reference genome used in BAM files | Specify the human reference genome used for aligning sequencing reads; hg19 or hg38 | (none, argument required)
| `-N`, `--chrnotation` | Notation of chromosomes, with or without "chr" | Specify "true" or "false" according to whether chromosomes are prefixed by "chr" | (none, argument required)
| `--st`, `--samtools` | Path to samtools executable |  | "" (look on PATH)

## Output

download-panel produces various files for phasing germline SNPs, and performing liftover in situations when the reference genome used to align sequencing reads is different from the reference genome version used in the reference panel.

| Name | Description |
|------|-------------|
| 1000GP_Phase3 | Directory containing reference panel of haplotypes for phasing mutations, in HAP/LEGEND/SAMPLE format |
| hg19_no_chr.fa | hg19 reference genome, chromosomes named without "chr" prefix |
| hg19_no_chr.dict | Sequence dictionary for the hg19 reference genome |
| hg38ToHg19.chr.chain | Chain file for liftover, going from hg38 to hg19 coordinates |
| hg38ToHg19.no_chr.chain | Chain file for liftover, going from hg38 to hg19 coordinates, where hg38 chromosomes are named without "chr" prefix |
| hg19ToHg38.chr.chain | Chain file for liftover, going from hg19 to hg38 |
| hg19ToHg38.no_chr.chain | Chain file for liftover, going form hg19 to hg38, where hg38 chromosomes are named without "chr" prefix |
| rename_chrs1.txt | File for renaming chromosomes with bcftools, removing the "chr" prefix |
| rename_chrs2.txt | File for renaming chromosomes with bcftools, adding the "chr" prefix |

