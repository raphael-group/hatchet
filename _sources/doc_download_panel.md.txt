# download-panel

This step of HATCHet downloads the 1000 genomes reference panel to phase germline mutations. It also downloads and creates other files necessary for phasing when the user has aligned their reads to a version of the human reference genome that is not the same as the version used in the 1000 genomes project (which was hg19, with no 'chr' prefix preceding chromosome names).

**Note:** This step requires access to the internet in order to download files. This step only needs to be run once per system.

## Input

`download-panel` takes no input files and the following options:

| Name | Description | Usage | Default |
|------|-------------|-------|---------|
| `-D`, `--refpaneldir` | Path to the Reference Panel | This is the location where the 1000 genome project reference panel will be downloaded | None
| `-R`, `--refpanel` | Reference Panel to use; only 1000 genomes project currently supported | specify "1000GP_Phase3" to automatically download and use the panel from the 1000 genomes project | "1000GP_Phase3"

## Output

download-panel populates the directory `-D, --refpaneldir` with several files for phasing germline SNPs, and for performing liftover in situations when the reference genome used to align sequencing reads is different from the reference genome version used in the reference panel.

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
