# `plot-cn`
`plot-cn` performs genome-wide copy-number profile and RDR-vs-BAF plotting for a single CN solution.

## Input

A per-bin BBC UCN table (`--bbc`) and per-cluster SEG UCN table (`--seg`) from `compute-cn`, the `--gamma_file` (`gammas.tsv`), and the reference `--genome_size` / `--region_bed`; `--ploidy` selects the gamma column. See [reference.md#output](../reference.md#output) for the `results/` UCN file layout.

## Usage

```console
$ hatchet plot-cn --help
usage: hatchet plot-cn [-h] --bbc BBC --seg SEG --genome_size GENOME_SIZE
                       --region_bed REGION_BED -g GAMMA_FILE [-s SOLFILE]
                       -O PLOT_DIR [--dpi DPI] [--img_type {pdf,png,svg}]
                       [--transparent] [--show_gap | --no-show_gap]
                       [--tail_alpha TAIL_ALPHA] [--center_alpha CENTER_ALPHA]
                       [--onetail_area ONETAIL_AREA] [--maxlim_fcn MAXLIM_FCN]
                       --ploidy {diploid,tetraploid} [--patient_id PATIENT_ID]
```

## Main parameters

### Solution inputs

- **CN solution (`--bbc`, `--seg`, `-g`/`--gamma_file`, `--ploidy`).** plot-cn reads a per-bin BBC UCN table and per-cluster SEG UCN table from compute-cn together with the `gammas.tsv` scaling file; `--ploidy` (`diploid`/`tetraploid`) selects which gamma column to apply.

- **Solution override (`-s`/`--solfile`).** Optional file that overrides the CN states in the BBC table, e.g. to plot a manually edited solution.

### Figure output

- **Destination and format (`-O`/`--plot_dir`, `--img_type`, `--dpi`, `--transparent`).** Figures are written to `--plot_dir` as `--img_type` (`png` default, or `pdf`/`svg`) at `--dpi` resolution, optionally on a transparent background.

- **Axis and layout (`--maxlim_fcn`, `--show_gap`).** `--maxlim_fcn` caps the FCN axis (default 30); `--show_gap` keeps gap regions instead of collapsing them.

### CN-state shading

`--tail_alpha`, `--center_alpha`, and `--onetail_area` control the transparency gradient that conveys per-CN-state confidence in the profile; a larger `--onetail_area` widens the low-confidence tails.

See [reference.md#plot-cn](../reference.md#plot-cn) for the full parameter table.

## Output

Copy-number figures written to `--plot_dir` in the `--img_type` format: genome-wide 1D profiles, the FCN A/B panel, and per-sample 2D RDR-vs-BAF scatter plots.
