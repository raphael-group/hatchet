# `plot-panel`
`plot-panel` performs multi-sample copy-number panel plotting, composed from several per-sample CN solutions.

## Input

A panel TSV (`--panel_file`) listing the per-sample BBC UCN solution paths, plus the reference `--genome_size` / `--region_bed`. See [reference.md#output](../reference.md#output) for the `results/` UCN file layout.

## Usage

```console
$ hatchet plot-panel --help
usage: hatchet plot-panel [-h] --panel_file PANEL_FILE
                          --genome_size GENOME_SIZE --region_bed REGION_BED
                          [--width WIDTH] [--height HEIGHT]
                          [--show_clone_name | --no-show_clone_name]
                          [--show_prop | --no-show_prop]
                          [--show_ploidy | --no-show_ploidy]
                          [--min_prop MIN_PROP] [--dpi DPI] [--transparent]
                          [--title TITLE] -o OUT_FILE [--plot_1d2d]
                          [--plot_summary]
```

## Main parameters

### Panel input

`--panel_file` is a TSV listing the per-sample BBC UCN solution paths to compose into one multi-sample panel; `-o`/`--out_file` sets the output figure (e.g. `panel.svg`).

### Layout and annotations

- **Layout (`--width`, `--height`, `--dpi`, `--transparent`, `--title`).** Overall panel width and per-row height in inches (defaults 20 and 1), output resolution, background, and title.

- **Annotations (`--show_clone_name`, `--show_prop`, `--show_ploidy`, `--min_prop`).** Toggle clone names, proportions, and per-clone ploidy on each CN profile; `--min_prop` (default 0.01) hides clones below that proportion.

### Extras

`--plot_1d2d` also runs `plot-cn` per panel row (needs a `PATH_TO_BBC` column in the panel file); `--plot_summary` emits per-sample purity and ploidy barplots.

See [reference.md#plot-panel](../reference.md#plot-panel) for the full parameter table.

## Output

The composed multi-sample panel figure written to `--out_file`. With `--plot_1d2d` it also emits per-row `plot-cn` figures, and with `--plot_summary` per-sample purity and ploidy barplots.
