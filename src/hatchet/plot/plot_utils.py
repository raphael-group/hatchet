"""Shared helpers for the plot-* command modules."""

import matplotlib.pyplot as plt

from cnplot import GenomeAxis, read_chr_sizes


def get_plot_style(args: dict) -> dict:
    """Bundle the ``plot_*`` styling knobs, with the ``plot_`` prefix stripped.

    ``normalize_args`` already merges the hatchet.yaml ``plot_*`` defaults into
    ``args``, so this is the single shared accessor for the plot commands: each
    reads the keys it needs (e.g. ``style["row_width"]``, ``style["diag_hist_bins"]``).
    Defaults therefore live only in hatchet.yaml. Non-styling ``plot_*`` keys
    (``plot_dir``, ``plot_1d2d``, ...) come along harmlessly and are simply unread.
    """
    return {k[len("plot_") :]: v for k, v in args.items() if k.startswith("plot_")}


def use_editable_fonts():
    """Embed editable fonts in vector output (keeps text selectable in pdf/ps/svg)."""
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42
    plt.rcParams["svg.fonttype"] = "none"


def build_genome_axis(region_bed, genome_size, keep_chroms=None, collapse_gaps=True):
    """Build a cnplot GenomeAxis, optionally limited to chromosomes in the data.

    Args:
        region_bed: Region BED file path.
        genome_size: Chromosome-sizes file path.
        keep_chroms: If given, drop every chromosome not in this iterable (e.g.
            sex chromosomes absent from an autosome-only run). None keeps all.
        collapse_gaps: Shrink uncovered stretches from the axis.

    Returns:
        A cnplot GenomeAxis.
    """
    excluded = ()
    if keep_chroms is not None:
        keep = set(keep_chroms)
        excluded = [c for c in read_chr_sizes(genome_size) if c not in keep]
    return GenomeAxis(
        region_bed, genome_size, excluded_chroms=excluded, collapse_gaps=collapse_gaps
    )
