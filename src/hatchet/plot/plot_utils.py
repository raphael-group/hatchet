"""Shared helpers for the plot-* command modules."""

import matplotlib.pyplot as plt

from cnplot import GenomeAxis, read_chr_sizes


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
