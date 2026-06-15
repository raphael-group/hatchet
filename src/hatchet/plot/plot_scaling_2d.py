import os
import logging
import contextlib

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from adjustText import adjust_text

from hatchet.plot.plot_utils import set_palette

logging.getLogger("adjustText").setLevel(logging.WARNING)


def plot_scaling_2d(
    samples: list,
    bbcs: pd.DataFrame,
    segs: pd.DataFrame,
    scaling: dict,
    out_file: str,
    markersize: float = 3.0,
    markersize_centroid: float = 14,
    marker_bd_width: float = 0.8,
    dpi: int = 300,
    transparent: bool = False,
    maxlim_rdr: int = 10,
):
    """2D RDR-vs-BAF scatter anchoring the scaling inference from get_scaling_factor.

    One PDF page per sample with up to two panels: ax0 = noWGD (diploid),
    ax1 = WGD (tetraploid, drawn only when a WGD scaling was inferred). Bins
    are colored per cluster; anchor clusters are circled and annotated with
    their inferred (a, b) clonal states.
    """
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42
    plt.rcParams["svg.fonttype"] = "none"

    clusters = sorted(bbcs["CLUSTER"].unique().tolist())
    palette = set_palette(num_colors=len(clusters))

    panels = [("noWGD", scaling["diploid"])]
    if scaling.get("tetraploid") is not None:
        panels.append(("WGD", scaling["tetraploid"]))

    cent = segs.set_index(["#ID", "SAMPLE"])[["RD", "BAF"]]

    pdf = PdfPages(out_file)
    for sample in samples:
        sub = bbcs[bbcs["SAMPLE"] == sample]
        bafs, rdrs = sub["BAF"].to_numpy(), sub["RD"].to_numpy()
        hue = sub["CLUSTER"].to_numpy()

        lim_baf = (0, 1) if np.max(bafs) > 0.5 else (0, 0.55)
        lim_rdr = (0, min(max(3, int(np.ceil(np.max(rdrs)))), maxlim_rdr))

        fig, axes = plt.subplots(
            1, len(panels), figsize=(7 * len(panels), 6), squeeze=False
        )
        for ax, (label, info) in zip(axes[0], panels):
            sns.scatterplot(
                x=bafs,
                y=rdrs,
                hue=hue,
                hue_order=clusters,
                palette=palette,
                ax=ax,
                s=markersize,
                legend=False,
                edgecolor="none",
                linewidth=0,
                rasterized=True,
            )
            ax.axvline(0.5, color="grey", linestyle=":", linewidth=0.8)

            texts, ax_x, ax_y = [], [], []
            for c, (a, b) in info["clonal"].items():
                if (c, sample) not in cent.index:
                    continue
                cx, cy = cent.loc[(c, sample), "BAF"], cent.loc[(c, sample), "RD"]
                ax.scatter(
                    cx,
                    cy,
                    facecolors="none",
                    edgecolors="black",
                    s=markersize_centroid,
                    linewidth=marker_bd_width,
                )
                texts.append(
                    ax.text(cx, cy, f"({a},{b})", fontsize=10, fontweight="bold")
                )
                ax_x.append(cx)
                ax_y.append(cy)
            if texts:
                with open(os.devnull, "w") as _dn, contextlib.redirect_stdout(_dn):
                    adjust_text(
                        texts,
                        x=ax_x,
                        y=ax_y,
                        ax=ax,
                        arrowprops=dict(arrowstyle="-", color="black", lw=0.5),
                    )

            p = (info.get("purities") or {}).get(sample)
            ptxt = f"  purity={p:.3f}" if p is not None else ""
            ax.set(xlim=lim_baf, ylim=lim_rdr, xlabel="BAF", ylabel="RDR")
            ax.set_title(f"{label}{ptxt}")

        fig.suptitle(f"sample={sample}")
        fig.tight_layout()
        pdf.savefig(fig, dpi=dpi, bbox_inches="tight", transparent=transparent)
        plt.close(fig)
    pdf.close()
    logging.info(f"scaling 2D scatter saved to {out_file}")
