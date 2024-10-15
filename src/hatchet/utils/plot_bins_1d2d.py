import os
import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib.colors as mcolors

from hatchet.utils.ArgParsing import parse_plot_bins_1d2d_args
from hatchet.utils import Supporting as sp

plt.rcParams["savefig.dpi"] = 300
plt.style.use("ggplot")
sns.set_style("white")


def main(args=None):
    sp.log("# Checking and parsing input arguments\n")
    args = parse_plot_bins_1d2d_args(args)
    sp.logArgs(args, 80)

    bbc = pd.read_table(args["bbc"])

    segf = args["seg"]
    if segf is not None:
        seg = pd.read_table(segf)

    show_centers = args["centers"]
    outdir = args["outdir"]
    alpha = args["alpha"]
    show_centromeres = args["centromeres"]

    if args["minbaf"] is not None:
        baf_lim = args["minbaf"], args["maxbaf"]
    else:
        baf_lim = None
    if args["minrdr"] is not None:
        rdr_lim = args["minrdr"], args["maxrdr"]
    else:
        rdr_lim = None

    colors1 = plt.cm.tab20b(np.arange(20))
    colors2 = plt.cm.tab20c(np.arange(20))

    # combine them and build a new colormap
    colors = np.vstack((colors1, colors2))
    global mymap
    mymap = mcolors.LinearSegmentedColormap.from_list("my_colormap", colors)

    sp.log(" Starting plotting\n", level="STEP")

    display = False  # setting display=True produces plots when run interactively
    plot_2d(
        bbc=bbc,
        seg=seg,
        show_centers=show_centers,
        xlim=baf_lim,
        ylim=rdr_lim,
        display=display,
        outdir=outdir,
        alpha=alpha,
    )

    plot_1d(
        bbc=bbc,
        baf_lim=baf_lim,
        rdr_lim=rdr_lim,
        display=display,
        outdir=outdir,
        alpha=alpha,
        show_centromeres=show_centromeres,
    )
    sp.log(" Done plotting\n", level="STEP")


def plot_2d(
    bbc,
    seg=None,
    show_centers=False,
    xlim=None,
    ylim=None,
    figsize=(4, 4),
    display=True,
    outdir=None,
    alpha=1,
):
    """
    For each sample, plot the mBAF and RDR of each bin colored by cluster.
    Colors will match clusters in corresponding 1D plots.
    """
    for s, df in bbc.groupby("SAMPLE"):
        plt.figure(figsize=figsize)
        plt.scatter(df.BAF, df.RD, c=df.CLUSTER, cmap=mymap, s=1, alpha=alpha)
        plt.xlabel("Mirrored BAF")
        plt.ylabel("Read-depth ratio")

        if ylim is not None:
            plt.ylim(ylim)
            ymin, ymax = ylim
        else:
            ymin = df.RD.min()
            ymax = df.RD.max()

        if xlim is not None:
            xlim = xlim[0], min(df.BAF.max() + 0.02, xlim[1])
            plt.xlim(xlim)

        plt.vlines(
            ymin=ymin,
            ymax=ymax,
            x=0.5,
            linestyle=":",
            linewidth=1,
            colors="grey",
        )

        if show_centers:
            my_seg = seg[seg.SAMPLE == s]
            plt.scatter(
                my_seg.BAF,
                my_seg.RD,
                c=my_seg["#ID"],
                cmap=mymap,
                s=10,
                linewidth=1,
                edgecolor="k",
            )

            for _, r in my_seg.iterrows():
                plt.annotate(r["#ID"], (r.BAF, r.RD))

        plt.title(s)
        plt.tight_layout()

        if outdir is not None:
            plt.savefig(os.path.join(outdir, f"2D_{s}.png"))

        if not display:
            plt.close()


def plot_1d(
    bbc,
    baf_lim=None,
    rdr_lim=None,
    display=False,
    outdir=None,
    alpha=1,
    show_centromeres=False,
):
    # Prepend 'chr' to #CHR column if not already present
    using_chr = [str(a).startswith("chr") for a in bbc["#CHR"].unique()]
    if any(using_chr):
        if not all(using_chr):
            raise ValueError(
                sp.error("Some chromosomes use 'chr' notation while others do not.")
            )
        use_chr = True
    else:
        use_chr = False
    if not use_chr:
        bbc["#CHR"] = "chr" + bbc["#CHR"].astype("str")

    ### code below assumes that chromosomes in BBC table are named with prefix 'chr' ###
    chrlengths = {str(c): df.END.max() for c, df in bbc.groupby("#CHR")}
    chr_ends = [0]
    for i in range(22):
        chr_ends.append(chr_ends[-1] + chrlengths.get(f"chr{i + 1}", 0))

    # NOTE: this implementation assumes that the only gaps between bins are centromeres
    # If this is not the case in the future this needs to be updated
    chr2centro = {}
    ss_bbc = bbc[bbc.SAMPLE == bbc.iloc[0].SAMPLE]
    break_idx = np.where(
        np.logical_and(
            ss_bbc.END.to_numpy()[:-1] != ss_bbc.START.to_numpy()[1:],
            ss_bbc["#CHR"].to_numpy()[:-1] == ss_bbc["#CHR"].to_numpy()[1:],
        )
    )[0]
    chr2centro = {}
    for i in break_idx:
        rows = ss_bbc.iloc[i : i + 2]
        chr2centro[rows.iloc[0]["#CHR"]] = (
            rows.iloc[0].END,
            rows.iloc[1].START,
        )

    for s, df in bbc.groupby("SAMPLE"):
        plt.figure(figsize=(8, 5.6))
        plt.title(s)
        plt.subplot(211)
        plot_track(
            df,
            chr_ends,
            chr2centro,
            yval="RD",
            ylabel="Read-depth ratio",
            color_field="CLUSTER",
            ylim=rdr_lim,
            display=display,
            alpha=alpha,
            title=s,
            show_centromeres=show_centromeres,
        )
        plt.gca().grid(False)

        plt.subplot(212)
        plot_track(
            df,
            chr_ends,
            chr2centro,
            yval="BAF",
            ylabel="mhBAF",
            color_field="CLUSTER",
            ylim=baf_lim,
            display=display,
            alpha=alpha,
            title=s,
            show_centromeres=show_centromeres,
        )
        plt.gca().grid(False)

        if outdir is not None:
            plt.savefig(os.path.join(outdir, f"1D_{s}.png"))

        if not display:
            plt.close()


def plot_track(
    bb,
    chr_ends,
    chr2centro,
    yval="RD",
    ylabel=None,
    display=True,
    ylim=None,
    alpha=1,
    color_field=None,
    title=None,
    show_centromeres=False,
):
    """
    NOTE: this function assumes that
    1) bb contains data for a single sample
    2) chromosomes are specified using "chr" notation
    """
    xs = [
        int(r["START"])
        + chr_ends[int(r["#CHR"][3:]) - 1]
        + 0.5 * (r["END"] - r["START"])
        for _, r in bb.iterrows()
    ]
    ys = bb[yval]

    markersize = int(max(1, 4 - np.floor(len(bb) / 500)))

    xtick_labels = [f"chr{i}" for i in range(1, 23)]
    xtick_locs = [(chr_ends[i] + chr_ends[i + 1]) / 2 for i in range(22)]

    if color_field is not None:
        # TODO: throw specific value error if "color_field" is not a column in "bb"
        plt.scatter(xs, ys, s=markersize, c=bb[color_field], alpha=alpha, cmap=mymap)
    else:
        plt.scatter(xs, ys, s=markersize, alpha=alpha)

    minx = np.min(xs)
    maxx = np.max(xs)

    if ylim is not None:
        plt.ylim(ylim)
        miny = ylim[0]
        maxy = ylim[1]
    else:
        miny = np.min(ys)
        maxy = np.max(ys)

    plt.vlines(x=chr_ends[1:-1], ymin=miny, ymax=maxy, colors="k", linewidth=0.5)
    if yval == "BAF":
        plt.hlines(
            y=0.5,
            xmin=minx,
            xmax=maxx,
            colors="grey",
            linestyle=":",
            linewidth=1,
        )

    if show_centromeres:
        chromosomes = sorted(bb["#CHR"].unique(), key=lambda x: int(x[3:]))
        ax = plt.gca()
        ax.grid(False)
        for i, ch in enumerate(chromosomes):
            chr_start = chr_ends[i]
            if ch in chr2centro:
                cent_start, cent_end = chr2centro[ch]
                ax.add_patch(
                    Rectangle(
                        xy=(cent_start + chr_start, miny),
                        width=cent_end - cent_start,
                        height=maxy - miny,
                        linewidth=0,
                        color=(0, 0, 0, 0.3),
                    )
                )

    plt.xlim([minx - 1, maxx + 1])
    plt.ylim([miny, maxy])
    plt.xticks(xtick_locs, xtick_labels, rotation=70)

    if title is not None:
        plt.title(title)
    if ylabel is None:
        plt.ylabel(yval)
    else:
        plt.ylabel(ylabel)

    plt.tight_layout()
