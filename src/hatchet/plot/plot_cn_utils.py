import os
import sys

import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Rectangle, Polygon

"""
Plot integer CN 1D heatmap and legend.
Only segments within predefined whitelist segments will be plotted
e.g., centromere regions are excluded in predefined segments and will
be marked as dotted line.

HATCHet will never only a CNV segments that interleaved with the complement
region of predefined segments.
"""


def plot_cnv_profile(
    ax: plt.Axes,
    bin_info: pd.DataFrame,
    regions: pd.DataFrame,
    width=20,
    height=1,
    title=None,
    ylabel=None,
    plot_chrname=True,
    show_prop=True,
    show_clone_name=True,
):
    """
    plot chrom-level integer CNV profile.
    regions outside white list segments are marked as dashed bar.
    bin_info: #CHR, START, END, CNP
    """
    state_style, tcn_states = get_cn_colors()

    num_clones = (
        len(str(bin_info.iloc[0]["CNP"]).split(";")) - 1
    )  # first column is normal
    h = height / num_clones

    bulk_props = np.array([float(v) for v in str(bin_info["PROPS"].iloc[0]).split(";")])

    regions_chs = regions.groupby(by="#CHR", sort=False)
    bins_chs = bin_info.groupby(by="#CHR", sort=False, observed=True)

    ch_offset = 0
    ch_coords = []
    seg_coords = []
    chs = bin_info["#CHR"].unique()
    for ch in chs:
        ch_coords.append(ch_offset)
        regions_ch = regions_chs.get_group(ch)
        bins_ch = bins_chs.get_group(ch)
        for si in range(len(regions_ch)):
            wl_segment = regions_ch.iloc[si]
            seg_start = ch_offset
            wl_start = wl_segment["START"]
            wl_end = wl_segment["END"]
            seg_end = ch_offset + (wl_end - wl_start)

            bins_seg = bins_ch.loc[
                (bins_ch["START"] >= wl_start) & (bins_ch["END"] <= wl_end), :
            ]
            if bins_seg.empty:
                ch_offset = seg_end
                continue

            # global bin coords
            bin_starts = (bins_seg["START"] - wl_start + ch_offset).to_numpy()
            bin_ends = (bins_seg["END"] - wl_start + ch_offset).to_numpy()
            # update global offsets
            ch_offset = seg_end
            if si < len(regions_ch) - 1:
                seg_coords.append(ch_offset)  # centromere offset

            # plot rectangles
            for bi in range(len(bins_seg)):
                x0, bin_end = bin_starts[bi], bin_ends[bi]
                w = bin_end - x0
                bin_cnvs = bins_seg["CNP"].iloc[bi].split(";")[1:]  # ignore normal cnp
                for k in range(num_clones):
                    cna, cnb = (
                        int(bin_cnvs[num_clones - k - 1].split("|")[0]),
                        int(bin_cnvs[num_clones - k - 1].split("|")[1]),
                    )
                    color = state_style.get((cna, cnb), ("white", None))
                    y0 = k * h
                    rect = Rectangle(
                        (x0, y0),
                        w,
                        h,
                        facecolor=color,
                        edgecolor=BLACK,
                        transform=ax.get_xaxis_transform(),
                        linewidth=0,
                    )
                    ax.add_patch(rect)

                    # plot mirrored events
                    if cnb > cna:
                        if cna < cnb:  # right triangle
                            trig = Polygon(
                                [
                                    [x0, y0 + h],
                                    [x0 + w, y0 + h / 2],
                                    [x0, y0],
                                ],  # top-left → mid-right → bottom-left
                                linewidth=0,
                                closed=True,
                                facecolor=BLACK,
                                transform=ax.get_xaxis_transform(),
                            )
                        else:
                            trig = Polygon(
                                [
                                    [x0 + w, y0 + h],
                                    [x0, y0 + h / 2],
                                    [x0 + w, y0],
                                ],  # top-left → mid-right → bottom-left
                                linewidth=0,
                                closed=True,
                                facecolor=BLACK,
                                transform=ax.get_xaxis_transform(),
                            )
                        ax.add_patch(trig)
            # plot segment bound (centromere) as dashed line
            if si < len(regions_ch) - 1:
                ax.vlines(
                    ch_offset,
                    ymin=0,
                    ymax=1,
                    transform=ax.get_xaxis_transform(),
                    linewidth=1,
                    colors=BLACK,
                    linestyles="dashed",
                )
        # add chromosome boundary
        ax.vlines(
            ch_offset,
            ymin=0,
            ymax=1,
            transform=ax.get_xaxis_transform(),
            linewidth=1,
            colors=BLACK,
        )
    ch_coords.append(ch_offset)  # genome end

    # plot clone separation
    if num_clones > 1:
        ax.hlines(
            y=[h * (i + 1) for i in range(num_clones - 1)],
            xmin=0,
            xmax=ch_offset,
            colors=BLACK,
            linewidth=1,
            transform=ax.get_xaxis_transform(),
        )
    ax.grid(False)
    ax.set_xlim(0, ch_offset)
    ax.set_xlabel("")
    if plot_chrname:
        ax.set_xticks(
            [
                ch_coords[i] + (ch_coords[i + 1] - ch_coords[i]) // 2
                for i in range(len(ch_coords) - 1)
            ]
        )
        ax.set_xticklabels(chs, rotation=60, fontsize=8)
        ax.tick_params(
            axis="x", labeltop=True, labelbottom=False, top=False, bottom=False
        )
    else:
        ax.set_xticks([])
        ax.set_xticklabels([])

    # plot clone name
    ax.set_yticks([h * (i + 0.5) for i in range(num_clones)])
    ylabels = []
    for ci in range(num_clones, 0, -1):
        prop = round(bulk_props[ci] * 100, 1)
        cname = str(ci)
        if show_clone_name:
            cname = f"clone {cname}"
        if show_prop:
            cname = f"{cname} ({prop}%)"
        ylabels.append(cname)
    ax.set_yticklabels(ylabels)
    ax.set_ylim(0, num_clones * h)
    ax.tick_params(axis="y", which="both", left=False, right=False)

    if not ylabel is None:
        ax.set_ylabel(ylabel, rotation=0, ha="right", va="center")
    if not title is None:
        ax.set_title(title)
    return ax


def plot_cnv_legend(ax: plt.Axes):
    state_style, tcn_states = get_cn_colors()
    ax.axis("off")

    pair_w = 2  # width of each (a,b) box
    pair_h = 0.6  # height of each box
    gap_pairs = 0.0  # horizontal gap between boxes in same group
    gap_groups = 1.5  # horizontal gap between total-CN groups

    leg_x = 0.0  # running x position (data coords)

    for total, states in sorted(tcn_states.items()):
        # merge mirrored pairs like (1,2),(2,1) -> (2,1)
        uniq_pairs = sorted({tuple(sorted(s, reverse=True)) for s in states})
        n_pairs = len(uniq_pairs)

        group_w = n_pairs * pair_w + (n_pairs - 1) * gap_pairs
        group_x0 = leg_x

        # draw boxes and per-pair labels
        for i, pair in enumerate(uniq_pairs):
            x0 = group_x0 + i * (pair_w + gap_pairs)
            color = state_style[pair]

            rect = Rectangle(
                (x0, 0.0),
                pair_w,
                pair_h,
                facecolor=color,
                edgecolor="black",
            )
            ax.add_patch(rect)

            ax.text(
                x0 + pair_w / 2.0,
                -0.2,
                f"{pair}",
                ha="center",
                va="top",
                fontsize=10,
            )

        # group title above the boxes
        ax.text(
            group_x0 + group_w / 2.0,
            pair_h + 0.1,
            f"Total CN={total}",
            ha="center",
            va="bottom",
            fontsize=10,
        )

        leg_x = group_x0 + group_w + gap_groups

    # --- CN>7 default box ----------------------------------------------------
    default_color = state_style["default"]
    group_x0 = leg_x
    rect = Rectangle(
        (group_x0, 0.0), pair_w, pair_h,
        facecolor=default_color, edgecolor="black",
    )
    ax.add_patch(rect)
    ax.text(
        group_x0 + pair_w / 2.0, pair_h + 0.1, "Total CN>7",
        ha="center", va="bottom", fontsize=10,
    )
    leg_x = group_x0 + pair_w + gap_groups

    # --- mirrored CNA symbol on the right ------------------------------------
    # leave a clear gap after the last CN group
    leg_x += 3.0

    mirror_y0 = 0.0
    mirror_w = 1.25

    # label text
    ax.text(
        leg_x,
        mirror_y0 + mirror_w / 2.0,
        "Mirrored CNA",
        ha="left",  # text grows to the right, away from CN=7 group
        va="center",
        fontsize=10,
    )

    # start symbol a bit to the right of the text
    leg_x += 4.0  # adjust if you change fontsize

    # right triangle box (a,b)
    rect = Rectangle(
        (leg_x, mirror_y0),
        mirror_w,
        mirror_w,
        facecolor="white",
        edgecolor="black",
    )
    ax.add_patch(rect)
    trig = Polygon(
        [
            [leg_x, mirror_y0 + mirror_w],
            [leg_x + mirror_w, mirror_y0 + mirror_w - mirror_w / 4.0],
            [leg_x, mirror_y0 + mirror_w - mirror_w / 2.0],
        ],
        linewidth=0,
        closed=True,
        facecolor="black",
    )
    ax.add_patch(trig)
    ax.text(
        leg_x + mirror_w / 2.0,
        mirror_y0 - 0.2,
        "(a,b)",
        ha="center",
        va="top",
        fontsize=10,
    )
    leg_x += mirror_w + 0.5

    # left triangle box (b,a)
    rect = Rectangle(
        (leg_x, mirror_y0),
        mirror_w,
        mirror_w,
        facecolor="white",
        edgecolor="black",
    )
    ax.add_patch(rect)
    trig = Polygon(
        [
            [leg_x, mirror_y0 + mirror_w / 4.0],
            [leg_x + mirror_w, mirror_y0 + mirror_w - mirror_w / 2.0],
            [leg_x + mirror_w, mirror_y0],
        ],
        linewidth=0,
        closed=True,
        facecolor="black",
    )
    ax.add_patch(trig)
    ax.text(
        leg_x + mirror_w / 2.0,
        mirror_y0 - 0.2,
        "(b,a)",
        ha="center",
        va="top",
        fontsize=10,
    )
    leg_x += mirror_w + 3.0

    # main title on the left
    ax.text(-0.5, pair_h / 2.0, "Copy numbers", fontsize=12, ha="right", va="center")

    # nice limits + aspect
    ax.set_xlim(-2.0, leg_x)
    ax.set_ylim(-0.8, pair_h + 0.8)
    ax.set_aspect("auto")  # was 'equal' → made legend look too thin

    return ax


WHITE = (1, 1, 1, 1)
BLACK = (0, 0, 0, 1)
RED = (1, 0, 0, 1)
BLUE = (0, 0, 1, 1)


def fcn_color(fcn: float, base=1, tol=1e-2):
    if fcn > base + tol:
        return RED
    elif fcn < base - tol:
        return BLUE
    else:
        return BLACK


def get_cn_colors():
    copy_states = [
        (1, 0),
        (0, 1),
        (0, 2),
        (1, 1),
        (2, 0),
        (0, 3),
        (1, 2),
        (2, 1),
        (3, 0),
        (0, 4),
        (1, 3),
        (2, 2),
        (3, 1),
        (4, 0),
        (0, 5),
        (1, 4),
        (2, 3),
        (3, 2),
        (4, 1),
        (5, 0),
        (0, 6),
        (1, 5),
        (2, 4),
        (3, 3),
        (4, 2),
        (5, 1),
        (6, 0),
        (0, 7),
        (1, 6),
        (2, 5),
        (3, 4),
        (4, 3),
        (5, 2),
        (6, 1),
        (7, 0),
    ]
    tcn_states = {}
    for a, b in copy_states:
        tcn_states.setdefault(int(a + b), []).append((a, b))

    fixed_colors = {
        1: ["#4a9bc4"],  # medium blue (loss)
        2: ["#a8a8a8", "#505050"],  # medium greys (diploid)
        3: ["#c8c000", "#a07800"],  # golden yellows (3-copy)
        4: ["#e09050", "#c06800", "#904000"],  # ambers/oranges (4-copy)
        5: ["#e84030", "#b00000", "#700000"],  # reds (5-copy)
        6: ["#cc80cc", "#a030a0", "#680068", "#400050"],  # purples
        7: [
            "#9080e0",
            "#6848c8",
            "#4c20a0",
            "#340070",
            "#1c0040",
        ],  # violet → indigo gradient
    }

    default_color = "#333333"
    state_style = {}
    for total, states in tcn_states.items():
        uniq_pairs = sorted({tuple(sorted(s, reverse=True)) for s in states})
        colors = fixed_colors[total]
        for pair, color in zip(uniq_pairs, colors):
            state_style[pair] = color
            state_style[pair[::-1]] = color
    state_style["default"] = default_color
    return state_style, tcn_states


def make_cnp_palette(clone_states, clone_props, state_style):
    """Map each CNP string to a color blended by clone proportions.

    Args:
        clone_states: iterable of unique CNP strings (e.g. ["1|1;2|1", "1|1;1|1"]).
        clone_props:  (n_clones,) proportions array (normal clone first).
        state_style:  dict (a, b) -> hex color, from get_cn_colors().

    Returns:
        dict mapping cnp_string -> RGB tuple (values in [0, 1]).
    """
    from matplotlib.colors import to_rgb

    palette = {}
    for cnp in clone_states:
        states = [
            (int(x.split("|")[0]), int(x.split("|")[1]))
            for x in cnp.split(";")
        ]
        rgbs = np.array([
            to_rgb(state_style.get((a, b), state_style["default"]))
            for a, b in states
        ])  # (n_clones, 3)
        blended = np.clip(rgbs.T @ clone_props, 0.0, 1.0)  # (3,)
        palette[cnp] = tuple(blended)
    return palette
