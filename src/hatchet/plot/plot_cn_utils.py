"""
Plot integer CN 1D heatmap and legend.
Only segments within predefined whitelist segments will be plotted
e.g., centromere regions are excluded in predefined segments and will
be marked as dotted line.

HATCHet will never only a CNV segments that interleaved with the complement
region of predefined segments.
"""

import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon


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
                    color = state_style.get((cna, cnb), state_style["default"])
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

                    # plot mirrored events: cnb > cna means B allele has more copies.
                    # The right-pointing triangle marks the (a, b) orientation.
                    if cnb > cna:
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

    if ylabel is not None:
        ax.set_ylabel(ylabel, rotation=0, ha="right", va="center")
    if title is not None:
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
        (group_x0, 0.0),
        pair_w,
        pair_h,
        facecolor=default_color,
        edgecolor="black",
    )
    ax.add_patch(rect)
    ax.text(
        group_x0 + pair_w / 2.0,
        pair_h + 0.1,
        "Total CN>7",
        ha="center",
        va="bottom",
        fontsize=10,
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
    _palette = {
        (0, 0): "darkblue",
        (1, 0): "lightblue",
        (1, 1): "lightgray",
        (2, 0): "dimgray",
        (2, 1): "lightgoldenrodyellow",
        (3, 0): "gold",
        (2, 2): "navajowhite",
        (3, 1): "orange",
        (4, 0): "darkorange",
        (3, 2): "salmon",
        (4, 1): "red",
        (5, 0): "darkred",
        (3, 3): "plum",
        (4, 2): "orchid",
        (5, 1): "purple",
        (6, 0): "indigo",
        (4, 3): "#c0b7f0",
        (5, 2): "#a485f4",
        (6, 1): "#6f42c1",
        (7, 0): "#4b0082",
    }

    copy_states = [
        (1, 0), (0, 1),
        (0, 2), (1, 1), (2, 0),
        (0, 3), (1, 2), (2, 1), (3, 0),
        (0, 4), (1, 3), (2, 2), (3, 1), (4, 0),
        (0, 5), (1, 4), (2, 3), (3, 2), (4, 1), (5, 0),
        (0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0),
        (0, 7), (1, 6), (2, 5), (3, 4), (4, 3), (5, 2), (6, 1), (7, 0),
    ]
    tcn_states = {}
    for a, b in copy_states:
        tcn_states.setdefault(int(a + b), []).append((a, b))

    default_color = "#333333"
    state_style = {}
    for (major, minor), color in _palette.items():
        state_style[(major, minor)] = color
        state_style[(minor, major)] = color
    state_style["default"] = default_color
    return state_style, tcn_states


_ASCN_COLORS = [
    "#2166AC",  # CN=0  blue (homozygous deletion)
    "#A8A8A8",  # CN=1  grey (neutral)
    "#FDBE85",  # CN=2  light peach
    "#FD8D3C",  # CN=3  light orange
    "#F46D43",  # CN=4  orange
    "#E8352E",  # CN=5  red
    "#CC1620",  # CN=6  medium-dark red
    "#B0000E",  # CN=7  dark red
    "#8B0000",  # CN=8  deeper red
    "#660000",  # CN=9  very dark red
    "#400000",  # CN=10 darkest red
]


def ascn_color(cn: int) -> str:
    """Return a hex color for an integer allele copy number (0–10)."""
    return _ASCN_COLORS[max(0, min(cn, 10))]


def plot_ascn_profile(
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
    Plot allele-specific CN profile with two sub-bars (A/B) per clone.

    Same genome layout as plot_cnv_profile but each clone row is split into
    an upper A-allele bar and a lower B-allele bar, colored by individual
    integer copy number on a blue→grey→red scale.
    """
    num_clones = len(str(bin_info.iloc[0]["CNP"]).split(";")) - 1
    h = height / num_clones
    h_sub = h / 2

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
            wl_start = wl_segment["START"]
            wl_end = wl_segment["END"]
            seg_end = ch_offset + (wl_end - wl_start)

            bins_seg = bins_ch.loc[
                (bins_ch["START"] >= wl_start) & (bins_ch["END"] <= wl_end), :
            ]
            if bins_seg.empty:
                ch_offset = seg_end
                continue

            bin_starts = (bins_seg["START"] - wl_start + ch_offset).to_numpy()
            bin_ends = (bins_seg["END"] - wl_start + ch_offset).to_numpy()
            ch_offset = seg_end
            if si < len(regions_ch) - 1:
                seg_coords.append(ch_offset)

            for bi in range(len(bins_seg)):
                x0, bin_end = bin_starts[bi], bin_ends[bi]
                w = bin_end - x0
                bin_cnvs = bins_seg["CNP"].iloc[bi].split(";")[1:]
                for k in range(num_clones):
                    clone_cn = bin_cnvs[num_clones - k - 1]
                    cna_str, cnb_str = clone_cn.split("|")
                    cna, cnb = int(cna_str), int(cnb_str)
                    y0 = k * h
                    # B allele (bottom half)
                    rect_b = Rectangle(
                        (x0, y0),
                        w,
                        h_sub,
                        facecolor=ascn_color(cnb),
                        edgecolor=BLACK,
                        transform=ax.get_xaxis_transform(),
                        linewidth=0,
                    )
                    ax.add_patch(rect_b)
                    # A allele (top half)
                    rect_a = Rectangle(
                        (x0, y0 + h_sub),
                        w,
                        h_sub,
                        facecolor=ascn_color(cna),
                        edgecolor=BLACK,
                        transform=ax.get_xaxis_transform(),
                        linewidth=0,
                    )
                    ax.add_patch(rect_a)

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
        ax.vlines(
            ch_offset,
            ymin=0,
            ymax=1,
            transform=ax.get_xaxis_transform(),
            linewidth=1,
            colors=BLACK,
        )
    ch_coords.append(ch_offset)

    # clone separation lines
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

    # Y-axis: clone name at center, A/B sub-labels
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

    # A/B sub-labels via minor ticks
    minor_positions = []
    minor_labels = []
    for k in range(num_clones):
        minor_positions.append(k * h + h_sub * 0.5)
        minor_labels.append("B")
        minor_positions.append(k * h + h_sub * 1.5)
        minor_labels.append("A")
    ax.set_yticks(minor_positions, minor=True)
    ax.set_yticklabels(minor_labels, minor=True, fontsize=6)
    ax.tick_params(axis="y", which="minor", left=False, right=False, pad=2)

    ax.set_ylim(0, num_clones * h)
    ax.tick_params(axis="y", which="major", left=False, right=False)

    if ylabel is not None:
        ax.set_ylabel(ylabel, rotation=0, ha="right", va="center")
    if title is not None:
        ax.set_title(title)
    return ax


def plot_ascn_legend(ax: plt.Axes):
    """Draw a horizontal color bar legend for allele CN 0–10."""
    ax.axis("off")

    box_w = 2.0
    box_h = 0.6
    x0 = 0.0

    for cn in range(11):
        rect = Rectangle(
            (x0 + cn * box_w, 0.0),
            box_w,
            box_h,
            facecolor=ascn_color(cn),
            edgecolor="black",
        )
        ax.add_patch(rect)
        ax.text(
            x0 + cn * box_w + box_w / 2.0,
            -0.15,
            str(cn),
            ha="center",
            va="top",
            fontsize=9,
        )

    total_w = 11 * box_w
    ax.text(
        -0.5, box_h / 2.0, "Allele copy number", fontsize=12, ha="right", va="center"
    )

    ax.set_xlim(-2.0, total_w + 1.0)
    ax.set_ylim(-0.6, box_h + 0.4)
    ax.set_aspect("auto")
    return ax


def make_cnp_palette(clone_states, clone_props, state_style, style="ascn"):
    """Map each CNP string to a color blended by clone proportions.

    Args:
        clone_states: iterable of unique CNP strings (e.g. ["1|1;2|1", "1|1;1|1"]).
        clone_props:  (n_clones,) proportions array (normal clone first).
        state_style:  dict (a, b) -> hex color, from get_cn_colors().
        style:        ``"ascn"`` to color tumor clones by total CN via
                      ascn_color(a+b), re-normalised to exclude the normal
                      clone; ``"cnv"`` to blend all clones using state_style.

    Returns:
        dict mapping cnp_string -> RGB tuple (values in [0, 1]).
    """
    from matplotlib.colors import to_rgb

    palette = {}
    for cnp in clone_states:
        states = [(int(x.split("|")[0]), int(x.split("|")[1])) for x in cnp.split(";")]
        if style == "ascn":
            # Blend only tumor clones (index >= 1), re-normalised proportions
            tumor_props = clone_props[1:]
            tumor_sum = tumor_props.sum()
            if tumor_sum > 0:
                tumor_props_norm = tumor_props / tumor_sum
            else:
                tumor_props_norm = np.ones(len(tumor_props)) / len(tumor_props)
            rgbs = np.array([to_rgb(ascn_color(a + b)) for a, b in states[1:]])
            blended = np.clip(rgbs.T @ tumor_props_norm, 0.0, 1.0)
        else:
            # Blend only tumor clones (index >= 1) to avoid dilution by normal gray
            tumor_props = clone_props[1:]
            tumor_sum = tumor_props.sum()
            if tumor_sum > 0:
                tumor_props_norm = tumor_props / tumor_sum
            else:
                tumor_props_norm = np.ones(len(tumor_props)) / len(tumor_props)
            rgbs = np.array(
                [
                    to_rgb(state_style.get((a, b), state_style["default"]))
                    for a, b in states[1:]
                ]
            )
            blended = np.clip(rgbs.T @ tumor_props_norm, 0.0, 1.0)
        palette[cnp] = tuple(blended)
    return palette
