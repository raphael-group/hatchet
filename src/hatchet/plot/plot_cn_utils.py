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
from matplotlib.patches import Rectangle


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
    clone_ploidies=None,
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
            has_pi_viol = "PI_VIOL" in bins_seg.columns
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
                        edgecolor="none",
                        transform=ax.get_xaxis_transform(),
                        linewidth=0,
                        antialiased=False,
                    )
                    ax.add_patch(rect)

                    # plot mirrored events: cnb > cna is the mirrored (b, a)
                    # orientation, marked with a "/" hatch over the colored cell.
                    if cnb > cna:
                        ax.add_patch(
                            Rectangle(
                                (x0, y0),
                                w,
                                h,
                                facecolor="none",
                                edgecolor=BLACK,
                                hatch="/",
                                linewidth=0,
                                transform=ax.get_xaxis_transform(),
                            )
                        )

                # PI violation indicator: colored line at top of segment
                if has_pi_viol:
                    viol = bool(bins_seg["PI_VIOL"].iloc[bi])
                    edge_color = "#d62728" if viol else "#2ca02c"
                    y_top = num_clones * h
                    ax.plot(
                        [x0, x0 + w],
                        [y_top, y_top],
                        color=edge_color,
                        linewidth=3,
                        solid_capstyle="butt",
                        transform=ax.get_xaxis_transform(),
                    )
            # plot segment bound (centromere) as dashed line
            if si < len(regions_ch) - 1:
                ax.vlines(
                    ch_offset,
                    ymin=0,
                    ymax=1,
                    transform=ax.get_xaxis_transform(),
                    linewidth=0.5,
                    colors=BLACK,
                    linestyles="dashed",
                )
        # add chromosome boundary (skip last chr end)
        if ch != chs[-1]:
            line = ax.vlines(
                ch_offset,
                ymin=0,
                ymax=1.15,
                transform=ax.get_xaxis_transform(),
                linewidth=1,
                colors=BLACK,
            )
            line.set_clip_on(False)
    ch_coords.append(ch_offset)  # genome end

    ax.grid(False)
    ax.set_xlim(0, ch_offset)
    ax.set_xlabel("")
    for spine in ax.spines.values():
        spine.set_visible(False)
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
        prop = round(bulk_props[ci] * 100, 2)
        lines = []
        if show_clone_name:
            lines.append(f"Clone {ci}")
            lines.append(f"({prop}%)")
        else:
            lines.append(str(ci))
        if clone_ploidies is not None:
            clone_key = f"clone{ci}"
            if clone_key in clone_ploidies:
                lines.append(f"ploidy {round(clone_ploidies[clone_key], 2)}")
        if show_prop:
            lines.append(f"prop {prop}%")
        ylabels.append("\n".join(lines))
    ax.set_yticklabels(ylabels, fontsize=8, va="center")
    ax.set_ylim(0, num_clones * h)
    ax.tick_params(axis="y", which="both", left=True, right=False, length=4)

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

    for _total, states in sorted(tcn_states.items()):
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
        -0.2,
        ">7",
        ha="center",
        va="top",
        fontsize=10,
    )
    leg_x = group_x0 + pair_w + gap_groups

    # --- mirrored CNA box: single CN-state-shaped box with a "/" hatch -------
    rect = Rectangle(
        (leg_x, 0.0),
        pair_w,
        pair_h,
        facecolor="white",
        edgecolor="black",
        hatch="/",
    )
    ax.add_patch(rect)
    ax.text(
        leg_x + pair_w / 2.0,
        -0.2,
        "mirrored",
        ha="center",
        va="top",
        fontsize=10,
    )
    leg_x += pair_w + gap_groups

    # main title on the left
    ax.text(-0.5, pair_h / 2.0, "CNA", fontsize=12, fontweight="bold", ha="right", va="center")

    # nice limits + aspect
    ax.set_xlim(-2.0, leg_x)
    ax.set_ylim(-0.8, pair_h + 0.8)
    ax.set_aspect("auto")  # was 'equal' → made legend look too thin

    return ax


BLACK = (0, 0, 0, 1)


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

    default_color = "#00cc99"  # teal-green, distinct from all CN<=7 colors
    state_style = {}
    for (major, minor), color in _palette.items():
        state_style[(major, minor)] = color
        state_style[(minor, major)] = color
    state_style["default"] = default_color
    return state_style, tcn_states


def get_ascn_colors():
    """Return (state_style, tcn_states) for allele-CN coloring.

    cn=0 → white, cn=1 → black, cn=2..6 + 7+ → inferno (orange→pale yellow)
    evenly spaced. Inferno is colorblind-safe and perceptually uniform.
    All non-zero CN colors are drawn at alpha=0.5 at the call site.
    Sample positions: np.linspace(0.65, 0.97, 6) on inferno → cn=2..6, 7+.
    """
    state_style = {
        0: "#FFFFFF",  # white
        1: "#000000",  # black (rendered at alpha=0.5 → mid gray)
        2: "#EA632A",  # inferno 0.65 (orange)
        3: "#F57D15",  # inferno 0.71
        4: "#FB9B06",  # inferno 0.78 (amber)
        5: "#FBBA1F",  # inferno 0.84
        6: "#F5D949",  # inferno 0.91 (yellow)
    }
    state_style["default"] = "#F3F68A"  # inferno 0.97 (pale yellow) for cn >= 7
    tcn_states = sorted(k for k in state_style if isinstance(k, int))
    return state_style, tcn_states


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
    clone_ploidies=None,
):
    """
    Plot allele-specific CN profile with two sub-bars (A/B) per clone.

    Same genome layout as plot_cnv_profile but each clone row is split into
    an upper A-allele bar and a lower B-allele bar, colored by individual
    integer copy number.
    """
    state_style, _ = get_ascn_colors()
    num_clones = len(str(bin_info.iloc[0]["CNP"]).split(";")) - 1
    h = height / num_clones
    clone_gap = 0.10 * h  # vertical gap between adjacent clone slots
    h_pair = h - clone_gap  # height taken by A+B pair within a slot
    h_sub = h_pair / 2
    y_gap = clone_gap / 2  # padding above & below A/B pair

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
                clone_states = [
                    (int(cn.split("|")[0]), int(cn.split("|")[1])) for cn in bin_cnvs
                ]
                # Mirrored LOH: every clone state is LOH (a==0 or b==0),
                # AND the bin contains BOTH an A-LOH clone (a>0,b=0) and a
                # B-LOH clone (a=0,b>0).
                any_non_loh = any(a > 0 and b > 0 for a, b in clone_states)
                dirs = [
                    (1 if (a > 0 and b == 0) else (-1 if (a == 0 and b > 0) else 0))
                    for a, b in clone_states
                ]
                has_mirror = (not any_non_loh) and (1 in dirs) and (-1 in dirs)
                for k in range(num_clones):
                    cna, cnb = clone_states[num_clones - k - 1]
                    direction = dirs[num_clones - k - 1]
                    y_b = k * h + y_gap  # B allele (bottom half within slot)
                    y_a = y_b + h_sub  # A allele (top half within slot)
                    rect_b = Rectangle(
                        (x0, y_b),
                        w,
                        h_sub,
                        facecolor=state_style.get(cnb, state_style["default"]),
                        edgecolor="none",
                        transform=ax.get_xaxis_transform(),
                        linewidth=0,
                        alpha=1.0 if cnb == 0 else 0.5,
                    )
                    ax.add_patch(rect_b)
                    rect_a = Rectangle(
                        (x0, y_a),
                        w,
                        h_sub,
                        facecolor=state_style.get(cna, state_style["default"]),
                        edgecolor="none",
                        transform=ax.get_xaxis_transform(),
                        linewidth=0,
                        alpha=1.0 if cna == 0 else 0.5,
                    )
                    ax.add_patch(rect_a)

                    # Mirrored-LOH symbol: ≪ / ≫ — two horizontal chevrons placed
                    # side-by-side, spanning the full clone height.
                    # A-LOH → ≫ (apex right), B-LOH → ≪ (apex left).
                    if has_mirror and direction != 0:
                        n_chev = 2
                        chev_unit = w * 0.12
                        gap = w * 0.04
                        total_w = n_chev * chev_unit + (n_chev - 1) * gap
                        x_start = x0 + (w - total_w) / 2.0
                        y_high = y_b + 2 * h_sub
                        y_low = y_b
                        y_mid = (y_high + y_low) / 2.0
                        for i in range(n_chev):
                            cx_left = x_start + i * (chev_unit + gap)
                            cx_right = cx_left + chev_unit
                            if direction > 0:
                                xs = [cx_left, cx_right, cx_left]
                            else:
                                xs = [cx_right, cx_left, cx_right]
                            ax.plot(
                                xs,
                                [y_high, y_mid, y_low],
                                color="black",
                                linewidth=1.2,
                                alpha=0.5,
                                solid_capstyle="round",
                                transform=ax.get_xaxis_transform(),
                            )

            if si < len(regions_ch) - 1:
                # Dashed centromere line, drawn per-clone so it doesn't cross gaps
                for k in range(num_clones):
                    ax.vlines(
                        ch_offset,
                        ymin=k * h + y_gap,
                        ymax=k * h + y_gap + h_pair,
                        transform=ax.get_xaxis_transform(),
                        linewidth=0.5,
                        colors=BLACK,
                        linestyles="dashed",
                    )
        if ch != chs[-1]:
            line = ax.vlines(
                ch_offset,
                ymin=0,
                ymax=1.15,
                transform=ax.get_xaxis_transform(),
                linewidth=1,
                colors=BLACK,
            )
            line.set_clip_on(False)
    ch_coords.append(ch_offset)

    # Outline each clone's A and B rows with a single black border
    for k in range(num_clones):
        y_b_k = k * h + y_gap
        y_a_k = y_b_k + h_sub
        for y0 in (y_b_k, y_a_k):
            ax.add_patch(
                Rectangle(
                    (0, y0),
                    ch_offset,
                    h_sub,
                    facecolor="none",
                    edgecolor="black",
                    linewidth=0.5,
                    transform=ax.get_xaxis_transform(),
                )
            )

    ax.grid(False)
    ax.set_xlim(0, ch_offset)
    ax.set_xlabel("")
    for spine in ax.spines.values():
        spine.set_visible(False)
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
        lines = []
        if show_clone_name:
            lines.append(f"Clone {ci}")
        else:
            lines.append(str(ci))
        if clone_ploidies is not None:
            clone_key = f"clone{ci}"
            if clone_key in clone_ploidies:
                lines.append(f"ploidy {round(clone_ploidies[clone_key], 2)}")
        if show_prop:
            lines.append(f"prop {prop}%")
        ylabels.append("\n".join(lines))
    ax.set_yticklabels(ylabels, fontsize=8, va="center")

    # A/B sub-labels via minor ticks (centered on each sub-bar)
    minor_positions = []
    minor_labels = []
    for k in range(num_clones):
        minor_positions.append(k * h + y_gap + h_sub * 0.5)
        minor_labels.append("B")
        minor_positions.append(k * h + y_gap + h_sub * 1.5)
        minor_labels.append("A")
    ax.set_yticks(minor_positions, minor=True)
    ax.set_yticklabels(minor_labels, minor=True, fontsize=6)
    ax.tick_params(axis="y", which="minor", left=False, right=False, pad=2)

    ax.set_ylim(0, num_clones * h)
    ax.tick_params(axis="y", which="major", left=True, right=False, length=4, pad=20)

    if ylabel is not None:
        ax.set_ylabel(ylabel, rotation=0, ha="right", va="center")
    if title is not None:
        ax.set_title(title)
    return ax


def plot_ascn_legend(
    ax: plt.Axes,
    box_w: float = 1.2,
    box_h: float = 0.4,
    tick_len: float = 0.08,
    label_fontsize: int = 12,
):
    """Draw a horizontal color bar legend for allele CN values."""
    state_style, tcn_states = get_ascn_colors()
    boxes = list(tcn_states) + ["7+"]
    ax.axis("off")
    x0 = 0.0

    for i, label in enumerate(boxes):
        color = state_style["default"] if label == "7+" else state_style[label]
        rect = Rectangle(
            (x0 + i * box_w, 0.0),
            box_w,
            box_h,
            facecolor=color,
            edgecolor="black",
            alpha=1.0 if label == 0 else 0.5,
        )
        ax.add_patch(rect)
        # tick mark below each box
        xc = x0 + i * box_w + box_w / 2.0
        ax.plot([xc, xc], [-tick_len, 0.0], color="black", linewidth=0.8)
        ax.text(
            xc,
            -tick_len - 0.04,
            str(label),
            ha="center",
            va="top",
            fontsize=label_fontsize,
            fontweight="bold",
        )

    total_w = len(boxes) * box_w
    ax.text(
        -0.3,
        box_h / 2.0,
        "Allele copy number",
        fontsize=label_fontsize,
        fontweight="bold",
        ha="right",
        va="center",
    )

    # Mirrored-LOH swatch: ≫ inside a small bordered box.
    swatch_w = box_w * 0.7
    chev_box_x = total_w + 1.0
    ax.add_patch(
        Rectangle(
            (chev_box_x, 0.0),
            swatch_w,
            box_h,
            facecolor="white",
            edgecolor="black",
        )
    )
    n_chev = 2
    chev_unit = swatch_w * 0.30
    gap = swatch_w * 0.10
    chev_total_w = n_chev * chev_unit + (n_chev - 1) * gap
    chev_x_start = chev_box_x + (swatch_w - chev_total_w) / 2.0
    y_high = box_h
    y_low = 0.0
    y_mid = (y_high + y_low) / 2.0
    for i in range(n_chev):
        cx_left = chev_x_start + i * (chev_unit + gap)
        cx_right = cx_left + chev_unit
        ax.plot(
            [cx_left, cx_right, cx_left],
            [y_high, y_mid, y_low],
            color="black",
            linewidth=1.2,
            alpha=0.5,
            solid_capstyle="round",
        )
    ax.text(
        chev_box_x + swatch_w / 2.0,
        -tick_len - 0.04,
        "Mirrored LOH",
        ha="center",
        va="top",
        fontsize=label_fontsize,
        fontweight="bold",
    )
    ax.set_xlim(-2.0, chev_box_x + swatch_w + 0.5)
    ax.set_ylim(-0.5, box_h + 0.2)
    ax.set_aspect("auto")
    return ax
