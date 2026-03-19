import os
import sys
import logging
import contextlib

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from collections import OrderedDict
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.collections import LineCollection
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Ellipse


from scipy.stats import norm, beta as beta_dist, gaussian_kde
from scipy.signal import find_peaks
from adjustText import adjust_text

logging.getLogger("adjustText").setLevel(logging.WARNING)

from hatchet.utils import *
from hatchet.plot.plot_utils import *


##################################################
def get_transparency(
    bbc: pd.DataFrame, by: str, one_tail=0.25, tail_alpha=0.2, nontail_alpha=1.0
) -> pd.Series:
    def get_alpha_one_row(row, per_grp_thres):
        rd, baf = row["RD"], row["BAF"]
        rd_l, rd_r, baf_l, baf_r = per_grp_thres[row[by]]
        if rd <= rd_l or rd >= rd_r or baf <= baf_l or baf >= baf_r:
            return tail_alpha
        else:
            return nontail_alpha

    get_left_threshold = lambda arr: arr[
        min(max(0, int(len(arr) * one_tail)), len(arr) - 1)
    ]
    get_right_threshold = lambda arr: arr[
        min(max(0, int(len(arr) * (1 - one_tail))), len(arr) - 1)
    ]

    per_grp_thres = {}
    for key in bbc[by].unique():
        df = bbc[bbc[by] == key]
        rd = sorted(df["RD"].tolist())
        baf = sorted(df["BAF"].tolist())
        per_grp_thres[key] = [
            get_left_threshold(rd),
            get_right_threshold(rd),
            get_left_threshold(baf),
            get_right_threshold(baf),
        ]

    alphas = bbc.apply(func=lambda r: get_alpha_one_row(r, per_grp_thres), axis=1)
    alphas.name = "transparency"

    return alphas


def get_abs_positions_ignore_gap(
    bin_info: pd.DataFrame,
    regions: pd.DataFrame,
    chrs: list,
    chr_shift=int(10e6),
):
    bin_info["abs_pos"] = 0
    bin_info["abs_start"] = 0
    bin_info["abs_end"] = 0

    ch_offset = chr_shift
    ch_coords = []  # chromosome bounderies
    seg_coords = []  # inter-segment gaps

    regions_chs = regions.groupby(by="#CHR", sort=False)
    bins_chs = bin_info.groupby(by="#CHR", sort=False, observed=True)
    for ch in chrs:
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
                (bins_ch["START"] >= wl_start) & (bins_ch["END"] < wl_end), :
            ]
            if bins_seg.empty:
                ch_offset = seg_end
                continue

            # global bin coords
            bin_starts = (bins_seg["START"] - wl_start + ch_offset).to_numpy()
            bin_ends = (bins_seg["END"] - wl_start + ch_offset).to_numpy()

            bin_info.loc[bins_seg.index, "abs_start"] = bin_starts
            bin_info.loc[bins_seg.index, "abs_end"] = bin_ends
            bin_info.loc[bins_seg.index, "abs_pos"] = (bin_starts + bin_ends) // 2
            # update global offsets
            ch_offset = seg_end
            if si < len(regions_ch) - 1:
                seg_coords.append(ch_offset)  # centromere offset
    ch_coords.append(ch_offset)  # genome end
    axis_start = chr_shift
    axis_end = ch_coords[-1] + chr_shift
    return bin_info, ch_coords, axis_start, axis_end, seg_coords


def get_abs_positions_keep_gap(
    bin_info: pd.DataFrame,
    chrom_sizes: dict,
    chrs: list,
    chr_shift=int(10e6),
):
    """
    chr_bounds: global chromosome boundary
    axis_end:
    """
    chr_offsets = OrderedDict()
    for i, ch in enumerate(chrs):
        if i == 0:
            chr_offsets[ch] = chr_shift
        else:
            prev_ch = chrs[i - 1]
            offset = chr_offsets[prev_ch] + chrom_sizes[prev_ch]
            chr_offsets[ch] = offset
    ch_coords = list(chr_offsets.values()) + [
        chr_offsets[chrs[-1]] + chrom_sizes[chrs[-1]]
    ]
    axis_start = chr_shift
    axis_end = ch_coords[-1] + chr_shift
    bin_info["abs_pos"] = bin_info.apply(
        func=lambda r: chr_offsets[r["#CHR"]] + (r.START + r.END) // 2, axis=1
    ).to_numpy()
    bin_info["abs_start"] = bin_info.apply(
        func=lambda r: chr_offsets[r["#CHR"]] + r.START, axis=1
    ).to_numpy()
    bin_info["abs_end"] = bin_info.apply(
        func=lambda r: chr_offsets[r["#CHR"]] + r.END, axis=1
    ).to_numpy()
    return bin_info, ch_coords, axis_start, axis_end, []


##################################################
def plot_1d(
    ax: plt.Axes,
    sample: str,
    bin_info: pd.DataFrame,
    vals: np.ndarray,
    regions: pd.DataFrame,
    chrom_sizes: dict,
    exp_colname=None,
    val_type="BAF",
    colors=None,
    hue=None,
    palette=None,
    ylim=(0, 1),
    ylab="value",
    title=None,
    plot_chrname=True,
    ignore_gap=True,
    chr_shift=int(10e6),
    markersize=2.0,
    exp_linewidth=1.5,
    bd_linewidth=1,
    linecolor=(0, 0, 0, 1),
    show_legend=True,
    rasterized=True,
):
    """
    plot chrom-level 1D BAF/RDR/FCN values.
    assert bin_info is sorted by chromosome and positions
    If ignore_gap=True, ignore region outside regions (like centromeres)
        and represented as dashed line.
    """
    bin_info = bin_info.copy(deep=True)
    chrs = bin_info["#CHR"].unique().tolist()
    if ignore_gap:
        bin_info, ch_coords, axis_start, axis_end, seg_coords = (
            get_abs_positions_ignore_gap(
                bin_info, regions, chrs, chr_shift=chr_shift
            )
        )
    else:
        bin_info, ch_coords, axis_start, axis_end, seg_coords = (
            get_abs_positions_keep_gap(bin_info, chrom_sizes, chrs, chr_shift)
        )
    ##################################################
    g = sns.scatterplot(
        x=bin_info["abs_pos"],
        y=vals,
        ax=ax,
        s=markersize,
        color=colors,
        hue=hue,
        palette=palette,
        legend=(hue is not None) and show_legend,
    )
    if rasterized:
        for coll in ax.collections:
            coll.set_rasterized(True)
    if val_type == "BAF":
        ax.hlines(
            y=0.5,
            xmin=axis_start,
            xmax=axis_end,
            colors="grey",
            linestyle=":",
            linewidth=bd_linewidth,
        )

    if not exp_colname is None:
        # plot expected values as hlines
        exp_lines = []
        exp_colors = [linecolor] * len(bin_info)
        for _, row in bin_info.iterrows():
            exp_lines.append(
                [
                    (row["abs_start"], row[exp_colname]),
                    (row["abs_end"], row[exp_colname]),
                ]
            )
        ax.add_collection(
            LineCollection(exp_lines, linewidth=exp_linewidth, colors=exp_colors)
        )

    ax.vlines(
        ch_coords,
        ymin=0,
        ymax=1,
        transform=ax.get_xaxis_transform(),
        linewidth=bd_linewidth,
        colors="k",
    )

    ax.set_ylim(ylim[0], ylim[1])
    ax.set_ylabel(ylab)
    ax.set_title(title)
    xticks = [
        ch_coords[i] + (ch_coords[i + 1] - ch_coords[i]) // 2
        for i in range(len(ch_coords) - 1)
    ]
    plt.setp(ax, xlim=(0, axis_end), xticks=xticks, xlabel="")
    if plot_chrname:
        ax.set_xticklabels(chrs, rotation=60, fontsize=8)

    return


##################################################
def plot_2d(
    sample: str,
    bin_info: pd.DataFrame,
    xvals: np.ndarray,
    yvals: np.ndarray,
    exp_xvals=None,
    exp_yvals=None,
    exp_labels=None,
    clone_props=None,
    alphas=None,
    hue=None,
    palette=None,
    xlab="BAF",
    ylab="RDR",
    xlim=None,
    ylim=None,
    title=None,
    label_clone=False,
    markersize=2.0,
    markersize_centroid=10,
    marker_bd_width=0.8,
    dpi=300,
    transparent=False,
    out_file=None,
    rasterized=True,
):
    """2D RDR-vs-BAF joint plot with per-cluster KDE marginals.

    Draws a scatter plot on the joint axes (colored by ``hue``) and
    unfilled KDE curves on the marginal axes — one line per cluster,
    colored to match the scatter palette.
    """
    g0 = sns.JointGrid(x=xvals, y=yvals, hue=hue, palette=palette, xlim=xlim, ylim=ylim)
    g0.refline(x=0.50)
    g0.plot_joint(sns.scatterplot, s=markersize, legend=False, edgecolors="none")
    g0.plot_marginals(sns.kdeplot, common_norm=False, linewidth=0.8, fill=False)
    scatter = g0.ax_joint.collections[0]
    scatter.set_rasterized(rasterized)
    scatter.set_antialiased(False)

    if alphas is not None:
        scatter.set_alpha(None)
        colors_ = scatter.get_facecolors()  # reused in 1D plot
        colors_[:, 3] = np.asarray(alphas)
        scatter.set_facecolors(colors_)
    g0_colors = scatter.get_facecolors()

    if not exp_labels is None:
        texts = []
        for ci, cid in enumerate(exp_labels):
            center_text = cid
            fontdict = {"fontsize": 10}
            if label_clone:
                states = [f"({x[0]},{x[1]})" for x in cid][1:]
                if len(set(states)) == 1:
                    fontdict["fontweight"] = "bold"
                    center_text = states[0]
                else:
                    center_text = ",".join(states)
            else:
                center_text = cid
            t = g0.ax_joint.text(exp_xvals[ci], exp_yvals[ci], center_text, **fontdict)
            texts.append(t)
        with open(os.devnull, "w") as _devnull, contextlib.redirect_stdout(_devnull):
            adjust_text(
                texts,
                x=exp_xvals,
                y=exp_yvals,
                ax=g0.ax_joint,
                arrowprops=dict(arrowstyle="-", color="black", lw=0.5),
            )
        g0.ax_joint.scatter(
            x=exp_xvals,
            y=exp_yvals,
            facecolors="none",
            edgecolors="black",
            s=markersize_centroid,
            linewidth=marker_bd_width,
        )

    if not clone_props is None:
        custom_handles = []
        for i, prop in enumerate(clone_props):
            lab = f"Normal: {prop:.3f}" if i == 0 else f"Clone {i}: {prop:.3f}"
            custom_handles.append(Line2D([0], [0], alpha=0, label=lab))
        g0.ax_joint.legend(
            handles=custom_handles,
            title="",
            loc="best",
            fontsize="small",
            fancybox=True,
            framealpha=0.7,
            handlelength=0,
            handletextpad=0,
        )

    g0.set_axis_labels(xlabel=xlab, ylabel=ylab)
    g0.figure.suptitle(title)
    plt.tight_layout()
    if out_file is not None:
        g0.savefig(out_file, dpi=dpi, bbox_inches="tight", transparent=transparent)
        plt.close(g0.figure)
    return g0.figure, g0_colors


##################################################
def plot_baf_normal_sample(
    bafs: np.ndarray,
    out_file: str,
    dpi=150,
    transparent=False,
):
    fig, ax = plt.subplots(1, 1)
    g0 = sns.histplot(x=bafs, ax=ax, bins=50, binrange=[0, 1])
    ax.vlines(
        0.5,
        ymin=0,
        ymax=1,
        transform=ax.get_xaxis_transform(),
        linewidth=0.5,
        colors="k",
    )
    mu_baf = np.mean(bafs)
    std_baf = np.std(bafs)
    med_baf = np.median(bafs)
    ax.set_title(
        f"Normal #bins={len(bafs)}\nmu={mu_baf:.3f} std={std_baf:.3f} med={med_baf:.3f}"
    )
    ax.set_xlabel(xlabel="mhBAF")
    ax.grid(False)
    plt.tight_layout()
    plt.savefig(out_file, dpi=dpi, transparent=transparent)
    plt.close(fig)
    return


##################################################
def _is_multimodal(obs, min_count=30):
    """Return True if KDE of obs has >1 prominent peak."""
    if len(obs) < min_count:
        return False
    try:
        kde = gaussian_kde(obs)
    except (np.linalg.LinAlgError, ValueError):
        return False
    grid = np.linspace(obs.min(), obs.max(), 200)
    kde_vals = kde(grid)
    peaks, _ = find_peaks(kde_vals, prominence=0.1 * kde_vals.max())
    return len(peaks) > 1


def plot_clusters(
    cluster_labels: np.ndarray,
    X_rdrs: np.ndarray,
    X_bafs: np.ndarray,
    rdr_means: np.ndarray,
    rdr_vars: np.ndarray,
    baf_means: np.ndarray,
    baf_taus: np.ndarray,
    tumor_samples: list,
    out_file: str = None,
    cluster_ids: np.ndarray = None,
    log_rdr: bool = False,
    pdf=None,
):
    """Plot per-cluster joint BAF-vs-RDR scatter with marginals and QQ plots.

    One PDF page per cluster with joint-grid panels (scatter + contour,
    marginal histograms) on top and BAF/RDR QQ plots below.

    Args:
        cluster_labels: (N,) 0-indexed cluster assignments.
        X_rdrs:         (N, M) RDR values (original scale).
        X_bafs:         (N, M) phased BAF values in [0, 1].
        rdr_means:      (K, M) fitted RDR means (log-space if log_rdr).
        rdr_vars:       (K, M) fitted RDR variances (log-space if log_rdr).
        baf_means:      (K, M) fitted BAF means (one mixture component).
        baf_taus:       (M,) fitted BAF dispersion per sample.
        tumor_samples:  list of sample names.
        out_file:       output PDF path (used only when pdf is None).
        log_rdr:        if True, plot log(RDR) and overlay Gaussian in log-space.
        pdf:            open PdfPages to append pages to; if None, opens out_file.
    """
    import warnings
    from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec

    K = rdr_means.shape[0]
    M = rdr_means.shape[1]
    N_total = len(cluster_labels)
    if cluster_ids is None:
        cluster_ids = np.arange(K)

    _close_pdf = pdf is None
    if _close_pdf:
        pdf = PdfPages(out_file)
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore", message="divide by zero", category=RuntimeWarning
        )
        for ki, k in enumerate(cluster_ids):
            mask = cluster_labels == k
            n_bins = int(mask.sum())
            if n_bins == 0:
                continue

            is_multimodal = False
            for m in range(M):
                baf_obs_check = X_bafs[mask, m]
                if log_rdr:
                    rdr_obs_check = np.log(np.clip(X_rdrs[mask, m], 1e-6, None))
                else:
                    rdr_obs_check = X_rdrs[mask, m]
                if _is_multimodal(baf_obs_check) or _is_multimodal(rdr_obs_check):
                    is_multimodal = True
                    break

            title = f"Cluster {k}  (n={n_bins}/{N_total})"
            title_color = "red" if is_multimodal else "black"

            fig = plt.figure(figsize=(7 * M, 12))
            fig.suptitle(title, fontsize=14, y=0.99, color=title_color)
            outer = GridSpec(
                2, M, figure=fig, height_ratios=[3, 2], wspace=0.35, hspace=0.35
            )

            for m in range(M):
                baf_obs = X_bafs[mask, m]
                if log_rdr:
                    rdr_obs = np.log(np.clip(X_rdrs[mask, m], 1e-6, None))
                    mu_k = float(rdr_means[ki, m])
                    var_k = float(rdr_vars[ki, m])
                    ylab = "log(RDR)"
                else:
                    rdr_obs = X_rdrs[mask, m]
                    mu_k = float(rdr_means[ki, m])
                    var_k = float(rdr_vars[ki, m])
                    ylab = "RDR"
                sigma_k = np.sqrt(max(var_k, 1e-12))

                p_k = float(baf_means[ki, m])
                tau_m = float(baf_taus[m])
                a_param = p_k * tau_m
                b_param = (1.0 - p_k) * tau_m

                inner = GridSpecFromSubplotSpec(
                    2,
                    2,
                    subplot_spec=outer[0, m],
                    height_ratios=[1, 4],
                    width_ratios=[4, 1],
                    hspace=0.05,
                    wspace=0.05,
                )
                ax_main = fig.add_subplot(inner[1, 0])
                ax_top = fig.add_subplot(inner[0, 0], sharex=ax_main)
                ax_right = fig.add_subplot(inner[1, 1], sharey=ax_main)
                fig.add_subplot(inner[0, 1]).axis("off")

                ax_main.scatter(
                    baf_obs, rdr_obs, s=1, alpha=0.15, color="0.3", rasterized=True
                )
                if a_param > 0 and b_param > 0:
                    xg = np.linspace(max(0.001, p_k - 0.3), min(0.999, p_k + 0.3), 150)
                    yg = np.linspace(mu_k - 4 * sigma_k, mu_k + 4 * sigma_k, 150)
                    Xg, Yg = np.meshgrid(xg, yg)
                    log_joint = np.log(
                        np.maximum(beta_dist.pdf(Xg, a_param, b_param), 1e-300)
                    ) + norm.logpdf(Yg, mu_k, sigma_k)
                    ax_main.contour(
                        Xg,
                        Yg,
                        np.exp(log_joint),
                        levels=6,
                        colors="red",
                        linewidths=0.8,
                        alpha=0.7,
                    )
                ax_main.set_xlabel("BAF", fontsize=9)
                ax_main.set_ylabel(ylab, fontsize=9)

                ax_top.hist(
                    baf_obs, bins=80, density=True, alpha=0.6, color="steelblue"
                )
                if a_param > 0 and b_param > 0:
                    x_baf = np.linspace(0.001, 0.999, 300)
                    ax_top.plot(
                        x_baf,
                        beta_dist.pdf(x_baf, a_param, b_param),
                        "r-",
                        lw=1.2,
                        label=f"p={p_k:.3f} tau={tau_m:.0f}",
                    )
                    ax_top.legend(fontsize=7, loc="upper right")
                ax_top.tick_params(labelbottom=False)
                ax_top.set_title(tumor_samples[m], fontsize=10)

                ax_right.hist(
                    rdr_obs,
                    bins=80,
                    density=True,
                    alpha=0.6,
                    color="salmon",
                    orientation="horizontal",
                )
                y_rdr = np.linspace(mu_k - 4 * sigma_k, mu_k + 4 * sigma_k, 300)
                ax_right.plot(
                    norm.pdf(y_rdr, mu_k, sigma_k),
                    y_rdr,
                    "r-",
                    lw=1.2,
                    label=f"mu={mu_k:.3f}\nvar={var_k:.4f}",
                )
                ax_right.legend(fontsize=7, loc="upper right")
                ax_right.tick_params(labelleft=False)

                inner_qq = GridSpecFromSubplotSpec(
                    1,
                    2,
                    subplot_spec=outer[1, m],
                    wspace=0.35,
                )
                n_obs = int(mask.sum())
                theoretical_q = np.linspace(1 / (n_obs + 1), n_obs / (n_obs + 1), n_obs)

                ax_baf_qq = fig.add_subplot(inner_qq[0, 0])
                if a_param > 0 and b_param > 0:
                    baf_sorted = np.sort(baf_obs)
                    baf_theo = beta_dist.ppf(theoretical_q, a_param, b_param)
                    ax_baf_qq.scatter(
                        baf_theo, baf_sorted, s=1, alpha=0.3, color="steelblue"
                    )
                    qq_lo = min(baf_theo.min(), baf_sorted.min())
                    qq_hi = max(baf_theo.max(), baf_sorted.max())
                    ax_baf_qq.plot([qq_lo, qq_hi], [qq_lo, qq_hi], "r-", lw=1)
                    ss_res = np.sum((baf_sorted - baf_theo) ** 2)
                    ss_tot = np.sum((baf_sorted - baf_sorted.mean()) ** 2)
                    r2_baf = 1 - ss_res / ss_tot if ss_tot > 0 else float("nan")
                    ax_baf_qq.text(
                        0.05,
                        0.95,
                        f"$R^2$={r2_baf:.4f}",
                        transform=ax_baf_qq.transAxes,
                        fontsize=9,
                        va="top",
                        ha="left",
                    )
                ax_baf_qq.set_xlabel("Theoretical (Beta)", fontsize=9)
                ax_baf_qq.set_ylabel("Observed", fontsize=9)
                ax_baf_qq.set_title(f"BAF QQ", fontsize=10)

                ax_rdr_qq = fig.add_subplot(inner_qq[0, 1])
                rdr_sorted = np.sort(rdr_obs)
                rdr_theo = norm.ppf(theoretical_q, mu_k, sigma_k)
                ax_rdr_qq.scatter(rdr_theo, rdr_sorted, s=1, alpha=0.3, color="salmon")
                qq_lo = min(rdr_theo.min(), rdr_sorted.min())
                qq_hi = max(rdr_theo.max(), rdr_sorted.max())
                ax_rdr_qq.plot([qq_lo, qq_hi], [qq_lo, qq_hi], "r-", lw=1)
                ss_res = np.sum((rdr_sorted - rdr_theo) ** 2)
                ss_tot = np.sum((rdr_sorted - rdr_sorted.mean()) ** 2)
                r2_rdr = 1 - ss_res / ss_tot if ss_tot > 0 else float("nan")
                ax_rdr_qq.text(
                    0.05,
                    0.95,
                    f"$R^2$={r2_rdr:.4f}",
                    transform=ax_rdr_qq.transAxes,
                    fontsize=9,
                    va="top",
                    ha="left",
                )
                ax_rdr_qq.set_xlabel("Theoretical (Gaussian)", fontsize=9)
                ax_rdr_qq.set_ylabel("Observed", fontsize=9)
                ax_rdr_qq.set_title(f"{ylab} QQ", fontsize=10)

            fig.subplots_adjust(top=0.95)
            pdf.savefig(fig)
            plt.close(fig)

    if _close_pdf:
        pdf.close()
        logging.info(f"cluster diagnostics saved to {out_file}")


##################################################
def plot_rdr_baf(
    samples: list,
    bin_info: pd.DataFrame,
    baf_mat: np.ndarray,
    rdr_mat: np.ndarray,
    genome_file: str,
    cluster_labels=None,
    expected_rdrs=None,
    expected_bafs=None,
    unique_labels=None,
    label_clone=False,
    xlab="BAF",
    ylab="RDR",
    out_dir=None,
    out_prefix="",
    dpi=300,
    transparent=False,
    row_width=20,
    row_height=4,
    maxlim_rdr=100,
    rasterized=True,
    rdr_means=None,
    rdr_vars=None,
    baf_taus=None,
    log_rdr=False,
):
    """
    Plot BAF-RDR scatter 1D and 2D.
    baf_mat and rdr_mat have shape (bins, samples).
    Saves all 1D plots to one PDF and all 2D plots to another PDF,
    one page per sample in the order given by samples.
    """
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42
    plt.rcParams["svg.fonttype"] = "none"

    bin_info = bin_info.copy(deep=True)

    chrom_sizes = read_genome_sizes(genome_file)
    palette = None
    _lbl_to_idx = None
    if cluster_labels is not None:
        palette = set_palette(num_colors=len(np.unique(cluster_labels)))
        if unique_labels is not None:
            _lbl_to_idx = np.empty(int(cluster_labels.max()) + 1, dtype=int)
            for i, lbl in enumerate(unique_labels):
                _lbl_to_idx[lbl] = i

    pdf_1d = PdfPages(os.path.join(out_dir, f"{out_prefix}1D.pdf"))
    pdf_2d = PdfPages(os.path.join(out_dir, f"{out_prefix}2D.pdf"))

    for si, sample in enumerate(samples):
        logging.info(f"plot {sample}")
        bafs = baf_mat[:, si]
        rdrs = rdr_mat[:, si]
        exp_bafs = expected_bafs[:, si] if expected_bafs is not None else None
        exp_rdrs = expected_rdrs[:, si] if expected_rdrs is not None else None

        lim_baf = (0, 1) if np.max(bafs) > 0.5 else (0, 0.55)
        max_rdr = np.round(np.max(rdrs)).astype(int)
        if max_rdr > maxlim_rdr:
            num_exceeded = np.sum(rdrs >= maxlim_rdr)
            logging.warning(
                f"there are {num_exceeded} bins having RDR exceed maxlim_rdr={maxlim_rdr}"
            )
        lim_rdr = (0, min(max(3, max_rdr), maxlim_rdr))

        # 2D plot
        fig_2d, g0_colors = plot_2d(
            sample,
            bin_info,
            bafs,
            rdrs,
            exp_xvals=exp_bafs,
            exp_yvals=exp_rdrs,
            exp_labels=unique_labels,
            alphas=None,
            hue=cluster_labels,
            palette=palette,
            label_clone=label_clone,
            xlab=xlab,
            ylab=ylab,
            xlim=lim_baf,
            ylim=lim_rdr,
            title=f"sample={sample}",
            dpi=dpi,
            transparent=transparent,
            rasterized=rasterized,
        )
        pdf_2d.savefig(fig_2d, dpi=dpi, bbox_inches="tight")
        plt.close(fig_2d)

        # 1D plot
        fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(row_width, row_height))
        for i, [ctype, ylim, vals, exp_vals] in enumerate(
            [[ylab, lim_rdr, rdrs, exp_rdrs], [xlab, lim_baf, bafs, exp_bafs]]
        ):
            exp_colname = None
            if cluster_labels is not None:
                exp_colname = f"exp-{ctype}"
                bin_info[exp_colname] = exp_vals[_lbl_to_idx[cluster_labels]]
            plot_1d(
                axes[i],
                sample,
                bin_info,
                vals,
                None,
                chrom_sizes,
                exp_colname=exp_colname,
                val_type=ctype,
                colors=g0_colors,
                hue=cluster_labels if i == 0 else None,
                palette=palette if i == 0 else None,
                ylim=ylim,
                ylab=ctype,
                plot_chrname=True,
                ignore_gap=False,
                rasterized=rasterized,
            )
        fig.suptitle(f"sample={sample}")
        axes[0].grid(True, axis="y")
        axes[0].grid(False, axis="x")
        axes[1].grid(True, axis="y")
        axes[1].grid(False, axis="x")
        if cluster_labels is not None:
            ncol = max(1, int(np.ceil(len(unique_labels) / 10)))
            sns.move_legend(
                axes[0],
                "center left",
                bbox_to_anchor=(1.01, 0.5),
                frameon=False,
                title=None,
                ncol=ncol,
                fontsize=8,
                markerscale=6.0,
            )
        fig.subplots_adjust(right=0.82)
        fig.tight_layout(rect=[0, 0, 0.82, 1])
        pdf_1d.savefig(fig, dpi=dpi, bbox_inches="tight")
        plt.close(fig)

    pdf_1d.close()
    if rdr_means is not None:
        plot_clusters(
            cluster_labels,
            rdr_mat,
            baf_mat,
            rdr_means,
            rdr_vars,
            expected_bafs,
            baf_taus,
            samples,
            cluster_ids=unique_labels,
            log_rdr=log_rdr,
            pdf=pdf_2d,
        )
    pdf_2d.close()
    return
