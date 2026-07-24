import os
import logging

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.backends.backend_pdf import PdfPages

from scipy.stats import norm, beta as beta_dist, gaussian_kde
from scipy.signal import find_peaks

from cnplot import (
    plot_scatter_1d,
    plot_scatter_2d,
    annotate_landmarks,
    set_palette,
)
from hatchet.plot.plot_utils import build_genome_axis, use_editable_fonts


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
    palette=None,
    bin_info=None,
    lim_baf: tuple = None,
    lim_rdr: tuple = None,
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
        baf_taus:       (K, M) fitted BAF dispersion per cluster per sample.
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
                tau_m = float(baf_taus[ki, m])
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

                cluster_color = palette[ki] if palette is not None else "0.3"
                ax_main.scatter(
                    baf_obs,
                    rdr_obs,
                    s=4,
                    alpha=0.3,
                    color=cluster_color,
                    rasterized=True,
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
                if lim_baf is not None:
                    ax_main.set_xlim(lim_baf)
                if lim_rdr is not None:
                    ax_main.set_ylim(lim_rdr)

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
                ax_baf_qq.set_title("BAF QQ", fontsize=10)

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
    region_bed: str,
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
    filtered_ids=None,
    balanced_ids=None,
):
    """Plot BAF-RDR scatter 1D and 2D pages, one page each per sample.

    Rendered with cnplot primitives. baf_mat and rdr_mat have shape
    (bins, samples). One page per sample: a 2D BAF-vs-RDR joint plot, then a
    1D genome-wide RDR/BAF panel. When rdr_means is given, per-cluster
    diagnostic pages are appended. Output is one PDF at
    out_dir/<out_prefix rstrip "_">.pdf.

    Args:
        genome_file: Chromosome-sizes file path for the genome axis.
        region_bed: Plotted-regions BED path for the genome axis.
    """
    use_editable_fonts()

    bin_info = bin_info.copy(deep=True)
    genome_axis = build_genome_axis(
        region_bed, genome_file, keep_chroms=bin_info["#CHR"].unique()
    )

    filtered_ids = filtered_ids or set()
    balanced_ids = balanced_ids or set()

    palette = None  # list palette indexed by cluster position, for plot_clusters
    pal_dict = None  # {str(cluster): color} for cnplot hue
    lbl_to_idx = None
    if cluster_labels is not None:
        palette = set_palette(num_colors=len(np.unique(cluster_labels)))
        lbl_to_idx = np.empty(int(cluster_labels.max()) + 1, dtype=int)
        for i, lbl in enumerate(unique_labels):
            lbl_to_idx[lbl] = i
        if filtered_ids:
            gray = (0.75, 0.75, 0.75)
            for i, lbl in enumerate(unique_labels):
                if lbl in filtered_ids:
                    palette[i] = gray
        pal_dict = {str(lbl): palette[i] for i, lbl in enumerate(unique_labels)}

    out_name = out_prefix.rstrip("_") if out_prefix else "plot"
    pdf = PdfPages(os.path.join(out_dir, f"{out_name}.pdf"))

    global_lim_baf = (0, 1) if np.max(baf_mat) > 0.5 else (0, 0.55)
    global_max_rdr = int(np.ceil(np.max(rdr_mat)))
    global_lim_rdr = (0, min(max(2, global_max_rdr), maxlim_rdr))

    for si, sample in enumerate(samples):
        logging.info(f"plot {sample}")
        bafs = baf_mat[:, si]
        rdrs = rdr_mat[:, si]

        lim_baf = (0, 1) if np.max(bafs) > 0.5 else (0, 0.55)
        max_rdr = int(np.ceil(np.max(rdrs)))
        if max_rdr > maxlim_rdr:
            num_exceeded = np.sum(rdrs >= maxlim_rdr)
            logging.warning(
                f"there are {num_exceeded} bins having RDR exceed maxlim_rdr={maxlim_rdr}"
            )
        lim_rdr = (0, min(max(2, max_rdr), maxlim_rdr))

        obs = bin_info[["#CHR", "START", "END"]].copy()
        obs["BAF"] = bafs
        obs["RD"] = rdrs
        exp = None
        hue = None
        if cluster_labels is not None:
            obs["cluster"] = [str(c) for c in cluster_labels]
            hue = "cluster"
            exp = bin_info[["#CHR", "START", "END"]].copy()
            exp[f"exp_RD_{sample}"] = expected_rdrs[lbl_to_idx[cluster_labels], si]
            exp[f"exp_BAF_{sample}"] = expected_bafs[lbl_to_idx[cluster_labels], si]

        # Page 1: 2D plot
        grid = plot_scatter_2d(
            obs,
            xcol="BAF",
            ycol="RD",
            expected_df=None,
            hue=hue,
            palette=pal_dict,
            xlim=lim_baf,
            ylim=lim_rdr,
            xlabel=xlab,
            ylabel=ylab,
            title=f"sample={sample}",
            show_marginals=True,
            show_props=False,
            rasterized=rasterized,
        )
        if cluster_labels is not None:
            landmarks = [
                {
                    "x": float(expected_bafs[i, si]),
                    "y": float(expected_rdrs[i, si]),
                    "label": str(lbl),
                    "clonal": False,
                    "balanced": lbl in balanced_ids,
                }
                for i, lbl in enumerate(unique_labels)
                if lbl not in filtered_ids
            ]
            if landmarks:
                annotate_landmarks(grid.ax_joint, landmarks)
        pdf.savefig(grid.figure, dpi=dpi, bbox_inches="tight", transparent=transparent)
        plt.close(grid.figure)

        # Page 2: 1D plot
        fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(row_width, row_height))
        plot_scatter_1d(
            axes[0],
            obs,
            genome_axis,
            "RD",
            expected_df=exp,
            group=sample,
            hue=hue,
            palette=pal_dict,
            ylim=lim_rdr,
            ylabel=ylab,
            plot_chrname=False,
            show_legend=False,
            rasterized=rasterized,
        )
        plot_scatter_1d(
            axes[1],
            obs,
            genome_axis,
            "BAF",
            expected_df=exp,
            group=sample,
            hue=None,
            ylim=lim_baf,
            ylabel=xlab,
            href=0.5,
            plot_chrname=True,
            rasterized=rasterized,
        )
        fig.suptitle(f"sample={sample}")
        if cluster_labels is not None:
            handles = [
                Line2D(
                    [0],
                    [0],
                    marker="o",
                    linestyle="",
                    markersize=5,
                    color=pal_dict[str(lbl)],
                    label=str(lbl),
                )
                for lbl in unique_labels
            ]
            ncol = max(1, int(np.ceil(len(unique_labels) / 10)))
            axes[0].legend(
                handles=handles,
                loc="center left",
                bbox_to_anchor=(1.01, 0.5),
                frameon=False,
                title=None,
                ncol=ncol,
                fontsize=8,
            )
        fig.subplots_adjust(right=0.82)
        fig.tight_layout(rect=[0, 0, 0.82, 1])
        pdf.savefig(fig, dpi=dpi, bbox_inches="tight")
        plt.close(fig)

    # Cluster-level diagnostic pages
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
            pdf=pdf,
            palette=palette,
            bin_info=bin_info,
            lim_baf=global_lim_baf,
            lim_rdr=global_lim_rdr,
        )
    pdf.close()
    return
