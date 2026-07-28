"""Plotting utilities for HMM initialization diagnostics."""

from __future__ import annotations

import logging
import os

import numpy as np

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from adjustText import adjust_text

from hatchet.utils import sort_chroms

logging.getLogger("adjustText").setLevel(logging.WARNING)


def _clamp_texts(texts, xlim, ylim):
    """Clamp text positions back inside axes limits after adjust_text."""
    for t in texts:
        x, y = t.get_position()
        t.set_position((np.clip(x, xlim[0], xlim[1]), np.clip(y, ylim[0], ylim[1])))


##################################################
def plot_init_sampling_probs(
    X_rdrs: np.ndarray,
    X_bafs: np.ndarray,
    probs_history: list[np.ndarray],
    centroids_history: list[tuple[np.ndarray, np.ndarray]],
    plot_dir: str,
    name: str,
    log_rdr: bool = False,
    bin_info=None,
    chrom_sizes=None,
    selected_bins_history=None,
    candidates_history=None,
    probs_entropy_history=None,
    final_baf_means: np.ndarray | None = None,
    final_rdr_means: np.ndarray | None = None,
):
    """Save a single multi-page PDF of init diagnostics.

    Each page corresponds to one K-seeding step.  Left portion: genome-wide
    1D tracks (sampling-prob bar + RDR/BAF per sample).  Right portion: M 2D
    scatter plots, one per sample, each aligned with its sample's 1D rows.
    An optional final page shows the post-screening cluster means.
    """
    # Always plot in natural RDR space.
    if log_rdr:
        X_rdrs = np.exp(X_rdrs)
    M = X_rdrs.shape[1]
    K_total = len(probs_history) + 1
    rdr_ylim_2d = int(np.round(X_rdrs.max()))
    nbins_baf, nbins_rdr = 200, 200
    baf_edges = np.linspace(0.0, 1.0, nbins_baf + 1)
    rdr_edges = np.linspace(0.0, float(rdr_ylim_2d), nbins_rdr + 1)
    cmap_2d = plt.get_cmap("YlOrRd")
    cmap_tab10 = plt.get_cmap("tab10")

    has_1d = bin_info is not None and chrom_sizes is not None
    bin_colors = None
    if has_1d:
        bin_chrs = bin_info["#CHR"].to_numpy()
        chroms_in_bb = set(bin_chrs)
        sorted_chroms = [
            c for c in sort_chroms(list(chrom_sizes.keys())) if c in chroms_in_bb
        ]
        chrom_offsets = {}
        offset = 0
        for chrom in sorted_chroms:
            chrom_offsets[chrom] = offset
            offset += chrom_sizes[chrom]
        cmap_chrom = plt.get_cmap("tab20")
        chrom_color = {c: cmap_chrom(i % 20) for i, c in enumerate(sorted_chroms)}
        bin_colors = np.array([chrom_color[c] for c in bin_chrs])
        total_genome = offset
        bin_mids = (bin_info["START"].to_numpy() + bin_info["END"].to_numpy()) / 2.0
        abs_mids = np.array(
            [chrom_offsets.get(c, 0) + mid for c, mid in zip(bin_chrs, bin_mids)]
        )
        n_plot_bins = 200
        bin_edges_1d = np.linspace(0, total_genome, n_plot_bins + 1)
        bin_centers = (bin_edges_1d[:-1] + bin_edges_1d[1:]) / 2
        bar_width = (bin_edges_1d[1] - bin_edges_1d[0]) * 0.9
        rdr_ylim_1d = float(np.nanpercentile(X_rdrs, 99)) * 1.1

    out = os.path.join(plot_dir, f"{name}_init.pdf")
    pdf_pages = PdfPages(out)
    prev_summed = None

    for step_idx, (probs, (baf_centroids, rdr_centroids_raw)) in enumerate(
        zip(probs_history, centroids_history)
    ):
        rdr_centroids = np.exp(rdr_centroids_raw) if log_rdr else rdr_centroids_raw
        k_current = baf_centroids.shape[0]
        ent_str = (
            f"  H={probs_entropy_history[step_idx]:.3f}"
            if probs_entropy_history is not None
            else ""
        )

        # --- 2D: compute histograms and vmin/vmax across samples ---
        prob_ts = []
        for m in range(M):
            prob_hist, _, _ = np.histogram2d(
                X_bafs[:, m],
                X_rdrs[:, m],
                bins=[baf_edges, rdr_edges],
                weights=probs,
            )
            prob_ts.append(prob_hist.T)
        populated = np.concatenate(
            [pt[pt > 0] for pt in prob_ts if np.any(pt > 0)], axis=None
        )
        if len(populated) > 0:
            vmin, vmax = populated.min(), populated.max()
            pad = (vmax - vmin) * 0.05
            vmin = max(0.0, vmin - pad)
            vmax = vmax + pad
        else:
            vmin, vmax = 0.0, 1.0

        # --- Layout ---
        # nrows = 1 + 2*M  (row 0: sampling-prob bar; rows 1+2m, 2+2m: RDR/BAF for sample m)
        # ncols = 2         (col 0: 1D tracks; col 1: 2D scatters stacked, each spanning 2 rows)
        # Each 2D scatter for sample m spans rows [1+2m : 3+2m], col 1.
        nrows = 1 + 2 * M
        ncols = 2
        height_ratios = [1.0] + [1.5] * (2 * M)
        width_ratios = [4.0, 1.5]
        fig = plt.figure(figsize=(26, 3 + 4 * M))
        gs = fig.add_gridspec(
            nrows,
            ncols,
            height_ratios=height_ratios,
            width_ratios=width_ratios,
            hspace=0.15,
            wspace=0.05,
        )
        fig.subplots_adjust(left=0.04, right=0.99, top=0.94, bottom=0.03)

        # --- Row 0, col 0: sampling probability bar chart ---
        ax_prob = fig.add_subplot(gs[0, 0])
        if has_1d:
            summed, _ = np.histogram(abs_mids, bins=bin_edges_1d, weights=probs)
            if prev_summed is not None:
                bar_colors = np.where(
                    summed < prev_summed,
                    "orange",
                    np.where(summed > prev_summed, "green", "steelblue"),
                )
            else:
                bar_colors = ["steelblue"] * len(bin_centers)
            prev_summed = summed
            ax_prob.bar(bin_centers, summed, width=bar_width, color=bar_colors)
            for chrom in sorted_chroms:
                ax_prob.axvline(
                    chrom_offsets[chrom], color="gray", linewidth=0.5, linestyle="--"
                )
            ymax = max(float(summed.max()) * 1.1, 0.01)
            for chrom in sorted_chroms:
                mid = chrom_offsets[chrom] + chrom_sizes[chrom] / 2
                ax_prob.text(
                    mid,
                    ymax,
                    chrom.replace("chr", ""),
                    ha="center",
                    va="top",
                    fontsize=7,
                )
            if selected_bins_history is not None:
                for bin_idx in selected_bins_history[:step_idx]:
                    ax_prob.plot(
                        abs_mids[bin_idx],
                        0,
                        marker="|",
                        markersize=10,
                        color="black",
                        markeredgewidth=1.5,
                        clip_on=False,
                        zorder=5,
                    )
                if candidates_history is not None and step_idx < len(
                    candidates_history
                ):
                    cands, best_idx_c = candidates_history[step_idx]
                    for cand in cands:
                        if cand != best_idx_c:
                            ax_prob.plot(
                                abs_mids[cand],
                                0,
                                marker="|",
                                markersize=10,
                                color="steelblue",
                                markeredgewidth=1.5,
                                clip_on=False,
                                zorder=4,
                            )
                if step_idx < len(selected_bins_history):
                    ax_prob.plot(
                        abs_mids[selected_bins_history[step_idx]],
                        0,
                        marker="|",
                        markersize=10,
                        color="red",
                        markeredgewidth=1.5,
                        clip_on=False,
                        zorder=6,
                    )
            ax_prob.set_xlim(0, total_genome)
            ax_prob.set_ylim(0, ymax)
            ax_prob.set_xticks([])
        ax_prob.set_ylabel("sampling prob (sum)", fontsize=8)
        fig.suptitle(f"K step {k_current}/{K_total}{ent_str}", fontsize=10)

        # --- Hard-assign bins to nearest centroid for 1D coloring ---
        if has_1d:
            dists = np.stack(
                [
                    np.mean(
                        np.minimum(
                            (X_bafs - baf_centroids[k]) ** 2,
                            (X_bafs - (1.0 - baf_centroids[k])) ** 2,
                        )
                        + (X_rdrs - rdr_centroids[k]) ** 2,
                        axis=1,
                    )
                    for k in range(k_current)
                ],
                axis=1,
            )  # (N, k_current)
            assignments = np.argmin(dists, axis=1)  # (N,)

        for m in range(M):
            # --- 1D tracks: col 0, rows 1+2m (RDR) and 2+2m (BAF) ---
            if has_1d:
                ax_rdr = fig.add_subplot(gs[1 + 2 * m, 0])
                ax_baf = fig.add_subplot(gs[2 + 2 * m, 0])
                for k in range(k_current):
                    is_new = k == k_current - 1
                    color = "red" if is_new else cmap_tab10(k % 10)
                    lw = 0.8 if is_new else 0.5
                    zorder = 5 if is_new else 3
                    mask = assignments == k
                    if not np.any(mask):
                        continue
                    ax_rdr.scatter(
                        abs_mids[mask],
                        X_rdrs[mask, m],
                        s=2.0,
                        alpha=1.0,
                        color="lightblue",
                        rasterized=True,
                        linewidths=0,
                        zorder=zorder,
                    )
                    ax_baf.scatter(
                        abs_mids[mask],
                        X_bafs[mask, m],
                        s=2.0,
                        alpha=1.0,
                        color="lightblue",
                        rasterized=True,
                        linewidths=0,
                        zorder=zorder,
                    )
                    if is_new:
                        ax_rdr.axhline(
                            float(rdr_centroids[k, m]),
                            color=color,
                            lw=lw,
                            zorder=zorder + 1,
                            label=f"k={k}",
                        )
                        ax_baf.axhline(
                            float(baf_centroids[k, m]),
                            color=color,
                            lw=lw,
                            zorder=zorder + 1,
                        )
                if candidates_history is not None and step_idx < len(
                    candidates_history
                ):
                    cands, best_idx_c = candidates_history[step_idx]
                    for cand in cands:
                        color = "red" if cand == best_idx_c else "steelblue"
                        zorder = 8 if cand == best_idx_c else 6
                        size = 20 if cand == best_idx_c else 12
                        ax_rdr.scatter(
                            abs_mids[cand],
                            X_rdrs[cand, m],
                            s=size,
                            color=color,
                            zorder=zorder,
                            linewidths=0,
                        )
                        ax_baf.scatter(
                            abs_mids[cand],
                            X_bafs[cand, m],
                            s=size,
                            color=color,
                            zorder=zorder,
                            linewidths=0,
                        )
                for chrom in sorted_chroms:
                    ax_rdr.axvline(
                        chrom_offsets[chrom],
                        color="gray",
                        linewidth=0.3,
                        linestyle="--",
                    )
                    ax_baf.axvline(
                        chrom_offsets[chrom],
                        color="gray",
                        linewidth=0.3,
                        linestyle="--",
                    )
                ax_rdr.set_xlim(0, total_genome)
                ax_baf.set_xlim(0, total_genome)
                ax_rdr.set_ylim(0, rdr_ylim_1d)
                ax_baf.set_ylim(-0.05, 1.05)
                ax_rdr.set_xticks([])
                ax_baf.set_xticks([])
                ax_rdr.set_ylabel(f"sample {m + 1}\nRDR", fontsize=8)
                ax_baf.set_ylabel("BAF", fontsize=8)
                if k_current <= 12:
                    ax_rdr.legend(fontsize=6, loc="upper right", ncol=2)

            # --- 2D scatter: col 1, rows 1+2m..2+2m (spans 2 rows) ---
            ax2d = fig.add_subplot(gs[1 + 2 * m : 3 + 2 * m, 1])
            ax2d.scatter(
                X_bafs[:, m],
                X_rdrs[:, m],
                c=bin_colors if bin_colors is not None else "gray",
                s=1.0,
                alpha=0.4,
                rasterized=True,
                zorder=1,
            )
            prob_t = prob_ts[m]
            has_data = prob_t > 0
            norm_grid = np.where(
                has_data, np.clip((prob_t - vmin) / (vmax - vmin), 0.0, 1.0), 0.0
            )
            rgba = cmap_2d(norm_grid)
            rgba[..., 3] = np.where(has_data, 0.3 + 0.7 * norm_grid, 0.0)
            ax2d.imshow(
                rgba,
                extent=[baf_edges[0], baf_edges[-1], rdr_edges[0], rdr_edges[-1]],
                origin="lower",
                aspect="auto",
                interpolation="nearest",
                zorder=5,
            )
            texts = []
            for k in range(k_current):
                color = "red" if k == k_current - 1 else "black"
                ax2d.plot(
                    float(baf_centroids[k, m]),
                    float(rdr_centroids[k, m]),
                    marker="+",
                    markersize=6,
                    markeredgewidth=1.0,
                    color=color,
                    zorder=10,
                )
                texts.append(
                    ax2d.text(
                        float(baf_centroids[k, m]),
                        float(rdr_centroids[k, m]),
                        f" {k}",
                        fontsize=7,
                        color=color,
                        va="bottom",
                        zorder=10,
                    )
                )
            adjust_text(
                texts, ax=ax2d, arrowprops=dict(arrowstyle="-", color="gray", lw=0.5)
            )
            _clamp_texts(texts, (-0.01, 1.01), (0.0, rdr_ylim_2d))
            if candidates_history is not None and step_idx < len(candidates_history):
                cands, best_idx = candidates_history[step_idx]
                cand_texts = []
                for ci, cand in enumerate(cands):
                    if cand == best_idx:
                        continue
                    label = chr(ord("a") + ci)
                    ax2d.scatter(
                        float(X_bafs[cand, m]),
                        float(X_rdrs[cand, m]),
                        marker="o",
                        s=25,
                        color="steelblue",
                        zorder=8,
                        alpha=0.9,
                        linewidths=0,
                    )
                    cand_texts.append(
                        ax2d.text(
                            float(X_bafs[cand, m]),
                            float(X_rdrs[cand, m]),
                            f" {label}",
                            fontsize=6,
                            color="steelblue",
                            va="bottom",
                            zorder=9,
                        )
                    )
                adjust_text(
                    cand_texts,
                    ax=ax2d,
                    arrowprops=dict(arrowstyle="-", color="gray", lw=0.4),
                )
                _clamp_texts(cand_texts, (-0.01, 1.01), (0.0, rdr_ylim_2d))
            ax2d.set_xlim(-0.01, 1.01)
            ax2d.set_ylim(0.0, rdr_ylim_2d)
            ax2d.set_xlabel("BAF", fontsize=8)
            ax2d.set_ylabel("RDR", fontsize=8)
            ax2d.set_box_aspect(1)

        pdf_pages.savefig(fig, dpi=72)
        plt.close(fig)

    # --- Final page: post-screening cluster means (2D only) ---
    if final_baf_means is not None and final_rdr_means is not None:
        rdr_final = np.exp(final_rdr_means) if log_rdr else final_rdr_means
        K_final = final_baf_means.shape[0]
        fig, axes = plt.subplots(1, M, figsize=(5 * M, 5), squeeze=False)
        for m in range(M):
            ax = axes[0, m]
            ax.scatter(
                X_bafs[:, m],
                X_rdrs[:, m],
                c=bin_colors if bin_colors is not None else "gray",
                s=1.0,
                alpha=0.4,
                rasterized=True,
                zorder=1,
            )
            ax.scatter(
                final_baf_means[:, m],
                rdr_final[:, m],
                facecolors="none",
                edgecolors="black",
                s=10,
                linewidths=0.8,
                zorder=9,
            )
            texts = []
            for k in range(K_final):
                ax.plot(
                    float(final_baf_means[k, m]),
                    float(rdr_final[k, m]),
                    marker="+",
                    markersize=6,
                    markeredgewidth=1.0,
                    color="black",
                    zorder=10,
                )
                texts.append(
                    ax.text(
                        float(final_baf_means[k, m]),
                        float(rdr_final[k, m]),
                        f" {k}",
                        fontsize=7,
                        color="black",
                        va="bottom",
                        zorder=10,
                    )
                )
            adjust_text(
                texts, ax=ax, arrowprops=dict(arrowstyle="-", color="gray", lw=0.5)
            )
            _clamp_texts(texts, (-0.01, 1.01), (0.0, rdr_ylim_2d))
            ax.set_xlim(-0.01, 1.01)
            ax.set_ylim(0.0, rdr_ylim_2d)
            ax.set_xlabel("BAF")
            ax.set_ylabel("log(RDR)" if log_rdr else "RDR")
            ax.set_title(f"sample {m + 1} — final init means")
        fig.suptitle("Final init means (post-screening)")
        plt.tight_layout()
        pdf_pages.savefig(fig, dpi=100)
        plt.close(fig)

    pdf_pages.close()
    logging.info(f"saved init diagnostics to {out}")


##################################################
def plot_2d_inits(
    X_rdrs: np.ndarray,
    X_bafs: np.ndarray,
    params_dict: dict,
    M: int,
    K: int,
    out_file: str,
    baf_taus: np.ndarray | None = None,
    log_rdr: bool = False,
    bbs=None,
    chrom_sizes: dict | None = None,
    init_method: str = "cna_plus_plus",
    sample_names: list | None = None,
):
    """Save a multi-page PDF showing 2-D RDR-vs-mhBAF scatter for every restart.

    Each page corresponds to one restart, sorted by log-likelihood (best first).
    The title of the best restart is highlighted in red.  Cluster centroids are
    annotated with numbered labels.  When *baf_taus* is provided, 1-sigma and
    2-sigma ellipses are drawn around each centroid to visualize emission spread.

    Args:
        X_rdrs:     (N, M) RDR observations (numpy).
        X_bafs:     (N, M) BAF observations (numpy).
        params_dict: Mapping restart_index → [baf_means, rdr_means, rdr_vars, potential].
        M:          Number of tumor samples (one subplot column per sample).
        K:          Number of clusters.
        out_file:   Output PDF file path.
        baf_taus:   (M,) Beta-Binomial dispersion per sample (optional; unused, kept for API compatibility).
    """
    ylim = (
        np.ceil(X_rdrs.max() * 1.1 * 10) / 10
    )  # 10% headroom, round up to nearest 0.1
    pdf_pages = PdfPages(out_file)

    if bbs is not None and chrom_sizes is not None:
        bin_chrs = bbs["#CHR"].to_numpy()
        sorted_chroms = sort_chroms([c for c in chrom_sizes if c in set(bin_chrs)])
        cmap_chrom = plt.get_cmap("tab20")
        chrom_color = {c: cmap_chrom(i % 20) for i, c in enumerate(sorted_chroms)}
        _colors = np.array([chrom_color[c] for c in bin_chrs])
    else:
        _colors = "gray"

    # --- Page 0: plain scatter, no cluster annotations ---
    fig, axes = plt.subplots(1, M, figsize=(5 * M, 5), squeeze=False)
    axes = axes.ravel()
    for m, ax in enumerate(axes):
        ax.scatter(
            X_bafs[:, m],
            X_rdrs[:, m],
            c=_colors,
            s=2.0,
            alpha=0.6,
            rasterized=True,
        )
        ax.set_rasterization_zorder(0)
        ax.vlines(
            0.5,
            ymin=0,
            ymax=1,
            transform=ax.get_xaxis_transform(),
            linewidth=0.5,
            colors="k",
        )
        ax.set_aspect("auto")
        ax.set_title(
            sample_names[m] if sample_names is not None else f"tumor sample {m + 1}"
        )
        ax.set_xlim(-0.01, 1.01)
        ax.set_ylim(0.0, ylim)
    fig.supxlabel("BAF")
    fig.supylabel("RDR")
    fig.suptitle(f"{init_method} initialization K={K}")
    plt.tight_layout()
    pdf_pages.savefig(fig, dpi=100)
    plt.close()

    # --- Pages 1..R: one page per restart, sorted by ll (best first) ---
    # Use restart_idx as the lambda param to avoid shadowing the cluster-index k below.
    sorted_keys = sorted(
        params_dict.keys(),
        key=lambda restart_idx: params_dict[restart_idx][-1],
        reverse=True,
    )
    for it in sorted_keys:
        [baf_means, rdr_means, rdr_vars, potential] = params_dict[it]
        rdr_means_plot = np.exp(rdr_means) if log_rdr else rdr_means
        fig, axes = plt.subplots(1, M, figsize=(5 * M, 5), squeeze=False)
        axes = axes.ravel()
        for m, ax in enumerate(axes):
            ax.scatter(
                X_bafs[:, m],
                X_rdrs[:, m],
                c=_colors,
                s=2.0,
                alpha=0.6,
                label="BB",
                rasterized=True,
            )
            ax.set_rasterization_zorder(0)
            ax.scatter(
                baf_means[:, m],
                rdr_means_plot[:, m],
                facecolors="none",
                edgecolors="black",
                s=10,
                linewidths=0.8,
                zorder=9,
            )
            texts = []
            for k in range(K):
                texts.append(
                    ax.text(
                        float(baf_means[k, m]),
                        float(rdr_means_plot[k, m]),
                        str(k),
                        fontsize=8,
                        color="black",
                        va="center",
                        ha="center",
                        zorder=10,
                        fontweight="bold",
                    )
                )
            adjust_text(
                texts, ax=ax, arrowprops=dict(arrowstyle="-", color="gray", lw=0.5)
            )
            _clamp_texts(texts, (-0.01, 1.01), (0.0, ylim))
            ax.vlines(
                0.5,
                ymin=0,
                ymax=1,
                transform=ax.get_xaxis_transform(),
                linewidth=0.5,
                colors="k",
            )
            ax.set_aspect("auto")
            ax.set_title(
                sample_names[m] if sample_names is not None else f"tumor sample {m + 1}"
            )
            ax.set_xlim(-0.01, 1.01)
            ax.set_ylim(0.0, ylim)
            ax.legend()
        fig.supxlabel("mhBAF")
        fig.supylabel("RDR")
        fig.suptitle(f"{init_method} initialization K={K}, restart={it}")
        plt.tight_layout()
        pdf_pages.savefig(fig, dpi=100)
        plt.close()
    pdf_pages.close()
