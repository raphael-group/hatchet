import logging

import numpy as np
import pandas as pd
from scipy.special import betaln, polygamma
from scipy.optimize import minimize_scalar
from scipy.signal import find_peaks
from scipy.stats import gaussian_kde
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


##################################################
def mle_BB_dispersion(
    a_counts: np.ndarray, b_counts: np.ndarray, p=0.5, min_tau=50, max_tau=100
):
    """MLE of Beta-Binomial dispersion tau via Brent's method in log-space.

    Optimises tau over [min_tau, max_tau] using scipy minimize_scalar.

    Args:
        a_counts, b_counts: (N,) arrays of A- and B-allele counts.
        p:        BAF mean — scalar (e.g. 0.5 for balanced bins) or (N,) array
                  of per-bin means (e.g. from cluster assignments).
        min_tau, max_tau: Search bounds for tau.

    Returns:
        tau: MLE dispersion estimate.
    """

    def neg_loglik_logw(logw):
        w = np.exp(logw)
        a0 = w * p
        b0 = w * (1 - p)
        a1 = a_counts + a0
        b1 = b_counts + b0
        ll = np.sum(betaln(a1, b1) - betaln(a0, b0))
        return -ll

    res = minimize_scalar(
        neg_loglik_logw,
        bounds=(np.log(min_tau), np.log(max_tau)),
    )
    tau = np.exp(res.x)
    return tau


def estimate_BB_dispersion_normal(
    X_alphas_normal: np.ndarray,
    X_betas_normal: np.ndarray,
    M: int,
    min_tau=50,
    max_tau=500,
):
    """Estimate BB dispersion tau from the normal sample, shared by all tumor samples.

    The normal sample is diploid (BAF ≈ 0.5 genome-wide), so it provides a
    clean estimate of the sequencing/technical dispersion without copy-number
    confounding.  A single tau is estimated from the normal and broadcast to
    all M tumor samples.

    Args:
        X_alphas_normal: (N,) A-allele counts for the normal sample.
        X_betas_normal:  (N,) B-allele counts for the normal sample.
        M:               Number of tumor samples (for output shape).
        min_tau, max_tau: Bounds passed to mle_BB_dispersion.

    Returns:
        bb_taus: (M,) float32 array with the same tau for every sample.
    """
    logging.info("estimate BB dispersion from normal sample (shared across tumor samples)")
    logging.info(f"tau bound=[{min_tau},{max_tau}]")
    tau = mle_BB_dispersion(
        X_alphas_normal, X_betas_normal, p=0.5, min_tau=min_tau, max_tau=max_tau
    )
    bb_taus = np.full(M, tau, dtype=np.float32)
    return bb_taus


def estimate_BB_dispersion_segment(
    X_alphas: np.ndarray,
    X_betas: np.ndarray,
    X_bafs: np.ndarray,
    X_lengths: np.ndarray,
    M: int,
    min_tau=50,
    max_tau=500,
):
    """Estimate BB dispersion tau from the most balanced segment.

    When no normal sample is available, picks the segment whose bins have
    the smallest mean |BAF - 0.5| across samples (most likely diploid),
    and estimates tau from that segment's bins.  A single tau is estimated
    per sample.

    Args:
        X_alphas, X_betas: (N, M) allele count arrays (tumor only).
        X_bafs:            (N, M) observed BAF values (tumor only).
        X_lengths:         (S,) number of bins per segment.
        M:                 Number of tumor samples.
        min_tau, max_tau:  Bounds passed to mle_BB_dispersion.

    Returns:
        bb_taus: (M,) float32 array of per-sample tau estimates.
    """
    # TODO: support panel-of-normals (PON) file as an alternative source
    logging.info("estimate BB dispersion from most balanced segment (no normal sample)")
    logging.info(f"tau bound=[{min_tau},{max_tau}]")

    seg_starts = np.concatenate(([0], np.cumsum(X_lengths[:-1])))
    best_seg = -1
    best_dev = np.inf
    for s, (start, length) in enumerate(zip(seg_starts, X_lengths)):
        seg_bafs = X_bafs[start : start + length]
        mean_dev = np.mean(np.abs(seg_bafs - 0.5))
        if mean_dev < best_dev:
            best_dev = mean_dev
            best_seg = s

    start = int(seg_starts[best_seg])
    length = int(X_lengths[best_seg])
    logging.info(
        f"selected segment {best_seg} ({length} bins, mean |BAF-0.5|={best_dev:.4f})"
    )

    bb_taus = np.zeros(M, dtype=np.float32)
    for si in range(M):
        bb_taus[si] = mle_BB_dispersion(
            X_alphas[start : start + length, si],
            X_betas[start : start + length, si],
            min_tau=min_tau,
            max_tau=max_tau,
        )
    return bb_taus


def estimate_BB_dispersion_balanced(
    X_alphas: np.ndarray,
    X_betas: np.ndarray,
    X_bafs: np.ndarray,
    M: int,
    min_tau=50,
    max_tau=100,
    bb_quantile=0.2,
):
    """Estimate per-sample BB dispersion tau from near-balanced bins.

    Selects the `bb_quantile` fraction of bins with the smallest mean
    |BAF - 0.5| deviation (i.e., most balanced) and fits tau via MLE
    independently for each tumor sample.

    Args:
        X_alphas, X_betas: (N, M) allele count arrays.
        X_bafs:            (N, M) observed BAF values.
        M:                 Number of tumor samples.
        min_tau, max_tau:  Bounds passed to mle_BB_dispersion.
        bb_quantile:       Fraction of balanced bins to use (default 0.2).

    Returns:
        bb_taus: (M,) float32 array of per-sample tau estimates.
    """
    logging.info("initialize BB dispersion parameter via MLE on balanced bins")
    logging.info(
        f"tau bound=[{min_tau},{max_tau}], balanced quantile={bb_quantile:.3%}"
    )
    mean_baf_dev = np.mean(np.abs(X_bafs - 0.5), axis=1)
    balanced_idx = np.where(mean_baf_dev <= np.quantile(mean_baf_dev, bb_quantile))[0]
    bb_taus = np.zeros(M, dtype=np.float32)
    for si in range(M):
        bb_taus[si] = mle_BB_dispersion(
            X_alphas[balanced_idx][:, si],
            X_betas[balanced_idx][:, si],
            min_tau=min_tau,
            max_tau=max_tau,
        )
    return bb_taus


##################################################
def estimate_rdr_vars(
    X_rdrs: np.ndarray,
    X_lengths: np.ndarray,
    min_var: float = 1e-3,
) -> np.ndarray:
    """Estimate RDR noise variance globally using first-difference MAD.

    First differences between adjacent bins within a segment are used to estimate
    noise variance robustly (immune to CN-state transitions). MAD * 1.4826 gives a
    consistent sigma estimate under normality; dividing by 2 corrects for the fact
    that Var(X_{i+1} - X_i) = 2 * sigma^2 for i.i.d. noise.

    Segment boundaries are masked out so cross-region diffs are excluded.

    Args:
        X_rdrs:    (N, M) RDR array (may be log-transformed).
        X_lengths: (S,) number of bins per segment.
        min_var:   Minimum variance floor (default 1e-3).

    Returns:
        rdr_vars0: (1, M) global variance estimate.
    """
    N, M = X_rdrs.shape
    seg_starts = np.cumsum(np.concatenate(([0], X_lengths[:-1])))

    # Global estimate: first-diff MAD across all within-segment adjacent pairs
    boundary_mask = np.zeros(N - 1, dtype=bool)
    for s in seg_starts[1:]:
        boundary_mask[s - 1] = True
    diff_rdr = np.diff(X_rdrs, axis=0)  # (N-1, M)
    diff_valid = diff_rdr[~boundary_mask]  # mask out cross-segment diffs
    mad = np.median(np.abs(diff_valid - np.median(diff_valid, axis=0)), axis=0)
    raw_global_var = (mad * 1.4826) ** 2 / 2.0
    logging.info(
        f"RDR global var (pre-clip): {np.array2string(raw_global_var, precision=4)}"
    )
    global_var = np.maximum(raw_global_var, min_var)  # (M,)
    return global_var[None, :]  # (1, M)


def compute_baf_se(k_labels, k_betas_phased, X_totals, k_baf_means, k_baf_taus, k_cids):
    """Observed Fisher information SE for per-cluster BAF means.

    Args:
        k_labels:        (N,) cluster assignment per bin.
        k_betas_phased:  (N, M) phased B-allele counts.
        X_totals:        (N, M) total allele counts.
        k_baf_means:     (K, M) fitted BAF means per cluster.
        k_baf_taus:      (M,)   fitted BAF dispersions per sample.
        k_cids:          (K,)   ordered active cluster IDs.

    Returns:
        baf_ses: (K, M) standard errors of BAF means.
    """
    n_clusters = len(k_cids)
    n_samples = k_baf_means.shape[1]
    baf_ses = np.full((n_clusters, n_samples), np.nan)
    for ci, c in enumerate(k_cids):
        mask = k_labels == c
        for m in range(n_samples):
            p = k_baf_means[ci, m]
            tau = k_baf_taus[m]
            a = tau * p
            b = tau * (1.0 - p)
            beta_bins = k_betas_phased[mask, m]
            alpha_bins = X_totals[mask, m] - beta_bins
            fisher = tau**2 * np.sum(
                polygamma(1, a)
                + polygamma(1, b)
                - polygamma(1, beta_bins + a)
                - polygamma(1, alpha_bins + b)
            )
            baf_ses[ci, m] = 1.0 / np.sqrt(fisher) if fisher > 0 else np.inf
    return baf_ses


def compute_rdr_se(k_labels, k_rdr_vars, k_cids):
    """Standard error of cluster-level RDR means.

    For the Gaussian emission model, SE = sqrt(var / n_bins).

    Args:
        k_labels:    (N,) cluster assignment per bin.
        k_rdr_vars:  (K, M) fitted RDR variances per cluster per sample.
        k_cids:      (K,) ordered active cluster IDs.

    Returns:
        rdr_ses: (K, M) standard errors of RDR means.
    """
    n_clusters, n_samples = k_rdr_vars.shape
    rdr_ses = np.full((n_clusters, n_samples), np.nan)
    for ci, c in enumerate(k_cids):
        n_bins = np.sum(k_labels == c)
        if n_bins > 0:
            rdr_ses[ci, :] = np.sqrt(k_rdr_vars[ci, :] / n_bins)
    return rdr_ses


##################################################
def mat2segs(
    bbcs: pd.DataFrame,
    tumor_samples: list,
    baf_means: np.ndarray,
    baf_taus: np.ndarray,
    k_baf_ses: np.ndarray,
    rdr_means: np.ndarray,
    rdr_vars: np.ndarray,
    k_rdr_ses: np.ndarray,
    cluster_ids: np.ndarray,
):
    """Build the SEG summary DataFrame from per-bin BBC data and cluster parameters.

    Groups bins by cluster, then for each (cluster, sample) pair records the
    number of bins, total SNP count, and the fitted emission parameters.

    Args:
        bbcs:         BBC DataFrame with columns CLUSTER, SAMPLE, START, END, #SNPS, ALPHA, BETA, COV, BAF, RD.
        tumor_samples: Ordered list of tumor sample names.
        baf_means:    (K, M) fitted BAF means.
        baf_taus:     (M,)   fitted BAF dispersions (Beta-Binomial tau per sample).
        rdr_means:    (K, M) fitted RDR means.
        rdr_vars:     (K, M) fitted RDR variances.
        cluster_ids:  Ordered array of active cluster IDs (0-indexed).
        baf_ses:      (K, M) BAF standard errors (optional; NaN if not provided).

    Returns:
        segs: DataFrame with columns #ID, SAMPLE, #BINS, #SNPS, LENGTH, ALPHA, BETA, COV, BAF, BAF-se, BAF-tau, RD, RD-var.
              LENGTH is the total base-pair span of bins in the cluster.
              COV is the bin-length-weighted mean depth across bins in the cluster.
    """
    bb_grps = bbcs.groupby(by="CLUSTER", sort=False)
    seg_rows = []
    for l, label in enumerate(cluster_ids):
        bb_grp = bb_grps.get_group(label)
        for s, sample in enumerate(tumor_samples):
            bb_sample = bb_grp.loc[bb_grp["SAMPLE"] == sample, :]
            bin_lengths = (bb_sample["END"] - bb_sample["START"]).to_numpy()
            total_len = bin_lengths.sum()
            cov = (bb_sample["COV"].to_numpy() * bin_lengths).sum() / total_len
            seg_rows.append(
                [
                    label,
                    sample,
                    len(bb_sample),
                    bb_sample["#SNPS"].sum(),
                    total_len,
                    bb_sample["ALPHA"].sum(),
                    bb_sample["BETA"].sum(),
                    cov,
                    baf_means[l, s],
                    k_baf_ses[l, s],
                    baf_taus[s],
                    rdr_means[l, s],
                    k_rdr_ses[l, s],
                    rdr_vars[l, s],
                ]
            )
    segs = pd.DataFrame(
        data=seg_rows,
        columns=[
            "#ID",
            "SAMPLE",
            "#BINS",
            "#SNPS",
            "LENGTH",
            "ALPHA",
            "BETA",
            "COV",
            "BAF",
            "BAF-se",
            "BAF-tau",
            "RD",
            "RD-se",
            "RD-var",
        ],
    )
    return segs


def plot_elbo_traces(traces_per_k: list, out_file: str):
    """Save ELBO traces for all K values into a single PDF, one page per K.

    Args:
        traces_per_k: list of (K, elbo_traces, best_it) tuples in K order.
        out_file: path to the output PDF file.
    """
    with PdfPages(out_file) as pdf:
        for K, elbo_traces, best_it in traces_per_k:
            fig, ax = plt.subplots(figsize=(4, 3))
            for it, trace in elbo_traces.items():
                vals = trace[1:]
                iters = list(range(len(vals)))
                if it == best_it:
                    ax.plot(
                        iters, vals, color="red", linewidth=1.25, zorder=3, label="best"
                    )
                else:
                    ax.plot(iters, vals, color="gray", linewidth=1.0, zorder=2)
            ax.set_xlabel("Iteration")
            ax.set_ylabel("ELBO")
            ax.set_title(f"ELBO trace (K={K})")
            plt.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)


def plot_score(scores_df: pd.DataFrame, score_method: str, out_file: str):
    """Save a two-panel plot of BIC/ICL (top) and log-likelihood (bottom) vs K.

    Each panel shows a line through the best restart per K and a shaded band
    over the min-to-max range across restarts.  A dashed vertical line marks
    the selected K (global minimum of the best-restart score line).
    """
    fig, (ax_score, ax_ll) = plt.subplots(
        2,
        1,
        figsize=(6, 5),
        sharex=True,
        layout="constrained",
        gridspec_kw={"hspace": 0.08},
    )

    Ks = sorted(scores_df["K"].unique())
    best_score, lo_score, hi_score = [], [], []
    best_ll, lo_ll, hi_ll = [], [], []
    for K in Ks:
        rows = scores_df.loc[scores_df["K"] == K]
        sv = rows[score_method].values
        lv = rows["ll"].values
        best_score.append(float(sv.min()))
        lo_score.append(float(sv.min()))
        hi_score.append(float(sv.max()))
        best_ll.append(float(lv.max()))  # higher LL is better
        lo_ll.append(float(lv.min()))
        hi_ll.append(float(lv.max()))

    best_score = np.array(best_score)
    best_ll = np.array(best_ll)
    best_K = Ks[int(np.argmin(best_score))]

    for ax, best, lo, hi, ylabel, color in [
        (ax_score, best_score, lo_score, hi_score, score_method.upper(), "steelblue"),
        (ax_ll, best_ll, lo_ll, hi_ll, "log-likelihood", "tomato"),
    ]:
        ax.fill_between(Ks, lo, hi, alpha=0.20, color=color, label="")
        ax.plot(Ks, best, color=color, marker="o", ms=5, lw=1.5, label="best restart")
        ax.axvline(
            best_K, color=color, linestyle="--", alpha=0.7, label=f"best K={best_K}"
        )
        ax.set_ylabel(ylabel)
        ax.legend(fontsize=7)

    ax_ll.set_xlabel("K")
    ax_ll.set_xticks(Ks)
    ax_score.set_title(f"Model selection ({score_method.upper()})")

    plt.savefig(out_file)
    plt.close()


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


def count_multimodal_clusters(labels, X_rdrs, X_bafs, log_rdr):
    """Count clusters whose marginal RDR or BAF distribution is multimodal.

    Args:
        labels:   (N,) 0-indexed cluster assignments.
        X_rdrs:   (N, M) RDR values (original scale).
        X_bafs:   (N, M) phased BAF values in [0, 1].
        log_rdr:  if True, check multimodality on log(RDR).

    Returns:
        (n_multimodal, multimodal_ids): count and list of multimodal cluster IDs.
    """
    cluster_ids = np.unique(labels)
    M = X_rdrs.shape[1]
    multimodal_ids = []
    for k in cluster_ids:
        mask = labels == k
        for m in range(M):
            baf_obs = X_bafs[mask, m]
            if log_rdr:
                rdr_obs = np.log(np.clip(X_rdrs[mask, m], 1e-6, None))
            else:
                rdr_obs = X_rdrs[mask, m]
            if _is_multimodal(baf_obs) or _is_multimodal(rdr_obs):
                multimodal_ids.append(int(k))
                break
    return len(multimodal_ids), multimodal_ids


def label_multimodal_bins(labels, X_rdrs, X_bafs, log_rdr):
    """Return (N,) bool array: True if bin belongs to a multimodal cluster."""
    _, multimodal_ids = count_multimodal_clusters(labels, X_rdrs, X_bafs, log_rdr)
    return np.isin(labels, multimodal_ids)
