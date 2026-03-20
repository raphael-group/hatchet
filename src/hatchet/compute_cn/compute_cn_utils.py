import os
import logging
import pandas as pd
import numpy as np

from hatchet.utils import read_region_bed, build_seg_from_bbc
from hatchet.plot import plot_cn as _plot_cn


def build_cluster_data(segs):
    """Pivot cluster-level SEG into (cluster x sample) DataFrames.

    Returns a dict with keys ``rdr``, ``baf``, ``rdr_se``, ``baf_se``,
    ``nbins``, ``weights`` — the same interface as ``build_segment_data``.
    """
    segs_sorted = segs.sort_values(["#ID", "SAMPLE"])
    rdr = segs_sorted.pivot(index="#ID", columns="SAMPLE", values="RD")
    baf = segs_sorted.pivot(index="#ID", columns="SAMPLE", values="BAF")
    rdr_se = segs_sorted.pivot(index="#ID", columns="SAMPLE", values="RD-se")
    baf_se = segs_sorted.pivot(index="#ID", columns="SAMPLE", values="BAF-se")
    nbins = segs_sorted.pivot(index="#ID", columns="SAMPLE", values="#BINS")
    first_sample = sorted(segs["SAMPLE"].unique())[0]
    lengths = (
        segs.loc[segs["SAMPLE"] == first_sample]
        .set_index("#ID")["LENGTH"]
        .sort_index()
    )
    weights = 100 * lengths / lengths.sum()
    return {
        "rdr": rdr,
        "baf": baf,
        "rdr_se": rdr_se,
        "baf_se": baf_se,
        "nbins": nbins,
        "weights": weights,
    }


def build_segment_data(bbcs, segs):
    """Build genomic-segment-level DataFrames from bin-level BBC and cluster-level SEG.

    A genomic segment is a maximal contiguous run of bins with the same CLUSTER
    on the same chromosome.  Each segment inherits its RD, BAF, and associated
    standard-error columns from the cluster-level ``segs`` DataFrame (no
    recomputation from bin values).

    Parameters
    ----------
    bbcs : pd.DataFrame
        Bin-level BBC DataFrame with columns ``#CHR``, ``START``, ``END``,
        ``SAMPLE``, ``CLUSTER``.
    segs : pd.DataFrame
        Cluster-level SEG DataFrame with columns ``#ID``, ``SAMPLE``, ``RD``,
        ``BAF``, ``RD-se``, ``BAF-se``, ``#BINS``, ``LENGTH``.

    Returns
    -------
    dict
        Keys: ``rdr``, ``baf``, ``rdr_se``, ``baf_se``, ``nbins`` (DataFrames
        of shape ``(num_segments, num_samples)``), ``weights`` (Series of
        length ``num_segments``), ``seg_to_cluster`` (Series mapping segment
        index to cluster ID).
    """
    bbcs = bbcs.sort_values(["#CHR", "START", "END", "SAMPLE"]).reset_index(drop=True)
    samples_sorted = sorted(bbcs["SAMPLE"].unique())
    first_sample = samples_sorted[0]

    # Identify segment boundaries using only the first sample's rows
    mask = bbcs["SAMPLE"] == first_sample
    first_df = bbcs.loc[mask].reset_index(drop=True)

    seg_boundary = (
        (first_df["CLUSTER"] != first_df["CLUSTER"].shift())
        | (first_df["#CHR"] != first_df["#CHR"].shift())
    )
    seg_ids_first = seg_boundary.cumsum() - 1  # 0-based segment IDs

    # Map segment IDs back to all rows (same positional order per sample)
    n_bins_per_sample = len(first_df)
    bbcs["_seg_id"] = np.tile(seg_ids_first.values, len(samples_sorted))

    # Aggregate per (segment, sample)
    agg = (
        bbcs.groupby(["_seg_id", "SAMPLE"])
        .agg(
            CHR=("#CHR", "first"),
            START=("START", "min"),
            END=("END", "max"),
            CLUSTER=("CLUSTER", "first"),
            NBINS=("CLUSTER", "count"),
        )
        .reset_index()
    )
    agg["LENGTH"] = agg["END"] - agg["START"]

    # Join cluster-level statistics from segs
    seg_cols = ["#ID", "SAMPLE", "RD", "BAF", "RD-se", "BAF-se"]
    agg = agg.merge(
        segs[seg_cols],
        left_on=["CLUSTER", "SAMPLE"],
        right_on=["#ID", "SAMPLE"],
        how="left",
    )

    # Pivot into (seg_id x sample) DataFrames
    rdr = agg.pivot(index="_seg_id", columns="SAMPLE", values="RD")
    baf = agg.pivot(index="_seg_id", columns="SAMPLE", values="BAF")
    rdr_se = agg.pivot(index="_seg_id", columns="SAMPLE", values="RD-se")
    baf_se = agg.pivot(index="_seg_id", columns="SAMPLE", values="BAF-se")
    nbins = agg.pivot(index="_seg_id", columns="SAMPLE", values="NBINS")

    # Weights: segment length as percentage of genome (from first sample)
    first_agg = agg.loc[agg["SAMPLE"] == first_sample].set_index("_seg_id")
    seg_lengths = first_agg["LENGTH"]
    weights = 100 * seg_lengths / seg_lengths.sum()

    # Cluster assignment per segment
    seg_to_cluster = first_agg["CLUSTER"]

    bbcs.drop(columns=["_seg_id"], inplace=True)

    n_segs = rdr.shape[0]
    n_samples = rdr.shape[1]
    logging.info(
        f"segment mode: {n_segs} genomic segments x {n_samples} samples "
        f"(from {len(segs['#ID'].unique())} clusters)"
    )

    return {
        "rdr": rdr,
        "baf": baf,
        "rdr_se": rdr_se,
        "baf_se": baf_se,
        "nbins": nbins,
        "weights": weights,
        "seg_to_cluster": seg_to_cluster,
    }


def filtering(
    bbc: pd.DataFrame,
    seg: pd.DataFrame,
    samples: list,
    clusters: list,
    fstd=2.0,
    min_nbins=10,
    ub_nbins=50,
):
    """Filter clusters before the optimization step using variance outlier detection.

    Steps:
        0. Remove any cluster with fewer than ``min_nbins`` bins.
        1. Compute per-sample, per-cluster RD and BAF variance (SCV).
        2. Compute per-sample mean variance (MV) and standard deviation (STDV)
           across clusters that passed step 0.
        3. Mark a cluster as an outlier if, across all samples, its SCV deviates
           from MV by more than ``fstd`` standard deviations, AND the cluster
           has at most ``ub_nbins`` bins.

    Args:
        bbc: Bin-level DataFrame with columns ``SAMPLE``, ``CLUSTER``, ``RD``, ``BAF``.
        seg: Segment-level DataFrame with columns ``SAMPLE``, ``#ID``, ``RD``, ``BAF``,
            ``#BINS``.
        samples: Ordered list of sample identifiers.
        clusters: Ordered list of cluster identifiers.
        fstd: Number of standard deviations used as the outlier threshold.
        min_nbins: Clusters with fewer bins than this are always removed.
        ub_nbins: Outlier detection only applies to clusters with at most this
            many bins (large clusters are kept regardless of variance).

    Returns:
        A tuple ``(good_clusters, bad_clusters)`` where each element is a list
        of cluster IDs.
    """
    logging.info("preprocessing, filtering clusters")

    var_rd_matrix = np.zeros((len(clusters), len(samples)), dtype=np.float64)
    var_baf_matrix = np.zeros((len(clusters), len(samples)), dtype=np.float64)

    for i, cluster in enumerate(clusters):
        for j, sample in enumerate(samples):
            bbc_ = bbc[(bbc["SAMPLE"] == sample) & (bbc["CLUSTER"] == cluster)]
            seg_ = seg[(seg["SAMPLE"] == sample) & (seg["#ID"] == cluster)]
            seg_baf = seg_["BAF"].iloc[0]
            seg_rdr = seg_["RD"].iloc[0]
            var_rd_matrix[i, j] = np.linalg.norm(bbc_["RD"] - seg_rdr, 2) / len(bbc_)
            var_baf_matrix[i, j] = np.linalg.norm(bbc_["BAF"] - seg_baf, 2) / len(bbc_)

    cluster_filtered = np.zeros(len(clusters), dtype=bool)
    for i, cluster in enumerate(clusters):
        cluster_filtered[i] = seg[seg["#ID"] == cluster]["#BINS"].iloc[0] < min_nbins

    mv_rd = np.mean(var_rd_matrix[~cluster_filtered, :], axis=0)
    stdv_rd = np.std(var_rd_matrix[~cluster_filtered, :], axis=0, ddof=1)
    mv_baf = np.mean(var_baf_matrix[~cluster_filtered, :], axis=0)
    stdv_baf = np.std(var_baf_matrix[~cluster_filtered, :], axis=0, ddof=1)
    for j, sample in enumerate(samples):
        lb_rd = mv_rd[j] - fstd * stdv_rd[j]
        ub_rd = mv_rd[j] + fstd * stdv_rd[j]
        lb_baf = mv_baf[j] - fstd * stdv_baf[j]
        ub_baf = mv_baf[j] + fstd * stdv_baf[j]
        logging.debug(
            f"{sample} RD-var bound=({lb_rd:.6f}, {ub_rd:.6f}) BAF-var bound=({lb_baf:.6f}, {ub_baf:.6f})"
        )

    good_clusters = []
    bad_clusters = []
    for i, cluster in enumerate(clusters):
        nbins = seg[seg["#ID"] == cluster]["#BINS"].iloc[0]
        dv_rd = np.abs(var_rd_matrix[i, :] - mv_rd)
        dv_baf = np.abs(var_baf_matrix[i, :] - mv_baf)
        z_rd = dv_rd / stdv_rd
        z_baf = dv_baf / stdv_baf
        is_outlier = (
            np.all(dv_rd > (fstd * stdv_rd))
            or np.all(dv_baf > (fstd * stdv_baf))
            or cluster_filtered[i]
        ) and (nbins <= ub_nbins)
        status = "REMOVED" if is_outlier else "kept"
        for j, sample in enumerate(samples):
            logging.debug(
                f"z={cluster} {sample} #bins={nbins} "
                f"RD-var={var_rd_matrix[i, j]:.6f} Z(RD)={z_rd[j]:.4f} "
                f"BAF-var={var_baf_matrix[i, j]:.6f} Z(BAF)={z_baf[j]:.4f} "
                f"{status}"
            )
        if is_outlier:
            bad_clusters.append(cluster)
        else:
            good_clusters.append(cluster)

    logging.info(f"remaining clusters: {good_clusters}")
    return good_clusters, bad_clusters


def compute_fractional_cn(rdr, baf, rdr_se, baf_se, gammas, alpha=0.05):
    """Compute fractional copy numbers and confidence intervals.

    fcn = gamma * rdr
    fb  = fcn * baf        (B-allele fractional CN)
    fa  = fcn * (1 - baf)  (A-allele fractional CN)

    SE propagation uses the delta method for products of independent MLEs:
        SE^2(X*Y) = mu_X^2 * SE_Y^2 + mu_Y^2 * SE_X^2 + SE_X^2 * SE_Y^2

    Returns a dict with keys: fcn, fa, fb, fa_lo, fa_hi, fb_lo, fb_hi.
    CI bounds use (1-alpha) normal approximation.
    """
    from scipy.stats import norm

    fcn = rdr * gammas
    fb = fcn * baf
    fa = fcn - fb

    g2 = gammas**2
    se_r2 = rdr_se**2
    se_b2 = baf_se**2
    r2 = rdr**2
    b2 = baf**2

    fb_se = np.sqrt(g2 * (r2 * se_b2 + b2 * se_r2 + se_r2 * se_b2))
    one_minus_b2 = (1 - baf) ** 2
    fa_se = np.sqrt(g2 * (r2 * se_b2 + one_minus_b2 * se_r2 + se_r2 * se_b2))

    z = norm.ppf(1 - alpha / 2)
    return {
        "fcn": fcn,
        "fa": fa,
        "fb": fb,
        "fa_lo": fa - z * fa_se,
        "fa_hi": fa + z * fa_se,
        "fb_lo": fb - z * fb_se,
        "fb_hi": fb + z * fb_se,
    }


def segmentation(
    cA,
    cB,
    u,
    cluster_ids,
    sample_ids,
    bbcs: pd.DataFrame,
    region_file: str,
    bbc_out_file=None,
    seg_out_file=None,
):
    """Annotate bins with inferred CN states and build a segment-level DataFrame.

    Merges the inferred allele-specific copy numbers (``cA``, ``cB``) and clone
    proportions (``u``) into the bin-level BBC DataFrame, then calls
    ``build_seg_from_bbc`` to merge adjacent bins with identical CN states into
    segments (respecting region boundaries).

    Args:
        cA: List of shape (num_clusters, num_clones) with allele-A CN integers.
        cB: List of shape (num_clusters, num_clones) with allele-B CN integers.
        u: List of shape (num_clones, num_samples) with clone proportions.
        cluster_ids: Ordered cluster identifiers matching the row order of cA/cB.
        sample_ids: Ordered sample identifiers matching the column order of u.
        bbcs: Bin-level DataFrame; must contain at least ``#CHR``, ``START``,
            ``END``, ``SAMPLE``, ``CLUSTER``.
        region_file: Path to the BED file of genomic regions used as segment
            merge barriers.
        bbc_out_file: If provided, write the annotated bin-level TSV to this path.
        seg_out_file: If provided, write the segment-level TSV to this path.

    Returns:
        The segment-level DataFrame (always returned regardless of whether
        ``seg_out_file`` is set).
    """
    df = bbcs.copy()

    n_clone = len(cA[0])
    cA = pd.DataFrame(cA, index=cluster_ids, columns=range(n_clone))
    cB = pd.DataFrame(cB, index=cluster_ids, columns=range(n_clone))
    u = pd.DataFrame(
        u, index=range(n_clone), columns=sample_ids
    ).T  # (n_sample, n_clone)

    cN = cA.astype(str) + "|" + cB.astype(str)
    cN.columns = ["cn_normal"] + [f"cn_clone{i}" for i in range(1, n_clone)]
    u.columns = ["u_normal"] + [f"u_clone{i}" for i in range(1, n_clone)]
    extra_columns = [col for pair in zip(cN.columns, u.columns) for col in pair]

    df = df.merge(cN, left_on="CLUSTER", right_index=True)
    df = df.merge(u, left_on="SAMPLE", right_index=True)
    df = df.sort_values(["#CHR", "START", "END", "SAMPLE"]).reset_index(drop=True)

    if bbc_out_file is not None:
        orig_cols = df.columns[: -2 * n_clone].tolist()
        df[orig_cols + extra_columns].to_csv(bbc_out_file, sep="\t", index=False)

    regions = read_region_bed(region_file)
    seg_df = build_seg_from_bbc(df, regions)
    if seg_out_file is not None:
        seg_df.to_csv(seg_out_file, sep="\t", index=False)

    return seg_df


def run_plot_cn(args, bbc, seg, gamma_file, plot_dir, ploidy):
    if not os.path.exists(bbc) or not os.path.exists(seg):
        return
    _plot_cn.run(
        {
            "bbc": bbc,
            "seg": seg,
            "genome_size": args["genome_size"],
            "region_bed": args["region_bed"],
            "gamma_file": gamma_file,
            "solfile": None,
            "plot_dir": plot_dir,
            "dpi": 150,
            "img_type": "png",
            "transparent": False,
            "keep_gap": False,
            "tail_alpha": 0.8,
            "center_alpha": 1.0,
            "onetail_area": 0.025,
            "maxlim_fcn": 30,
            "ploidy": ploidy,
        }
    )
