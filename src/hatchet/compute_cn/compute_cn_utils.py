import os
import re
import glob
import logging
import pandas as pd
import numpy as np

from hatchet.utils import read_region_bed, build_seg_from_bbc
from hatchet.plot import plot_cn as _plot_cn
from hatchet.compute_cn.solve.utils import (
    model_selection_instance,
    compute_individual_objs,
    compute_pairwise_cnt,
    filter_non_pareto,
    dedup_solutions,
)


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
        segs.loc[segs["SAMPLE"] == first_sample].set_index("#ID")["LENGTH"].sort_index()
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
        index to cluster ID), ``chr_boundaries`` (boolean array of length
        ``num_segments`` where True marks the first segment of each
        chromosome).
    """
    bbcs = bbcs.sort_values(["#CHR", "START", "END", "SAMPLE"]).reset_index(drop=True)
    samples_sorted = sorted(bbcs["SAMPLE"].unique())
    first_sample = samples_sorted[0]

    # Identify segment boundaries using only the first sample's rows
    mask = bbcs["SAMPLE"] == first_sample
    first_df = bbcs.loc[mask].reset_index(drop=True)

    seg_boundary = (first_df["CLUSTER"] != first_df["CLUSTER"].shift()) | (
        first_df["#CHR"] != first_df["#CHR"].shift()
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

    # Chromosome boundaries: True at the first segment of each chromosome
    chr_boundaries = (first_agg["CHR"] != first_agg["CHR"].shift()).values

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
        "chr_boundaries": chr_boundaries,
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


def compute_fractional_cn(rdr, baf, bbcs, gammas, alpha=0.5):
    """Compute fractional copy numbers and prediction intervals.

    FCN point estimates from cluster-level rdr/baf.
    Prediction interval width from std of per-bin FA/FB values.

    Returns a dict with keys: fcn, fa, fb, fa_lo, fa_hi, fb_lo, fb_hi.
    """
    from scipy.stats import norm

    fcn = rdr * gammas
    fb = fcn * baf
    fa = fcn - fb

    z = norm.ppf(1 - alpha / 2)
    cluster_ids = rdr.index.tolist()
    sample_ids = rdr.columns.tolist()

    fa_std = pd.DataFrame(0.0, index=rdr.index, columns=rdr.columns)
    fb_std = pd.DataFrame(0.0, index=rdr.index, columns=rdr.columns)

    for cid in cluster_ids:
        for sid in sample_ids:
            bins = bbcs[(bbcs["CLUSTER"] == cid) & (bbcs["SAMPLE"] == sid)]
            if len(bins) == 0:
                continue
            g = gammas[sid] if hasattr(gammas, "__getitem__") else gammas
            fcn_bins = g * bins["RD"].to_numpy()
            fa_bins = fcn_bins * (1 - bins["BAF"].to_numpy())
            fb_bins = fcn_bins * bins["BAF"].to_numpy()
            fa_std.loc[cid, sid] = np.std(fa_bins, ddof=1) if len(fa_bins) > 1 else 0.0
            fb_std.loc[cid, sid] = np.std(fb_bins, ddof=1) if len(fb_bins) > 1 else 0.0

    return {
        "fcn": fcn,
        "fa": fa,
        "fb": fb,
        "fa_lo": fa - z * fa_std,
        "fa_hi": fa + z * fa_std,
        "fb_lo": fb - z * fb_std,
        "fb_hi": fb + z * fb_std,
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
            "img_type": "pdf",
            "transparent": False,
            "keep_gap": False,
            "tail_alpha": 0.8,
            "center_alpha": 1.0,
            "onetail_area": 0.025,
            "maxlim_fcn": 30,
            "ploidy": ploidy,
            "style": args.get("style", "cnv"),
        }
    )


def plot_pareto_pdf(summary_df, plot_dir, reg_term, elbow_fig=None):
    """Plot REG vs IMF Pareto curves + elbow/BIC page as a multi-page PDF."""
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    outfile = os.path.join(plot_dir, "model_selection.pdf")
    reg_col = reg_term if reg_term in summary_df.columns else "REG"
    groups = sorted(summary_df.groupby(["ploidy", "n_clones"]))

    with PdfPages(outfile) as pdf:
        for (ploidy, n_clones), grp in groups:
            fig, ax = plt.subplots(figsize=(7, 5))

            pareto_mask = grp["is_pareto"] == True
            non_pareto = grp[~pareto_mask]
            pareto = grp[pareto_mask].sort_values(reg_col)

            # Non-pareto: gray
            if len(non_pareto) > 0:
                ax.scatter(
                    non_pareto[reg_col],
                    non_pareto["IMF"],
                    c="0.75",
                    s=25,
                    zorder=2,
                    alpha=0.5,
                    edgecolors="white",
                    linewidths=0.3,
                )

            # Pareto points: colored by CNT feasibility
            if len(pareto) > 0:
                if "CNT_from_c1" in pareto.columns:
                    has_inf = pareto["CNT_from_c1"].apply(
                        lambda v: (
                            v == "inf" or (isinstance(v, float) and not np.isfinite(v))
                        )
                    )
                else:
                    has_inf = pd.Series(False, index=pareto.index)

                finite_p = pareto[~has_inf]
                inf_p = pareto[has_inf]

                if len(finite_p) > 0:
                    ax.scatter(
                        finite_p[reg_col],
                        finite_p["IMF"],
                        c="#1f77b4",
                        s=50,
                        zorder=4,
                        label="Pareto",
                        edgecolors="white",
                        linewidths=0.5,
                    )
                if len(inf_p) > 0:
                    ax.scatter(
                        inf_p[reg_col],
                        inf_p["IMF"],
                        c="#d62728",
                        s=50,
                        marker="x",
                        zorder=4,
                        linewidths=1.5,
                        label="Pareto (CNT inf)",
                    )
                ax.plot(
                    pareto[reg_col],
                    pareto["IMF"],
                    c="black",
                    linewidth=1.5,
                    alpha=0.4,
                    zorder=3,
                )

            sel = grp[grp["is_instance_selected"] == True]
            if len(sel) > 0:
                ax.scatter(
                    sel[reg_col],
                    sel["IMF"],
                    c="gold",
                    marker="*",
                    s=250,
                    zorder=5,
                    edgecolors="black",
                    linewidths=1,
                    label="selected",
                )

            ax.set_xlabel(reg_col, fontsize=11)
            ax.set_ylabel("IMF", fontsize=11)
            ax.set_title(
                f"{ploidy} n={n_clones} ({len(grp)} solutions)",
                fontsize=13,
                fontweight="bold",
            )
            ax.legend(fontsize=9)
            ax.grid(True, alpha=0.3)
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)

        # Append elbow/BIC figure as last page
        if elbow_fig is not None:
            pdf.savefig(elbow_fig)
            plt.close(elbow_fig)

    logging.info(f"wrote {outfile} ({len(groups) + (1 if elbow_fig else 0)} pages)")


def pool_entries_for_plot(pool):
    """Extract plot-ready tuples from pool dict, dropping cnt_pairs."""
    return [
        (tag, seg_df, imf, pareto, selected)
        for tag, (seg_df, imf, _reg, pareto, selected, _cnt) in pool.items()
    ]


def load_pool_from_disk(sol_dir, cluster_ids, sample_ids):
    """Read pool solution TSVs from sol_dir into {pparam: [(obj, cA, cB, u), ...]}."""
    pool = {}
    for path in sorted(glob.glob(os.path.join(sol_dir, "*_sol*_pool*.tsv"))):
        m = re.match(r".*_sol([\d.]+)_pool(\d+)\.tsv", os.path.basename(path))
        if not m:
            continue
        pparam = float(m.group(1)) if "." in m.group(1) else int(m.group(1))

        sol = pd.read_csv(path, sep="\t")
        cn_cols = sorted(
            [c for c in sol.columns if c.startswith("cn_")],
            key=lambda c: (0 if c == "cn_normal" else 1, c),
        )
        u_cols = sorted(
            [c for c in sol.columns if c.startswith("u_")],
            key=lambda c: (0 if c == "u_normal" else 1, c),
        )

        sol_s = (
            sol[sol["SAMPLE"] == sample_ids[0]]
            .sort_values("CLUSTER")
            .reset_index(drop=True)
        )
        cA, cB = [], []
        for _, row in sol_s.iterrows():
            ca, cb = zip(
                *(
                    (int(a), int(b))
                    for a, b in (str(row[c]).split("|") for c in cn_cols)
                )
            )
            cA.append(list(ca))
            cB.append(list(cb))

        u = [
            [float(sol[sol["SAMPLE"] == sid].iloc[0][uc]) for sid in sample_ids]
            for uc in u_cols
        ]
        pool.setdefault(pparam, []).append((0.0, cA, cB, u))

    if pool:
        logging.info(
            f"loaded {sum(len(v) for v in pool.values())} pool solutions from {sol_dir}"
        )
    return pool


def build_pool_output(
    pool_instances,
    f_a,
    f_b,
    fcn_data,
    weights,
    nbins,
    args,
    cluster_ids,
    sample_ids,
    bbcs,
    out_bbc,
    out_seg,
    sol_dir,
):
    """Run model selection on pool_instances and build the pool output dict."""
    reg_term = args["reg_term"]

    best_instance, imf_obj, selected_key = model_selection_instance(
        f_a,
        f_b,
        weights,
        pool_instances,
        reg_term,
        args["mode"],
        sol_dir,
        fcn_data,
        nbins,
    )
    if best_instance is None:
        return 0.0, 0.0, {}

    obj, cA, cB, u = best_instance
    if not os.path.exists(out_bbc) or not os.path.exists(out_seg):
        segmentation(
            cA,
            cB,
            u,
            cluster_ids,
            sample_ids,
            bbcs=bbcs,
            region_file=args["region_bed"],
            bbc_out_file=out_bbc,
            seg_out_file=out_seg,
        )

    all_pool = {}
    pool_objs = []
    pool_tags = []
    pool_keys = []
    for pparam, sols in pool_instances.items():
        for pidx, (pobj, pcA, pcB, pu) in enumerate(sols):
            tag = f"pool_p{pparam}_s{pidx}"
            seg_df = segmentation(
                pcA,
                pcB,
                pu,
                cluster_ids,
                sample_ids,
                bbcs=bbcs,
                region_file=args["region_bed"],
            )
            p_imf, p_reg = compute_individual_objs(
                reg_term, weights, f_a, f_b, pcA, pcB, pu
            )
            cnt_pairs = compute_pairwise_cnt(pcA, pcB, bbcs, cluster_ids)
            pool_objs.append([p_imf, p_reg])
            pool_tags.append(tag)
            pool_keys.append((pparam, pidx))
            all_pool[tag] = (seg_df, p_imf, p_reg, False, False, cnt_pairs)

    if pool_objs:
        is_pareto = filter_non_pareto(np.array(pool_objs))
        for tag, key, pareto in zip(pool_tags, pool_keys, is_pareto):
            seg_df_, imf_, reg_, _, _, cnt_ = all_pool[tag]
            is_selected = key == selected_key
            all_pool[tag] = (seg_df_, imf_, reg_, bool(pareto), is_selected, cnt_)

    return obj, imf_obj, all_pool


def dedup_pool(pool_instances):
    """Deduplicate pool_instances dict across all pparam values."""
    flat_sols, flat_keys = [], []
    for pparam, sols in pool_instances.items():
        for pidx, sol in enumerate(sols):
            flat_sols.append(sol)
            flat_keys.append((pparam, pidx))

    deduped = dedup_solutions(flat_sols)
    deduped_ids = {id(s) for s in deduped}

    out = {}
    for sol, (pparam, _) in zip(flat_sols, flat_keys):
        if id(sol) in deduped_ids:
            out.setdefault(pparam, []).append(sol)
    return out
