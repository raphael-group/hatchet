import os
import re
import glob
import logging
import pandas as pd
import numpy as np

from hatchet.utils import read_region_bed, build_seg_from_bbc
from hatchet.plot import plot_cn as _plot_cn


def store_gammas(out_file, scaling, samples):
    """Write per-sample gamma values for all ploidies.

    Args:
        out_file: output TSV path.
        scaling: dict from get_scaling_factor with 'diploid', 'tetraploid' keys.
        samples: ordered sample list.
    """
    with open(out_file, "w") as fd:
        for sample in samples:
            g_dip = scaling["diploid"]["gammas"].get(sample, 0)
            g_tet = (
                scaling["tetraploid"]["gammas"].get(sample, 0)
                if scaling["tetraploid"]
                else 0
            )
            fd.write(f"{sample}\t{g_dip}\t{g_tet}\n")


def build_data(bbcs, segs, segment=False):
    """Build (cluster or segment) x sample DataFrames for the solver.

    Args:
        bbcs: bin-level BBC DataFrame (needed only when segment=True).
        segs: cluster-level SEG DataFrame.
        segment: if True, expand clusters into genomic segments (maximal
            contiguous runs of same cluster on same chromosome).

    Returns dict with keys: rdr, baf, rdr_se, baf_se, nbins, weights,
    cluster_ids, sample_ids.  When segment=True, also includes
    seg_to_cluster and chr_boundaries.
    """
    if not segment:
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
            "cluster_ids": rdr.index.tolist(),
            "sample_ids": rdr.columns.tolist(),
        }

    # Segment mode
    bbcs = bbcs.sort_values(["SAMPLE", "#CHR", "START", "END"]).reset_index(drop=True)
    samples_sorted = sorted(bbcs["SAMPLE"].unique())
    first_sample = samples_sorted[0]

    mask = bbcs["SAMPLE"] == first_sample
    first_df = bbcs.loc[mask].reset_index(drop=True)
    seg_boundary = (first_df["CLUSTER"] != first_df["CLUSTER"].shift()) | (
        first_df["#CHR"] != first_df["#CHR"].shift()
    )
    seg_ids_first = seg_boundary.cumsum() - 1

    bbcs["_seg_int"] = np.tile(seg_ids_first.values, len(samples_sorted))
    agg = (
        bbcs.groupby(["_seg_int", "SAMPLE"])
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

    # String seg IDs: "<cluster>:chr<chrom>:<global_idx_within_cluster>"
    first_rows = agg.loc[agg["SAMPLE"] == first_sample].sort_values("_seg_int")
    cluster_counters = {}
    sid_map = {}
    for _, row in first_rows.iterrows():
        cid, chrom = row["CLUSTER"], row["CHR"]
        idx = cluster_counters.get(cid, 0)
        cluster_counters[cid] = idx + 1
        sid_map[row["_seg_int"]] = f"{cid}:{chrom}:{idx}"
    agg["_seg_id"] = agg["_seg_int"].map(sid_map)

    seg_cols = ["#ID", "SAMPLE", "RD", "BAF", "RD-se", "BAF-se"]
    agg = agg.merge(
        segs[seg_cols],
        left_on=["CLUSTER", "SAMPLE"],
        right_on=["#ID", "SAMPLE"],
        how="left",
    )

    rdr = agg.pivot(index="_seg_id", columns="SAMPLE", values="RD")
    baf = agg.pivot(index="_seg_id", columns="SAMPLE", values="BAF")
    rdr_se = agg.pivot(index="_seg_id", columns="SAMPLE", values="RD-se")
    baf_se = agg.pivot(index="_seg_id", columns="SAMPLE", values="BAF-se")
    nbins = agg.pivot(index="_seg_id", columns="SAMPLE", values="NBINS")

    # Sort by original segment order (not lexicographic string order)
    seg_order = [sid_map[i] for i in sorted(sid_map)]
    rdr = rdr.loc[seg_order]
    baf = baf.loc[seg_order]
    rdr_se = rdr_se.loc[seg_order]
    baf_se = baf_se.loc[seg_order]
    nbins = nbins.loc[seg_order]

    first_agg = (
        agg.loc[agg["SAMPLE"] == first_sample].set_index("_seg_id").loc[seg_order]
    )
    seg_lengths = first_agg["LENGTH"]
    weights = 100 * seg_lengths / seg_lengths.sum()
    seg_to_cluster = first_agg["CLUSTER"]
    chr_boundaries = (first_agg["CHR"] != first_agg["CHR"].shift()).values

    bbcs.drop(columns=["_seg_int"], inplace=True)

    n_segs, n_samples = rdr.shape
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
        "cluster_ids": rdr.index.tolist(),
        "sample_ids": rdr.columns.tolist(),
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


def store_solve_input(out_file, input_data):
    """Write solver input (FCN + weights) to a TSV."""
    weights = input_data["weights"]
    nbins = input_data["nbins"]
    fcn_cols = ["fcn", "fa", "fb", "fa_lo", "fa_hi", "fb_lo", "fb_hi"]
    header = "CLUSTER\tSAMPLE\t#BINS\t" + "\t".join(fcn_cols) + "\tweight"
    fa = input_data["fa"]
    with open(out_file, "w") as fd:
        fd.write(header + "\n")
        for sample in fa.columns:
            for cid in fa.index:
                nb = int(nbins.loc[cid, sample])
                vals = [str(input_data[c].loc[cid, sample]) for c in fcn_cols]
                fd.write(
                    f"{cid}\t{sample}\t{nb}\t" + "\t".join(vals) + f"\t{weights[cid]}\n"
                )


def _write_solution_tsv(fd, input_data, cA, cB, u, n, cluster_ids, sample_ids, header):
    """Write a single solution's per-cluster/sample details to an open file."""
    cA_ = np.array(cA)
    cB_ = np.array(cB)
    u_ = np.array(u)
    exp_a = cA_ @ u_
    exp_b = cB_ @ u_

    fd.write(header + "\n")
    for ci, cid in enumerate(cluster_ids):
        for si, sample in enumerate(sample_ids):
            fa_lo = input_data["fa_lo"].iloc[ci, si]
            fa_hi = input_data["fa_hi"].iloc[ci, si]
            fb_lo = input_data["fb_lo"].iloc[ci, si]
            fb_hi = input_data["fb_hi"].iloc[ci, si]
            ea, eb = exp_a[ci, si], exp_b[ci, si]
            accepted = ea >= fa_lo and ea <= fa_hi and eb >= fb_lo and eb <= fb_hi
            n_bins = int(input_data["nbins"].loc[cid, sample])
            fields = [
                cid,
                sample,
                n_bins,
                input_data["fa"].loc[cid, sample],
                input_data["fb"].loc[cid, sample],
                ea,
                eb,
                fa_lo,
                fa_hi,
                fb_lo,
                fb_hi,
            ]
            for oi in range(n):
                fields.extend([f"{cA[ci][oi]}|{cB[ci][oi]}", u[oi][si]])
            fields.append(accepted)
            fd.write("\t".join(str(v) for v in fields) + "\n")


def store_instance_tofile(pool_instances, input_data, sol_dir, solve_mode):
    """Store all solution detail TSVs (and Newick/JSON for cnt_cd)."""
    import json as _json

    n = len(pool_instances[next(iter(pool_instances))]["cA"][0])
    cluster_ids = input_data["cluster_ids"]
    sample_ids = input_data["sample_ids"]
    clone_cols = ["cn_normal\tu_normal"] + [
        f"cn_clone{i}\tu_clone{i}" for i in range(1, n)
    ]
    cols = (
        [
            "CLUSTER",
            "SAMPLE",
            "#BINS",
            "f_a",
            "f_b",
            "exp_f_a",
            "exp_f_b",
            "fa_lo",
            "fa_hi",
            "fb_lo",
            "fb_hi",
        ]
        + clone_cols
        + ["ci_accepted"]
    )
    header = "\t".join(cols)

    for sol_id, sol in pool_instances.items():
        path = os.path.join(sol_dir, f"{solve_mode}_{sol_id}.tsv")
        with open(path, "w") as fd:
            _write_solution_tsv(
                fd,
                input_data,
                sol["cA"],
                sol["cB"],
                sol["u"],
                n,
                cluster_ids,
                sample_ids,
                header,
            )

        if solve_mode == "cnt_cd":
            from hatchet.compute_cn.solve.cnt_tree import LabeledCloneTree

            tree = sol.get("tree")
            if tree is not None and isinstance(tree, LabeledCloneTree):
                prefix = os.path.join(sol_dir, f"{solve_mode}_{sol_id}")
                with open(f"{prefix}.nwk", "w") as f:
                    f.write(tree.to_newick() + "\n")
                d = tree.to_dict()
                d["imf_obj"] = sol.get("imf_obj")
                d["tree_obj"] = sol.get("tree_obj")
                d["u"] = sol.get("u")
                with open(f"{prefix}.json", "w") as f:
                    _json.dump(d, f, indent=2)


def compute_fractional_cn(input_data, gammas, alpha=0.05, min_ci_margin=0.1):
    """Compute fractional copy numbers and CI.

    Returns a new dict containing all fields from input_data plus
    fcn, fa, fb, fa_lo, fa_hi, fb_lo, fb_hi. input_data is not modified.

    Args:
        input_data: dict from build_data with rdr, baf, rdr_se.
        gammas: dict or Series of per-sample gamma values.
        alpha: significance level (default 0.05 → 95% CI).
        min_ci_margin: hard minimum CI half-width in FCN space.
    """
    from scipy.stats import norm

    rdr = input_data["rdr"]
    baf = input_data["baf"]
    rdr_se = input_data["rdr_se"]

    gammas = pd.Series(gammas).sort_index()
    fcn = rdr * gammas
    fb = fcn * baf
    fa = fcn - fb

    z = norm.ppf(1 - alpha / 2)
    margin_fa = np.maximum(z * gammas * (1 - baf) * rdr_se, min_ci_margin)
    margin_fb = np.maximum(z * gammas * baf * rdr_se, min_ci_margin)

    return {
        **input_data,
        "fcn": fcn,
        "fa": fa,
        "fb": fb,
        "fa_lo": fa - margin_fa,
        "fa_hi": fa + margin_fa,
        "fb_lo": fb - margin_fb,
        "fb_hi": fb + margin_fb,
    }


def segmentation(
    cA,
    cB,
    u,
    input_data: dict,
    bbcs: pd.DataFrame,
    region_file: str,
    bbc_out_file=None,
    seg_out_file=None,
):
    """Annotate bins with inferred CN states and build a segment-level DataFrame.

    Args:
        cA: (num_clusters, num_clones) allele-A CN.
        cB: (num_clusters, num_clones) allele-B CN.
        u: (num_clones, num_samples) clone proportions.
        input_data: dict from build_data with cluster_ids, sample_ids.
        bbcs: Bin-level DataFrame.
        region_file: Path to region BED file.
        bbc_out_file: If provided, write annotated bin-level TSV.
        seg_out_file: If provided, write segment-level TSV.

    Returns:
        Segment-level DataFrame.
    """
    cluster_ids = input_data["cluster_ids"]
    sample_ids = input_data["sample_ids"]
    seg_to_cluster = input_data.get("seg_to_cluster")
    df = bbcs.copy()

    n_clone = len(cA[0])
    cN = pd.DataFrame(
        np.array(cA).astype(str) + "|" + np.array(cB).astype(str),
        index=cluster_ids,
        columns=["cn_normal"] + [f"cn_clone{i}" for i in range(1, n_clone)],
    )
    u_df = pd.DataFrame(u, index=range(n_clone), columns=sample_ids).T
    u_df.columns = ["u_normal"] + [f"u_clone{i}" for i in range(1, n_clone)]
    extra_columns = [col for pair in zip(cN.columns, u_df.columns) for col in pair]

    if seg_to_cluster is not None:
        # Segment mode: assign each bin its segment ID, merge CN by segment.
        df = df.sort_values(["SAMPLE", "#CHR", "START", "END"]).reset_index(drop=True)
        samples_sorted = sorted(df["SAMPLE"].unique())
        first_mask = df["SAMPLE"] == samples_sorted[0]
        first_df = df.loc[first_mask].reset_index(drop=True)
        seg_boundary = (first_df["CLUSTER"] != first_df["CLUSTER"].shift()) | (
            first_df["#CHR"] != first_df["#CHR"].shift()
        )
        seg_int = (seg_boundary.cumsum() - 1).values
        df["_seg_int"] = np.tile(seg_int, len(samples_sorted))
        # Map integer seg index to string seg ID (same logic as build_data)
        seg_id_list = cluster_ids  # already in segment order
        int_to_seg = {i: seg_id_list[i] for i in range(len(seg_id_list))}
        df["_seg_id"] = df["_seg_int"].map(int_to_seg)
        df = df.merge(cN, left_on="_seg_id", right_index=True)
        df = df.drop(columns=["_seg_int", "_seg_id"])
    else:
        df = df.merge(cN, left_on="CLUSTER", right_index=True)

    df = df.merge(u_df, left_on="SAMPLE", right_index=True)
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

            pareto_mask = grp["is_pareto"]
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

            # Pareto points
            if len(pareto) > 0:
                ax.scatter(
                    pareto[reg_col],
                    pareto["IMF"],
                    c="#1f77b4",
                    s=50,
                    zorder=4,
                    label="Pareto",
                    edgecolors="white",
                    linewidths=0.5,
                )
                ax.plot(
                    pareto[reg_col],
                    pareto["IMF"],
                    c="black",
                    linewidth=1.5,
                    alpha=0.4,
                    zorder=3,
                )

            sel = grp[grp["selected"] == "*"]
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


def annotate_seg_pi_violations(seg_df, cA, cB, u, fcn_data, cluster_ids, sample_ids):
    """Add PI_VIOL column to seg_df using the prediction interval bounds in fcn_data.

    Each segment belongs to a cluster (CLUSTER column). A segment violates
    the PI when the expected FCN (cA @ u or cB @ u) falls outside the
    [fa_lo, fa_hi] / [fb_lo, fb_hi] bounds for that (cluster, sample).
    """
    if "CLUSTER" not in seg_df.columns or fcn_data is None or "fa_lo" not in fcn_data:
        return seg_df

    exp_a = np.array(cA) @ np.array(u)
    exp_b = np.array(cB) @ np.array(u)
    fa_lo = fcn_data["fa_lo"].loc[cluster_ids, sample_ids].to_numpy()
    fa_hi = fcn_data["fa_hi"].loc[cluster_ids, sample_ids].to_numpy()
    fb_lo = fcn_data["fb_lo"].loc[cluster_ids, sample_ids].to_numpy()
    fb_hi = fcn_data["fb_hi"].loc[cluster_ids, sample_ids].to_numpy()
    violations = (exp_a < fa_lo) | (exp_a > fa_hi) | (exp_b < fb_lo) | (exp_b > fb_hi)
    viol_df = pd.DataFrame(violations, index=cluster_ids, columns=sample_ids)

    seg_df = seg_df.copy()
    seg_df["PI_VIOL"] = seg_df.apply(
        lambda r: (
            bool(viol_df.loc[r["CLUSTER"], r["SAMPLE"]])
            if r["CLUSTER"] in viol_df.index
            else False
        ),
        axis=1,
    )
    return seg_df


def load_pool_from_disk(sol_dir, cluster_ids, sample_ids):
    """Read pool solution TSVs from sol_dir into {sol_id: {"imf_obj": ..., "cA": ..., ...}}."""
    pool = {}
    for path in sorted(glob.glob(os.path.join(sol_dir, "*.tsv"))):
        basename = os.path.basename(path)
        # Match old format (sol*_pool*) or new format (mode_solid)
        m = re.match(r".*_sol([\d.]+)_pool(\d+)\.tsv", basename)
        if m:
            sol_id = f"p{m.group(1)}_s{m.group(2)}"
        else:
            m2 = re.match(r"(?:cd|ilp|cnt_cd)_(.+)\.tsv", basename)
            if m2:
                sol_id = m2.group(1)
            else:
                continue

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
        pool[sol_id] = {"imf_obj": 0.0, "reg_obj": 0.0, "cA": cA, "cB": cB, "u": u}

    if pool:
        logging.info(f"loaded {len(pool)} pool solutions from {sol_dir}")
    return pool
