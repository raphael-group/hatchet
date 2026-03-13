import os
import logging
import pandas as pd
import numpy as np

from hatchet.utils import read_region_bed, sort_df_chr, build_seg_from_bbc
from hatchet.plot import plot_cn as _plot_cn


def filtering(
    bbc: pd.DataFrame,
    seg: pd.DataFrame,
    samples: list,
    clusters: list,
    fstd=2.0,
    min_nbins=10,
    ub_nbins=50,
):
    """
    filter&merge clusters before optimization step
    0. filter any cluster has #bins<min_nbins
    1. compute per-sample per-cluster variance SCV,
    2. compute per-sample MV and STDV
    3. filter a cluster if it has |SCV - MV| >= 2 * STDV and below ub_nbins for all samples.

    Returns:
    1. list of remaining cluster IDs
    """
    logging.info(f"preprocessing, filtering clusters")

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
            f"{sample}\tRD-variance bound={(lb_rd, ub_rd)}\tBAF-variance bound={(lb_baf, ub_baf)}"
        )

    good_clusters = []
    bad_clusters = []
    for i, cluster in enumerate(clusters):
        nbins = seg[seg["#ID"] == cluster]["#BINS"].iloc[0]
        dv_rd = np.abs(var_rd_matrix[i, :] - mv_rd)
        dv_baf = np.abs(var_baf_matrix[i, :] - mv_baf)
        logging.debug(
            f"\t#ID={cluster}\t#bins={nbins}\tRD-variance={var_rd_matrix[i, :]}\tBAF-variance={var_baf_matrix[i, :]}"
        )
        logging.debug(f"\tZ(RD)={dv_rd / stdv_rd}\tZ(BAF)={dv_baf / stdv_baf}")
        if (
            np.all(dv_rd > (fstd * stdv_rd))
            or np.all(dv_baf > (fstd * stdv_baf))
            or cluster_filtered[i]
        ) and (nbins <= ub_nbins):
            logging.info(f"cluster {cluster} is outlier, removed")
            bad_clusters.append(cluster)
        else:
            good_clusters.append(cluster)

    logging.info(f"remaining clusters: {good_clusters}")
    return good_clusters, bad_clusters


##################################################
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

    if seg_out_file is not None:
        regions = read_region_bed(region_file)
        out = build_seg_from_bbc(df, regions)
        out.to_csv(seg_out_file, sep="\t", index=False)


##################################################
def run_plot_cn(args, out_dir, plot_dir, gamma_file, ploidy, n):
    bbc = os.path.join(out_dir, f"results.{ploidy}.n{n}.bbc.ucn.tsv")
    seg = os.path.join(out_dir, f"results.{ploidy}.n{n}.seg.ucn.tsv")
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
            "plot_dir": os.path.join(plot_dir, f"{ploidy}_n{n}"),
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
