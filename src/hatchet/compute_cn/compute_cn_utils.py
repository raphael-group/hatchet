import os
import sys
import logging
import pandas as pd
import numpy as np

from hatchet.utils import read_region_bed, sort_df_chr
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


# clonal cluster z (a,b)
def get_purity_by_baf(bafz: float, a: int, b: int):
    """
    compute tumor purity by BAF and copy-number (a,b)
    """
    num = 2 * bafz - 1
    dom = (b - 1) - bafz * (a + b - 2)
    return -1 if dom == 0 else num / dom


def get_purity_by_rrd(rrdz: float, a: int, b: int, is_wgd=True):
    """
    compute tumor purity by RRD and copy-number (a,b)
    assumes a + b > 2 when no WGD
    """
    num = 2 * rrdz - 2
    if is_wgd:
        dom = a + b - 2 - 2 * rrdz
    else:
        dom = a + b - 2
    return -1 if dom == 0 else num / dom


def purity_est_err(bafz: float, rrdz: float, a: int, b: int, is_wgd: bool):
    """
    estimate tumor purity by 1) BAF and 2) RRD with given cn (a,b)
    return: estimation error
    """
    pbaf = get_purity_by_baf(bafz, a, b)
    prrd = get_purity_by_rrd(rrdz, a, b, is_wgd)
    if pbaf <= 0.0 or pbaf > 1.0:
        return np.inf, pbaf, prrd
    if prrd <= 0.0 or prrd > 1.0:
        return np.inf, pbaf, prrd
    return abs(pbaf - prrd), pbaf, prrd


def get_gamma_WGD(rds: float, rdz: float, cz: int):
    dom = (cz - 2) * rds - 2 * rdz
    if dom == 0.0:
        return -1
    return (2 * cz - 8) / dom


def _fit_cn(z, cn_candidates, is_wgd, samples, baf, rrdr, tol_err, is_cand_loh):
    """
    Try each (a,b) in cn_candidates using purity_est_err across all samples.
    Returns ((a,b), purities_dict, min_err) on first match, or (None, None, min_err).
    """
    min_err = np.inf
    for a, b in cn_candidates:
        if b == 0 and not is_cand_loh[z]:
            continue
        purities = {}
        for sample in samples:
            perr, pbaf, prrd = purity_est_err(
                baf.loc[z, sample], rrdr.loc[z, sample], a, b, is_wgd
            )
            min_err = min(min_err, perr)
            if perr <= tol_err:
                purities[sample] = (pbaf + prrd) / 2
            else:
                purities = None
                break
        if purities is not None:
            return (a, b), purities, min_err
    return None, None, min_err


def get_scaling_factor(
    samples: list,
    seg: pd.DataFrame,
    balanced_s: list,
    unbalanced_z: list,
    tol_rd_ratio: float,
    tol_err: float,
    maxcn: int,
    maxcn_wgd: int,
):
    """
    Compute scaling factors
    """
    gammas_noWGD = {}
    purities_noWGD = {}
    pair_noWGD = None

    gammas_WGD = {}
    purities_WGD = {}
    pair_WGD = None

    rdr = seg.pivot(index="#ID", columns="SAMPLE", values="RD")
    baf = seg.pivot(index="#ID", columns="SAMPLE", values="BAF")
    balanced_s = sorted(balanced_s, key=lambda s: rdr.loc[s, :].mean())

    if len(balanced_s) >= 2:
        # TODO also reason about (0,0) or (2,2) base?
        s0, s1 = balanced_s[0], balanced_s[1]
        pair_noWGD = (s0, s1, (1, 1), (2, 2))
        logging.info(
            f"found >1 balanced clusters, assign (1,1) and (2,2) to {s0} and {s1}"
        )
        for sample in samples:
            gamma = 2 / rdr.loc[s0, sample]
            gammas_noWGD[sample] = gamma
            purity = (0.5 * rdr.loc[s1, sample] / gamma) - 1
            purities_noWGD[sample] = purity
        logging.info(
            f"estimated purities (noWGD): { {s: f'{p:.4f}' for s, p in purities_noWGD.items()} }"
        )
        return (
            s0,
            pair_noWGD,
            gammas_noWGD,
            purities_noWGD,
            pair_WGD,
            gammas_WGD,
            purities_WGD,
        )

    # found one unbalanced cluster z that pairs with s0 across all samples
    s0 = balanced_s[0]
    assert np.all(rdr.loc[s0, :] > 0), f"balanced cluster {s0} has RD<=0, invalid value"
    rrdr = rdr / rdr.loc[s0, :]

    is_cand_loh = {z: True for z in unbalanced_z}
    # LOH cluster cannot have any lower-RD cluster also has lower BAF.
    for z in unbalanced_z:
        for _z in unbalanced_z:
            if np.all(rdr.loc[_z] < rdr.loc[z]):
                if np.any(baf.loc[_z] < baf.loc[z]):
                    is_cand_loh[z] = False
                    break

    cn_nowgd_above = []
    cn_wgd_above = []
    # (3,0),(2,1),(4,0),(3,1),(5,0),(4,1),(3,2)...
    for c in range(3, max(maxcn, maxcn_wgd) + 1):
        for b in range(0, (c // 2) + 1):
            a = c - b
            if a == b:
                continue
            if c <= maxcn:
                cn_nowgd_above.append((a, b))
            if c <= maxcn_wgd and c >= 5:
                cn_wgd_above.append((a, b))

    cn_nowgd_below = [(1, 0)]
    cn_wgd_below = [(1, 0), (2, 0), (3, 0), (2, 1)]

    logging.debug(
        f"candiate un-paired clonal states (no WGD): {cn_nowgd_below},{cn_nowgd_above}"
    )
    logging.debug(
        f"candiate un-paired clonal states (WGD): {cn_wgd_below},{cn_wgd_above}"
    )

    err_noWGD = np.inf
    err_WGD = np.inf
    for z in sorted(unbalanced_z, key=lambda z: baf.loc[z, :].mean()):
        rd_ratio_zs = rdr.loc[z] / rdr.loc[s0]
        rd_dist_zs = rdr.loc[z] - rdr.loc[s0]
        if np.all(np.abs(rd_ratio_zs - 1) <= tol_rd_ratio):
            logging.debug(f"-----------z={z} pair {s0}")
            # case 1: s0=(1,1), z=(2,0) — BAF only (RD ≈ same, RRD uninformative)
            if pair_noWGD is None:
                purities = {}
                for sample in samples:
                    pbaf = get_purity_by_baf(baf.loc[z, sample], 2, 0)
                    if 0.0 < pbaf <= 1.0:
                        purities[sample] = pbaf
                    else:
                        purities = None
                        break
                if purities is not None:
                    pair_noWGD = (s0, z, (1, 1), (2, 0))
                    purities_noWGD = purities

            # case 2: s0=(2,2), z=(4,0)
            if pair_WGD is None:
                cn, purs, err = _fit_cn(
                    z, [(4, 0)], True, samples, baf, rrdr, tol_err, is_cand_loh
                )
                err_WGD = min(err_WGD, err)
                if cn is not None:
                    pair_WGD = (s0, z, (2, 2), cn)
                    purities_WGD = purs
        elif np.all(rd_dist_zs > 0):
            logging.debug(f"-----------z={z} above {s0}")
            if pair_noWGD is None:
                cn, purs, err = _fit_cn(
                    z, cn_nowgd_above, False, samples, baf, rrdr, tol_err, is_cand_loh
                )
                err_noWGD = min(err_noWGD, err)
                if cn is not None:
                    pair_noWGD = (s0, z, (1, 1), cn)
                    purities_noWGD = purs

            if pair_WGD is None:
                cn, purs, err = _fit_cn(
                    z, cn_wgd_above, True, samples, baf, rrdr, tol_err, is_cand_loh
                )
                err_WGD = min(err_WGD, err)
                if cn is not None:
                    pair_WGD = (s0, z, (2, 2), cn)
                    purities_WGD = purs
        elif np.all(rd_dist_zs < 0):
            logging.debug(f"-----------z={z} below {s0}")
            if pair_noWGD is None:
                cn, purs, err = _fit_cn(
                    z, cn_nowgd_below, False, samples, baf, rrdr, tol_err, is_cand_loh
                )
                err_noWGD = min(err_noWGD, err)
                if cn is not None:
                    pair_noWGD = (s0, z, (1, 1), cn)
                    purities_noWGD = purs

            if pair_WGD is None:
                cn, purs, err = _fit_cn(
                    z, cn_wgd_below, True, samples, baf, rrdr, tol_err, is_cand_loh
                )
                err_WGD = min(err_WGD, err)
                if cn is not None:
                    pair_WGD = (s0, z, (2, 2), cn)
                    purities_WGD = purs
        else:
            logging.debug(
                f"cluster {z} has inconsistent relative position to {s0} across samples"
            )
            logging.debug(f"RD-distance(z,s)={rd_dist_zs}")
            logging.debug(f"RD-ratio(z,s)={rd_ratio_zs}")

        if pair_noWGD is not None and pair_WGD is not None:
            break

    logging.debug(f"lowest purity-esimation error (noWGD)={err_noWGD}, bound={tol_err}")
    logging.debug(f"lowest purity-esimation error (WGD)={err_WGD}, bound={tol_err}")

    # in noWGD case, pair is not required.
    for sample in samples:
        gamma = 2 / rdr.loc[s0, sample]
        gammas_noWGD[sample] = gamma

    if pair_WGD is not None:
        (_, z, (_, _), (za, zb)) = pair_WGD
        for sample in samples:
            gamma = get_gamma_WGD(rdr.loc[s0, sample], rdr.loc[z, sample], za + zb)
            gammas_WGD[sample] = gamma

    if purities_noWGD:
        logging.info(
            f"estimated purities (noWGD): { {s: f'{p:.4f}' for s, p in purities_noWGD.items()} }"
        )
    if purities_WGD:
        logging.info(
            f"estimated purities (WGD): { {s: f'{p:.4f}' for s, p in purities_WGD.items()} }"
        )
    return (
        s0,
        pair_noWGD,
        gammas_noWGD,
        purities_noWGD,
        pair_WGD,
        gammas_WGD,
        purities_WGD,
    )


##################################################
def annotate_balanced_clusters(
    segs: pd.DataFrame, balanced_baf_tol: float, colname="balanced"
):
    baf_mat = segs.pivot(index="#ID", columns="SAMPLE", values="BAF")
    clusters = baf_mat.index.tolist()
    balanced_s = []
    for cid in clusters:
        # balance cluster over all samples
        if np.all(np.abs(baf_mat.loc[cid, :] - 0.5) <= balanced_baf_tol):
            balanced_s.append(cid)

    segs[colname] = segs["#ID"].isin(balanced_s)
    num_balanced = len(balanced_s)
    if num_balanced == 0:
        potential_baf_tol = min(np.max(np.abs(baf_mat - 0.5), axis=1))
        logging.info(
            f"failed to locate balanced clusters, increase balanced_baf_tol to at least {potential_baf_tol:.5f}"
        )
        sys.exit(1)
    return segs


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
        df["all_copy_numbers"] = df[cN.columns].apply(",".join, axis=1)
        first_sample = df["SAMPLE"].iloc[0]
        df["segment"] = (
            (df["SAMPLE"] == first_sample)
            & (
                (df["#CHR"] != df["#CHR"].shift())
                | (df["all_copy_numbers"] != df["all_copy_numbers"].shift())
                | (df["START"] != df["END"].shift())
            )
        ).cumsum()

        agg = {"#CHR": "first", "START": "min", "END": "max", "SAMPLE": "first"}
        agg.update({c: "first" for c in extra_columns})
        seg = df.groupby(["segment", "SAMPLE"]).agg(agg)

        regions = read_region_bed(region_file)

        # Merge adjacent same-CN segments within each region (regions act as merge boundaries).
        # All segments are preserved; two adjacent same-CN rows are only merged if they fall
        # within the same region.
        seg = seg.reset_index(drop=True)
        cn_cols = [c for c in extra_columns if c.startswith("cn_")]
        out_cols = ["#CHR", "START", "END", "SAMPLE"] + extra_columns

        # Assign each segment to a region index (-1 = outside all regions)
        seg["_region"] = -1
        for r_idx, region in regions.iterrows():
            mask = (
                (seg["#CHR"] == region["#CHR"])
                & (seg["START"] >= region["START"])
                & (seg["END"] <= region["END"])
            )
            seg.loc[mask, "_region"] = r_idx

        merged_rows = []
        for _, grp in seg.groupby(["SAMPLE", "#CHR", "_region"], sort=False):
            grp = grp.sort_values("START").reset_index(drop=True)
            state_key = grp[cn_cols].apply(tuple, axis=1)
            grp["_run"] = (state_key != state_key.shift()).cumsum()
            for _, run_grp in grp.groupby("_run"):
                row = run_grp.iloc[0].copy()
                row["START"] = run_grp["START"].min()
                row["END"] = run_grp["END"].max()
                merged_rows.append(row[out_cols])

        out = sort_df_chr(pd.DataFrame(merged_rows, columns=out_cols), pos="START")
        out = out.sort_values(["#CHR", "START", "SAMPLE"]).reset_index(drop=True)
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
