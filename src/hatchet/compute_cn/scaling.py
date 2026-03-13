import sys
import logging
import numpy as np
import pandas as pd


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
    balanced_baf_tol: float,
    tol_rd_ratio: float,
    tol_err: float,
    maxcn: int,
    maxcn_wgd: int,
):
    """Compute scaling factors.

    First identifies balanced/unbalanced clusters via balanced_baf_tol,
    then estimates gamma and purity for noWGD and WGD scenarios.

    Returns (s0, pair_noWGD, gammas_noWGD, purities_noWGD,
             pair_WGD, gammas_WGD, purities_WGD,
             balanced_s, unbalanced_z).
    """
    # identify balanced / unbalanced clusters
    rdr = seg.pivot(index="#ID", columns="SAMPLE", values="RD")
    baf = seg.pivot(index="#ID", columns="SAMPLE", values="BAF")
    baf_mat = seg.pivot(index="#ID", columns="SAMPLE", values="BAF")
    clusters = baf_mat.index.tolist()
    balanced_s = []
    for cid in clusters:
        if np.all(np.abs(baf_mat.loc[cid, :] - 0.5) <= balanced_baf_tol):
            balanced_s.append(cid)

    if len(balanced_s) == 0:
        potential_baf_tol = min(np.max(np.abs(baf_mat - 0.5), axis=1))
        logging.error(
            f"failed to locate balanced clusters, increase balanced_baf_tol to at least {potential_baf_tol:.5f}"
        )
        sys.exit(1)

    balanced_set = set(balanced_s)
    unbalanced_z = [cid for cid in clusters if cid not in balanced_set]

    gammas_noWGD = {}
    purities_noWGD = {}
    pair_noWGD = None

    gammas_WGD = {}
    purities_WGD = {}
    pair_WGD = None

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
            balanced_s,
            unbalanced_z,
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
        balanced_s,
        unbalanced_z,
    )
