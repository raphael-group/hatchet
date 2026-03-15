import math
import sys
import logging
import numpy as np
import pandas as pd
from scipy.stats import norm


def _build_cn_candidates(maxcn):
    """Build imbalanced CN candidates: (a, b) with a > b, a + b <= maxcn."""
    return [(1, 0)] + [
        (c - b, b) for c in range(2, maxcn + 1) for b in range(c // 2 + 1) if c - b != b
    ]


def _score_pair(
    purities,
    gammas,
    clusters,
    samples,
    baf,
    baf_tau,
    rdr,
    rd_std,
    nbins_total,
    maxcn,
    k,
):
    """Score a (s0, z) pair by total #BINS of clusters it can explain.

    A cluster j is "explained" if there exists integer CN (a, b) with a >= 0, b >= 0,
    and a + b <= maxcn such that for every sample p, the expected BAF and RDR
    (given purities[p] and gammas[p]) fall within k * std of the observed values.

    BAF std is computed from the expected BAF: sqrt(p*(1-p) / (tau_p + 1)).
    RDR std is provided directly (sqrt of within-cluster RDR variance).
    """
    score = 0
    for j in clusters:
        for c in range(maxcn + 1):
            explained = False
            for b in range(c + 1):
                a = c - b
                ok = True
                for si, s in enumerate(samples):
                    tau = purities[si]
                    denom = 2 * (1 - tau) + (a + b) * tau
                    if denom <= 0:
                        ok = False
                        break
                    exp_baf = (1 - tau + b * tau) / denom
                    exp_rdr = denom / gammas[si]
                    baf_std = math.sqrt(
                        exp_baf * (1 - exp_baf) / (baf_tau.loc[j, s] + 1)
                    )
                    if abs(baf.loc[j, s] - exp_baf) > k * baf_std:
                        ok = False
                        break
                    if abs(rdr.loc[j, s] - exp_rdr) > k * rd_std.loc[j, s]:
                        ok = False
                        break
                if ok:
                    score += nbins_total[j]
                    explained = True
                    break
            if explained:
                break
    return score


def get_scaling_factor(
    samples: list,
    seg: pd.DataFrame,
    bal_tost_alpha: float,
    bal_tost_margin: float,
    tol_nstd: float,
    tolerance: float,
    maxcn: int,
    maxcn_wgd: int,
):
    """Infer RDR scaling factors (gamma) and tumor purities.

    Steps:
      1. Classify clusters as balanced (BAF ≈ 0.5) or imbalanced via TOST
         equivalence test. Select base cluster s0 = largest balanced cluster.
      2. Compute noWGD gamma per sample: gamma_p = 2 / RDR(s0, p).
      3. For each imbalanced cluster z and candidate CN (a, b), estimate purity
         from BAF and RRD. When a+b equals base ploidy (2 noWGD, 4 WGD), RRD is
         degenerate so BAF-only purity is used. For WGD, also compute gamma from
         the (s0, z) RDR pair.
      4. Score each valid (z, a, b) pair by counting total #BINS of clusters
         explainable by some integer CN state within tol_nstd * std of observed
         BAF/RDR. Select the highest-scoring pair for noWGD and WGD separately.

    Returns (s0, pair_noWGD, gammas_noWGD, purities_noWGD,
             pair_WGD, gammas_WGD, purities_WGD).
    """
    logging.info("Infer scaling factors & tumor purity")

    rdr = seg.pivot(index="#ID", columns="SAMPLE", values="RD")
    baf = seg.pivot(index="#ID", columns="SAMPLE", values="BAF")
    baf_se = seg.pivot(index="#ID", columns="SAMPLE", values="BAF-se")
    baf_tau = seg.pivot(index="#ID", columns="SAMPLE", values="BAF-tau")
    rd_var = seg.pivot(index="#ID", columns="SAMPLE", values="RD-var")
    rd_std = np.sqrt(rd_var)
    nbins = seg.pivot(index="#ID", columns="SAMPLE", values="#BINS")
    clusters = baf.index.tolist()

    balanced_s = []
    imbalanced_z = []
    for cid in clusters:
        is_balanced = True
        for sample in samples:
            diff = baf.loc[cid, sample] - 0.5
            se = max(baf_se.loc[cid, sample], 1e-12)
            p_tost = max(
                norm.cdf((diff - bal_tost_margin) / se),
                norm.sf((diff + bal_tost_margin) / se),
            )
            if p_tost >= bal_tost_alpha:
                is_balanced = False
                break
        if is_balanced:
            balanced_s.append(cid)
        else:
            imbalanced_z.append(cid)

    if len(balanced_s) == 0:
        logging.error(
            "failed to locate balanced clusters via TOST equivalence test "
            f"(alpha={bal_tost_alpha}, margin={bal_tost_margin}). "
            "Consider increasing bal_tost_margin or checking clustering results."
        )
        sys.exit(1)

    s0 = max(balanced_s, key=lambda s: nbins.loc[s, samples].sum())
    assert np.all(rdr.loc[s0, :] > 0), f"balanced cluster {s0} has RD<=0"
    logging.info(f"balanced clusters={balanced_s}, base cluster s0={s0}")

    gammas_noWGD = {s: 2.0 / rdr.loc[s0, s] for s in samples}
    logging.info(f"inferred scaling factor for no-WGD case {gammas_noWGD}")
    assert all(g > 0 for g in gammas_noWGD.values()), (
        "at least one sample has invalid scaling factor"
    )

    if len(imbalanced_z) == 0:
        logging.warning(
            "no unbalanced clusters found, skip purity estimation & WGD case"
        )
        return (s0, None, gammas_noWGD, None, None, None, None)

    cn_nowgd_all = _build_cn_candidates(maxcn)
    cn_wgd_all = _build_cn_candidates(maxcn_wgd)

    nbins_total = nbins[samples].sum(axis=1)
    gammas_nowgd_arr = np.array([gammas_noWGD[s] for s in samples])

    valid_nowgd = {}
    valid_wgd = {}

    rrdr = rdr / rdr.loc[s0, :]
    for z in imbalanced_z:
        if not (
            np.all(baf.loc[z, :] < 0.5 + bal_tost_margin)
            or np.all(baf.loc[z, :] > 0.5 - bal_tost_margin)
        ):
            logging.debug(
                f"cluster {z} has inconsistent BAF side across samples, skipping"
            )
            continue
        is_major = np.all(baf.loc[z, :] > 0.5)

        baf_z = baf.loc[z]
        rrdr_z = rrdr.loc[z]

        for is_wgd, cn_list in [(False, cn_nowgd_all), (True, cn_wgd_all)]:
            mc = maxcn_wgd if is_wgd else maxcn
            for a_orig, b_orig in cn_list:
                a, b = (b_orig, a_orig) if is_major else (a_orig, b_orig)

                rrd_degenerate = (a + b == 4) if is_wgd else (a + b == 2)
                purities = np.empty(len(samples))
                valid = True
                for si, s in enumerate(samples):
                    baf_val = baf_z[s]
                    dom_baf = (b - 1) - baf_val * (a + b - 2)
                    pbaf = -1 if dom_baf == 0 else (2 * baf_val - 1) / dom_baf
                    if pbaf <= 0.0 or pbaf > 1.0:
                        valid = False
                        break
                    if rrd_degenerate:
                        purities[si] = pbaf
                    else:
                        rrd_val = rrdr_z[s]
                        num_rrd = 2 * rrd_val - 2
                        if is_wgd:
                            dom_rrd = a + b - 2 - 2 * rrd_val
                        else:
                            dom_rrd = a + b - 2
                        prrd = -1 if dom_rrd == 0 else num_rrd / dom_rrd
                        if prrd <= 0.0 or prrd > 1.0 or abs(pbaf - prrd) > tolerance:
                            valid = False
                            break
                        purities[si] = (pbaf + prrd) / 2
                if not valid:
                    continue

                if is_wgd:
                    gammas_arr = np.empty(len(samples))
                    gamma_valid = True
                    for si, s in enumerate(samples):
                        cz = a + b
                        dom = (cz - 2) * rdr.loc[s0, s] - 2 * rdr.loc[z, s]
                        if dom == 0.0:
                            gamma_valid = False
                            break
                        gammas_arr[si] = (2 * cz - 8) / dom
                    if not gamma_valid or np.any(gammas_arr <= 0):
                        continue
                else:
                    gammas_arr = gammas_nowgd_arr

                score = _score_pair(
                    purities,
                    gammas_arr,
                    clusters,
                    samples,
                    baf,
                    baf_tau,
                    rdr,
                    rd_std,
                    nbins_total,
                    mc,
                    tol_nstd,
                )
                if is_wgd:
                    valid_wgd[(z, a, b)] = (score, purities, gammas_arr)
                else:
                    valid_nowgd[(z, a, b)] = (score, purities)

    pair_nowgd, purities_nowgd = None, None
    best_score = -1
    for (z, a, b), (score, purs) in valid_nowgd.items():
        logging.debug(f"  noWGD z={z} cn=({a},{b}) score={score}")
        if score > best_score:
            best_score = score
            pair_nowgd = (s0, z, (1, 1), (a, b))
            purities_nowgd = dict(zip(samples, purs))
    if pair_nowgd is not None:
        logging.info(
            f"best noWGD pair: z={pair_nowgd[1]} cn={pair_nowgd[3]} score={best_score}"
        )

    pair_wgd, purities_wgd, gammas_wgd = None, None, None
    best_score = -1
    for (z, a, b), (score, purs, gams) in valid_wgd.items():
        logging.debug(f"  WGD z={z} cn=({a},{b}) score={score}")
        if score > best_score:
            best_score = score
            pair_wgd = (s0, z, (2, 2), (a, b))
            purities_wgd = dict(zip(samples, purs))
            gammas_wgd = dict(zip(samples, gams))
    if pair_wgd is not None:
        logging.info(
            f"best WGD pair: z={pair_wgd[1]} cn={pair_wgd[3]} score={best_score}"
        )

    return (
        s0,
        pair_nowgd,
        gammas_noWGD,
        purities_nowgd,
        pair_wgd,
        gammas_wgd,
        purities_wgd,
    )
