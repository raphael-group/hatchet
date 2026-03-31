import sys
import logging
import numpy as np
import pandas as pd
from scipy.stats import chi2, norm


def _build_cn_candidates(maxcn):
    """Build imbalanced CN candidates: (a, b) with a > b, a + b <= maxcn."""
    return [(1, 0)] + [
        (c - b, b) for c in range(2, maxcn + 1) for b in range(c // 2 + 1) if c - b != b
    ]


def _compute_chi2_pweight(
    cluster, a, b, purities, gammas, samples, baf, rdr, baf_tau, rd_std, nbins_total, df
):
    """Compute chi2 p-value weighted by cluster bin count under CN (a, b).

    Returns ``chi2.sf(D², df) * nbins_total[cluster]``, where D² is the
    sum of squared BAF and RDR z-scores across samples.  Returns None if
    the expected values are invalid (denom <= 0).
    """
    d2 = 0.0
    for si, s in enumerate(samples):
        denom = 2 * (1 - purities[si]) + (a + b) * purities[si]
        if denom <= 0:
            return None
        exp_baf = (1 - purities[si] + b * purities[si]) / denom
        exp_rdr = denom / gammas[si]
        baf_std = np.sqrt(exp_baf * (1 - exp_baf) / (baf_tau.loc[cluster, s] + 1))
        d2 += ((baf.loc[cluster, s] - exp_baf) / baf_std) ** 2
        d2 += ((rdr.loc[cluster, s] - exp_rdr) / rd_std.loc[cluster, s]) ** 2
    return chi2.sf(d2, df) * nbins_total[cluster]


def get_scaling_factor(
    samples: list,
    seg: pd.DataFrame,
    bal_tost_alpha: float,
    bal_tost_margin: float,
    maxcn: int,
    maxcn_wgd: int,
    maxcn_z: int = 4,
    maxcn_wgd_z: int = 6,
):
    """Infer RDR scaling factors (gamma) and tumor purities.

    Score each (z, a_z, b_z) by forming a clonal triple (s0, z, j): pick the
    single best-supporting cluster j. Select highest-scoring triple.

    Returns (s0, pair_noWGD, gammas_noWGD, purities_noWGD,
             pair_WGD, gammas_WGD, purities_WGD).
    """
    logging.info("Infer scaling factors & tumor purity")

    rdr = seg.pivot(index="#ID", columns="SAMPLE", values="RD")
    baf = seg.pivot(index="#ID", columns="SAMPLE", values="BAF")
    baf_se = seg.pivot(index="#ID", columns="SAMPLE", values="BAF-se")
    baf_tau = seg.pivot(index="#ID", columns="SAMPLE", values="BAF-tau")
    rd_std = np.sqrt(seg.pivot(index="#ID", columns="SAMPLE", values="RD-var"))
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
        logging.warning("no unbalanced clusters found, skip purity estimation")
        return (s0, None, gammas_noWGD, None, None, None, None, balanced_s)

    cn_nowgd_all = _build_cn_candidates(maxcn_z)
    cn_wgd_all = _build_cn_candidates(maxcn_wgd_z)
    nbins_total = nbins[samples].sum(axis=1)
    gammas_nowgd_arr = np.array([gammas_noWGD[s] for s in samples])
    rdr_s0 = rdr.loc[s0, :]
    df = 2 * len(samples)
    chi2_args = dict(
        samples=samples,
        baf=baf,
        rdr=rdr,
        baf_tau=baf_tau,
        rd_std=rd_std,
        nbins_total=nbins_total,
        df=df,
    )

    valid_nowgd = {}
    valid_wgd = {}

    for z in imbalanced_z:
        if not (
            np.all(baf.loc[z, :] < 0.5 + bal_tost_margin)
            or np.all(baf.loc[z, :] > 0.5 - bal_tost_margin)
        ):
            logging.debug(f"cluster {z} inconsistent BAF side across samples, skip")
            continue
        if not (np.all(rdr.loc[z, :] > rdr_s0) or np.all(rdr.loc[z, :] < rdr_s0)):
            logging.debug(f"cluster {z} inconsistent RDR side vs s0, skip")
            continue

        z_above_s0 = np.all(rdr.loc[z, :] > rdr_s0)
        is_major = np.all(baf.loc[z, :] > 0.5)
        baf_z = baf.loc[z]

        for is_wgd, cn_list, base_ploidy, maxcn_grid in [
            (False, cn_nowgd_all, 2, maxcn),
            (True, cn_wgd_all, 4, maxcn_wgd),
        ]:
            best_z_score = -1.0
            best_z_log = None
            for a_orig, b_orig in cn_list:
                a, b = (b_orig, a_orig) if is_major else (a_orig, b_orig)

                if z_above_s0 and (a + b) < base_ploidy:
                    continue
                if not z_above_s0 and (a + b) > base_ploidy:
                    continue

                # --- estimate purity from BAF ---
                purities_baf = np.empty(len(samples))
                valid = True
                for si, s in enumerate(samples):
                    dom = (b - 1) - baf_z[s] * (a + b - 2)
                    pbaf = -1 if dom == 0 else (2 * baf_z[s] - 1) / dom
                    if pbaf <= 0.0 or pbaf > 1.0:
                        valid = False
                        break
                    purities_baf[si] = pbaf
                if not valid:
                    continue

                # --- compute gamma ---
                if is_wgd:
                    gammas_arr = np.empty(len(samples))
                    gamma_valid = True
                    if (a + b) == base_ploidy:
                        # a + b equals WGD base ploidy — derive gamma from BAF purity
                        for si, s in enumerate(samples):
                            gammas_arr[si] = (2 + 2 * purities_baf[si]) / rdr.loc[s0, s]
                    else:
                        for si, s in enumerate(samples):
                            dom = (a + b - 2) * rdr.loc[s0, s] - 2 * rdr.loc[z, s]
                            if dom == 0.0:
                                gamma_valid = False
                                break
                            gammas_arr[si] = (2 * (a + b) - 8) / dom
                    if not gamma_valid or np.any(gammas_arr <= 0):
                        continue
                else:
                    gammas_arr = gammas_nowgd_arr

                # --- RDR purity concordance check (when a+b != 2) ---
                purities_rdr = None
                if (a + b) != 2:
                    purities_rdr = np.empty(len(samples))
                    for si, s in enumerate(samples):
                        prdr = (gammas_arr[si] * rdr.loc[z, s] - 2) / (a + b - 2)
                        if prdr <= 0.0 or prdr > 1.0:
                            valid = False
                            break
                        purities_rdr[si] = prdr

                        dom = (b - 1) - baf_z[s] * (a + b - 2)
                        var_pbaf = ((b - a) / dom**2) ** 2 * baf_se.loc[z, s] ** 2
                        var_prdr = (
                            gammas_arr[si] / (a + b - 2)
                        ) ** 2 * rd_std.loc[z, s] ** 2
                        se_diff = np.sqrt(var_pbaf + var_prdr)
                        if se_diff > 0:
                            z_stat = abs(purities_baf[si] - prdr) / se_diff
                            if 2 * norm.sf(z_stat) < bal_tost_alpha:
                                valid = False
                                break
                    if not valid:
                        continue
                purities = purities_baf

                # --- check all clusters within clonal grid bounds ---
                grid_ok = True
                for si, s in enumerate(samples):
                    p = purities[si]
                    denom_max = 2 * (1 - p) + maxcn_grid * p
                    min_baf = (1 - p) / denom_max
                    max_baf = (1 - p + maxcn_grid * p) / denom_max
                    min_rdr = 2 * (1 - p) / gammas_arr[si]
                    max_rdr = denom_max / gammas_arr[si]
                    for cid in clusters:
                        if (
                            baf.loc[cid, s] < min_baf
                            or baf.loc[cid, s] > max_baf
                            or rdr.loc[cid, s] < min_rdr
                            or rdr.loc[cid, s] > max_rdr
                        ):
                            grid_ok = False
                if not grid_ok:
                    continue

                # --- score anchor z ---
                score_z = _compute_chi2_pweight(
                    z, a, b, purities, gammas_arr, **chi2_args
                )
                if score_z is None:
                    continue

                # --- find single best-supporting cluster j ---
                best_j, best_j_cn, best_j_score = None, None, 0.0
                for j in imbalanced_z:
                    if j == z:
                        continue
                    j_above = np.all(rdr.loc[j, :] > rdr_s0)
                    j_below = np.all(rdr.loc[j, :] < rdr_s0)
                    for aa, bb in cn_list:
                        if j_above and (aa + bb) < base_ploidy:
                            continue
                        if j_below and (aa + bb) > base_ploidy:
                            continue
                        s_j = _compute_chi2_pweight(
                            j, aa, bb, purities, gammas_arr, **chi2_args
                        )
                        if s_j is not None:
                            if s_j > best_j_score:
                                best_j_score = s_j
                                best_j = j
                                best_j_cn = (aa, bb)

                score = score_z + best_j_score

                wgd_tag = "WGD" if is_wgd else "noWGD"
                if score > best_z_score:
                    best_z_score = score
                    j_info = (
                        f"j={best_j} jcn={best_j_cn}"
                        if best_j is not None
                        else "j=None"
                    )
                    rdr_info = (
                        f"pRDR={purities_rdr}"
                        if purities_rdr is not None
                        else "pRDR=N/A"
                    )
                    best_z_log = (
                        f"  {wgd_tag} z={z} cn=({a},{b}) pBAF={purities_baf} {rdr_info} "
                        f"score={score:.2f} {j_info}"
                    )
                if is_wgd:
                    valid_wgd[(z, a, b)] = (score, purities, gammas_arr)
                else:
                    valid_nowgd[(z, a, b)] = (score, purities)

            if best_z_log is not None:
                logging.debug(best_z_log)

    pair_nowgd, purities_nowgd = None, None
    best_score = -1.0
    for (z, a, b), (score, purs) in valid_nowgd.items():
        if score > best_score:
            best_score = score
            pair_nowgd = (s0, z, (1, 1), (a, b))
            purities_nowgd = dict(zip(samples, purs))
    if pair_nowgd is not None:
        logging.info(
            f"best noWGD pair: z={pair_nowgd[1]} cn={pair_nowgd[3]} score={best_score:.2f}"
        )

    pair_wgd, purities_wgd, gammas_wgd = None, None, None
    best_score = -1.0
    for (z, a, b), (score, purs, gams) in valid_wgd.items():
        if score > best_score:
            best_score = score
            pair_wgd = (s0, z, (2, 2), (a, b))
            purities_wgd = dict(zip(samples, purs))
            gammas_wgd = dict(zip(samples, gams))
    if pair_wgd is not None:
        logging.info(
            f"best WGD pair: z={pair_wgd[1]} cn={pair_wgd[3]} score={best_score:.2f}"
        )

    return (
        s0,
        pair_nowgd,
        gammas_noWGD,
        purities_nowgd,
        pair_wgd,
        gammas_wgd,
        purities_wgd,
        balanced_s,
    )
