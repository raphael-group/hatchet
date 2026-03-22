import math
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

    Steps:
      1. Classify clusters as balanced (BAF ≈ 0.5) or imbalanced via TOST
         equivalence test. Select base cluster s0 = largest balanced cluster.
      2. Compute noWGD gamma per sample: gamma_p = 2 / RDR(s0, p).
      3. For each imbalanced cluster z and candidate CN (a, b), estimate purity
         from BAF and RRD. When a+b equals base ploidy (2 noWGD, 4 WGD), RRD is
         degenerate so BAF-only purity is used. For WGD, also compute gamma from
         the (s0, z) RDR pair.
      4. Score each valid (z, a, b) pair: for each cluster j, find the best-fit
         CN state minimizing D² (sum of z-score² over samples for BAF and RDR),
         then weight by chi2.sf(D², df=2S) * nbins(j). Select highest-scoring pair.

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

    cn_nowgd_all = _build_cn_candidates(maxcn_z)
    cn_wgd_all = _build_cn_candidates(maxcn_wgd_z)

    nbins_total = nbins[samples].sum(axis=1)
    gammas_nowgd_arr = np.array([gammas_noWGD[s] for s in samples])

    valid_nowgd = {}
    valid_wgd = {}

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

        wgd_configs = [
            (False, maxcn, cn_nowgd_all),
            (True, maxcn_wgd, cn_wgd_all),
        ]
        for is_wgd, mc, cn_list in wgd_configs:
            for a_orig, b_orig in cn_list:
                a, b = (b_orig, a_orig) if is_major else (a_orig, b_orig)

                purities = np.empty(len(samples))
                valid = True
                for si, s in enumerate(samples):
                    baf_val = baf_z[s]
                    dom_baf = (b - 1) - baf_val * (a + b - 2)
                    pbaf = -1 if dom_baf == 0 else (2 * baf_val - 1) / dom_baf
                    if pbaf <= 0.0 or pbaf > 1.0:
                        valid = False
                        break
                    purities[si] = pbaf
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

                # Score (z, a, b) in two phases:
                # Phase 1: score anchor z under its hypothesized (a, b); if the
                #   expected-value math fails (denom <= 0) skip this candidate.
                # Phase 2: for every other imbalanced cluster j, search all CN
                #   states and pick the best-fit one.
                # score += chi2.sf(D², df=2S) * nbins(j) for each cluster.
                df = 2 * len(samples)
                score = 0.0
                matches = {}

                # 1. Score the anchor cluster z under its selected state (a, b)
                d2_z = 0.0
                z_ok = True
                for si, s in enumerate(samples):
                    num = (1 - purities[si] + b * purities[si])
                    denom = 2 * (1 - purities[si]) + (a + b) * purities[si]
                    if denom <= 0:
                        z_ok = False
                        break
                    exp_baf = num / denom
                    exp_rdr = denom / gammas_arr[si]
                    baf_std = math.sqrt(
                        exp_baf * (1 - exp_baf) / (baf_tau.loc[z, s] + 1)
                    )
                    d2_z += ((baf.loc[z, s] - exp_baf) / baf_std) ** 2
                    d2_z += ((rdr.loc[z, s] - exp_rdr) / rd_std.loc[z, s]) ** 2

                if z_ok:
                    weight_z = chi2.sf(d2_z, df)
                    score += weight_z * nbins_total[z]
                    matches[z] = ((a, b), int(nbins_total[z]), d2_z, weight_z)
                else:
                    continue  # anchor math failed — skip this candidate entirely

                # 2. Now loop over the rest of the imbalanced clusters
                for j in imbalanced_z:
                    if j == z:
                        continue  # already scored above
                    best_d2 = float("inf")
                    best_cn = None
                    cn_iter = [
                        (c - bb, bb)
                        for c in range(mc + 1)
                        for bb in range(c + 1)
                        if c - bb != bb
                    ]
                    for aa, bb in cn_iter:
                        d2 = 0.0
                        ok = True
                        for si, s in enumerate(samples):
                            num = (1 - purities[si] + bb * purities[si])
                            denom = 2 * (1 - purities[si]) + (aa + bb) * purities[si]
                            if denom <= 0:
                                ok = False
                                break
                            exp_baf = num / denom
                            exp_rdr = denom / gammas_arr[si]
                            baf_std = math.sqrt(
                                exp_baf * (1 - exp_baf) / (baf_tau.loc[j, s] + 1)
                            )
                            d2 += ((baf.loc[j, s] - exp_baf) / baf_std) ** 2
                            d2 += ((rdr.loc[j, s] - exp_rdr) / rd_std.loc[j, s]) ** 2
                        if ok and d2 < best_d2:
                            best_d2 = d2
                            best_cn = (aa, bb)
                    if best_cn is not None:
                        weight = chi2.sf(best_d2, df)
                        score += weight * nbins_total[j]
                        matches[j] = (best_cn, int(nbins_total[j]), best_d2, weight)

                wgd_tag = "WGD" if is_wgd else "noWGD"
                assign = " ".join(
                    f"({j},({cn[0]},{cn[1]}))" for j, (cn, nb, d2, w) in matches.items()
                )
                logging.debug(
                    f"  {wgd_tag} z={z} cn=({a},{b}) purity={purities} score={score:.2f} {assign}"
                )
                if is_wgd:
                    valid_wgd[(z, a, b)] = (score, purities, gammas_arr)
                else:
                    valid_nowgd[(z, a, b)] = (score, purities)

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
    )
