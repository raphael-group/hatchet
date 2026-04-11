import sys
import logging
import numpy as np
import pandas as pd
from scipy.stats import norm


def get_scaling_factor(
    samples: list,
    segs: pd.DataFrame,
    bbcs: pd.DataFrame,
    fix_cn_dip: dict,
    fix_cn_tet: dict,
    maxcn: int,
    maxcn_wgd: int,
    maxcn_z: int = 4,
    maxcn_wgd_z: int = 6,
):
    """Infer RDR scaling (gamma), tumor purity, and clonal CN anchors.

    Steps:
      1. Pick s0 = largest balanced cluster (user pins first). gamma_noWGD = 2/RD(s0).
      2. noWGD: try user-pinned fix_cn_dip → search imbalanced z by RDR MSE
         → fallback to s0-only with purity 0.
      3. WGD: try user-pinned fix_cn_tet → search imbalanced z → fallback to
         another balanced s1 interpreted as (1,1), giving gamma=2/RD(s1),
         p=RD(s0)/RD(s1)-1.

    Returns (clonal_dip, gammas_noWGD, purities_noWGD,
             clonal_tet, gammas_WGD, purities_WGD, balanced_s).
    """
    logging.info("Infer scaling factors & tumor purity")

    def _pivot(col):
        return segs.pivot(index="#ID", columns="SAMPLE", values=col)[samples]

    rdr, baf = _pivot("RD"), _pivot("BAF")
    baf_se, rd_var = _pivot("BAF-se"), _pivot("RD-var")
    nbins = _pivot("#BINS")
    clusters = baf.index.tolist()
    baf_mat, rdr_mat = baf.values, rdr.values

    bin_rdrs = {
        key: g["RD"].values.astype(np.float64)
        for key, g in bbcs.groupby(["CLUSTER", "SAMPLE"])
    }

    user_balanced = {c for c, cn in fix_cn_dip.items() if cn == (1, 1)} | {
        c for c, cn in fix_cn_tet.items() if cn == (2, 2)
    }
    seg_balanced = segs.drop_duplicates("#ID").set_index("#ID")["is_balanced"]
    balanced_s = [c for c in clusters if c in user_balanced or seg_balanced.loc[c]]
    imbalanced_z = [c for c in clusters if c not in balanced_s]
    if not balanced_s:
        logging.error(
            "no balanced clusters found. Use --fix_cn_dip/--fix_cn_tet to specify a baseline balanced cluster."
        )
        sys.exit(1)
    balanced_s.sort(
        key=lambda c: (c in user_balanced, nbins.loc[c].sum()),
        reverse=True,
    )
    s0 = balanced_s[0]
    rdr_s0 = rdr.loc[s0].values
    gammas_noWGD = dict(zip(samples, 2.0 / rdr_s0))
    logging.info(f"balanced={balanced_s}, s0={s0}")
    logging.info(f"gamma noWGD: {gammas_noWGD}")

    def _purity_from_baf(z, a, b):
        bz = baf.loc[z].values
        dom = (b - 1) - bz * (a + b - 2)
        with np.errstate(divide="ignore", invalid="ignore"):
            p = np.where(dom == 0, -1.0, (2 * bz - 1) / dom)
        return None if np.any((p <= 0) | (p > 1)) else p

    def _wgd_gamma(z, a, b, purs):
        c = a + b
        if c == 4:
            return (2 + 2 * purs) / rdr_s0
        dom = (c - 2) * rdr_s0 - 2 * rdr.loc[z].values
        return None if np.any(dom == 0) else (2 * c - 8) / dom

    def _cn_candidates(mx):
        return [(1, 0)] + [
            (c - b, b)
            for c in range(2, mx + 1)
            for b in range(c // 2 + 1)
            if c - b != b
        ]

    def _is_concordant(z, a, b, purs, gams, base_ploidy):
        c = a + b
        if c == base_ploidy or c == 2:
            return True
        prdr = (gams * rdr.loc[z].values - 2) / (c - 2)
        if np.any((prdr <= 0) | (prdr > 1)):
            return False
        dom = (b - 1) - baf.loc[z].values * (c - 2)
        var_p = ((b - a) / dom**2) ** 2 * baf_se.loc[z].values ** 2
        var_r = (gams / (c - 2)) ** 2 * rd_var.loc[z].values
        se = np.sqrt(var_p + var_r)
        with np.errstate(divide="ignore", invalid="ignore"):
            pvals = 2 * norm.sf(np.abs(purs - prdr) / se)
        return not np.any((se > 0) & (pvals < 0.05))

    def _fits_grid(purs, gams, mx):
        d = 2 * (1 - purs) + mx * purs
        lo_baf, hi_baf = (1 - purs) / d, (1 - purs + mx * purs) / d
        lo_rdr, hi_rdr = 2 * (1 - purs) / gams, d / gams
        return bool(
            ((baf_mat >= lo_baf) & (baf_mat <= hi_baf)).all()
            and ((rdr_mat >= lo_rdr) & (rdr_mat <= hi_rdr)).all()
        )

    def _rdr_mse(z, a, b, purs, gams):
        exp_rdr = (2 * (1 - purs) + (a + b) * purs) / gams
        ss, n = 0.0, 0
        for si, s in enumerate(samples):
            rdrs = bin_rdrs.get((z, s), np.empty(0))
            ss += np.sum((rdrs - exp_rdr[si]) ** 2)
            n += len(rdrs)
        return ss / n if n >= 2 else None

    def _score_candidates(cn_list, base_ploidy, maxcn_grid, is_wgd, gams_default):
        best = {}
        for z in imbalanced_z:
            baf_z, rdr_z = baf.loc[z].values, rdr.loc[z].values
            if not ((baf_z <= 0.5).all() or (baf_z >= 0.5).all()):
                continue
            if not ((rdr_z <= rdr_s0).all() or (rdr_z >= rdr_s0).all()):
                continue
            gain = bool((rdr_z >= rdr_s0).all())
            flip = bool((baf_z > 0.5).all())
            for a0, b0 in cn_list:
                a, b = (b0, a0) if flip else (a0, b0)
                c = a + b
                if (gain and c < base_ploidy) or (not gain and c > base_ploidy):
                    continue
                purs = _purity_from_baf(z, a, b)
                if purs is None:
                    continue
                gams = _wgd_gamma(z, a, b, purs) if is_wgd else gams_default
                if gams is None:
                    continue
                if not _is_concordant(z, a, b, purs, gams, base_ploidy):
                    continue
                if not _fits_grid(purs, gams, maxcn_grid):
                    continue
                mse = _rdr_mse(z, a, b, purs, gams)
                if mse is None:
                    continue
                if z not in best or mse < best[z][0]:
                    best[z] = (mse, (a, b), purs, gams)
        if not best:
            return None
        for z in sorted(best):
            logging.debug(f"  z={z} cn={best[z][1]} mse={best[z][0]:.6f}")
        z = min(best, key=lambda k: best[k][0])
        _, cn, purs, gams = best[z]
        return z, cn, purs, gams

    clonal_dip, purities_dip = None, None
    for z, cn in fix_cn_dip.items():
        if cn == (1, 1):
            continue
        purs = _purity_from_baf(z, *cn)
        if purs is None:
            continue
        clonal_dip = {s0: (1, 1), z: cn, **fix_cn_dip}
        purities_dip = dict(zip(samples, purs))
        logging.info(f"user-specified noWGD: z={z} cn={cn} purity={purities_dip}")
        break

    if purities_dip is None and imbalanced_z:
        result = _score_candidates(
            _cn_candidates(maxcn_z), 2, maxcn, False, 2.0 / rdr_s0
        )
        if result is not None:
            z, cn, purs, _ = result
            clonal_dip = {s0: (1, 1), z: cn, **fix_cn_dip}
            purities_dip = dict(zip(samples, purs))
            logging.info(f"best noWGD: z={z} cn={cn}")

    if purities_dip is None:
        clonal_dip = {s0: (1, 1), **fix_cn_dip}
        purities_dip = {s: 0.0 for s in samples}
        logging.warning("no valid noWGD pair, using s0 only")

    clonal_tet, purities_tet, gammas_wgd = None, None, None
    for z, cn in fix_cn_tet.items():
        if cn == (2, 2):
            continue
        purs = _purity_from_baf(z, *cn)
        if purs is None:
            continue
        gams = _wgd_gamma(z, *cn, purs)
        if gams is None:
            continue
        clonal_tet = {s0: (2, 2), z: cn, **fix_cn_tet}
        purities_tet = dict(zip(samples, purs))
        gammas_wgd = dict(zip(samples, gams))
        logging.info(f"user-specified WGD: z={z} cn={cn} purity={purities_tet}")
        break

    if purities_tet is None and imbalanced_z:
        result = _score_candidates(
            _cn_candidates(maxcn_wgd_z), 4, maxcn_wgd, True, None
        )
        if result is not None:
            z, cn, purs, gams = result
            clonal_tet = {s0: (2, 2), z: cn, **fix_cn_tet}
            purities_tet = dict(zip(samples, purs))
            gammas_wgd = dict(zip(samples, gams))
            logging.info(f"best WGD: z={z} cn={cn}")

    if purities_tet is None:
        for s1 in balanced_s[1:]:
            if fix_cn_tet.get(s1, (1, 1)) != (1, 1):
                continue
            r1 = rdr.loc[s1].values
            if np.any(r1 <= 0):
                continue
            purs = rdr_s0 / r1 - 1
            if np.any((purs <= 0) | (purs > 1)):
                continue
            gams = 2.0 / r1
            if not _fits_grid(purs, gams, maxcn_wgd):
                continue
            clonal_tet = {s0: (2, 2), s1: (1, 1), **fix_cn_tet}
            purities_tet = dict(zip(samples, purs))
            gammas_wgd = dict(zip(samples, gams))
            logging.info(f"balanced-fallback WGD: s1=c{s1}=(1,1)")
            break

    if purities_tet is None:
        logging.warning("no valid WGD pair; tetraploid solution unavailable")

    return (
        clonal_dip,
        gammas_noWGD,
        purities_dip,
        clonal_tet,
        gammas_wgd,
        purities_tet,
        balanced_s,
    )
