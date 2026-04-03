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
    """Infer RDR scaling factors (gamma), tumor purities, and clonal CN dicts.

    Returns (clonal_dip, gammas_noWGD, purities_noWGD,
             clonal_tet, gammas_WGD, purities_WGD, balanced_s).
    """
    logging.info("Infer scaling factors & tumor purity")

    # ----------------------------------------------------------------
    # pivot seg-level data
    # ----------------------------------------------------------------
    rdr = segs.pivot(index="#ID", columns="SAMPLE", values="RD")
    baf = segs.pivot(index="#ID", columns="SAMPLE", values="BAF")
    baf_se = segs.pivot(index="#ID", columns="SAMPLE", values="BAF-se")
    rd_var = segs.pivot(index="#ID", columns="SAMPLE", values="RD-var")
    nbins = segs.pivot(index="#ID", columns="SAMPLE", values="#BINS")
    clusters = baf.index.tolist()

    # pre-group per-bin RDR by (cluster, sample) for MSE scoring
    bin_rdrs = {}
    for cid in clusters:
        bin_rdrs[cid] = {}
        for s in samples:
            bc = bbcs[(bbcs["CLUSTER"] == cid) & (bbcs["SAMPLE"] == s)]
            bin_rdrs[cid][s] = bc["RD"].values.astype(np.float64)

    # ----------------------------------------------------------------
    # identify balanced (s0) and imbalanced clusters
    # ----------------------------------------------------------------
    user_balanced = {cid for cid, cn in fix_cn_dip.items() if cn == (1, 1)}
    seg_balanced = segs.drop_duplicates("#ID").set_index("#ID")["is_balanced"]

    balanced_s, imbalanced_z = [], []
    for cid in clusters:
        if cid in user_balanced or seg_balanced.loc[cid]:
            balanced_s.append(cid)
        else:
            imbalanced_z.append(cid)

    if not balanced_s:
        logging.error(
            "no balanced clusters found. Use --fix_cn_dip to specify a (1,1) cluster."
        )
        sys.exit(1)

    if len(user_balanced) > 0:
        s0 = max(user_balanced, key=lambda c: nbins.loc[c, samples].sum())
    else:
        s0 = max(balanced_s, key=lambda c: nbins.loc[c, samples].sum())
    logging.info(f"balanced={balanced_s}, s0={s0}")

    gammas_noWGD = {s: 2.0 / rdr.loc[s0, s] for s in samples}
    logging.info(f"gamma noWGD: {gammas_noWGD}")

    def _purity_from_baf(z, a, b):
        """Per-sample purity from cluster z BAF at CN (a,b). None if invalid."""
        baf_z = baf.loc[z]
        purs = np.empty(len(samples))
        for si, s in enumerate(samples):
            dom = (b - 1) - baf_z[s] * (a + b - 2)
            p = -1.0 if dom == 0 else (2 * baf_z[s] - 1) / dom
            if p <= 0.0 or p > 1.0:
                return None
            purs[si] = p
        return purs

    def _wgd_gamma(z, a, b, purs):
        """Per-sample WGD gamma. None if invalid."""
        gams = np.empty(len(samples))
        if (a + b) == 4:
            for si, s in enumerate(samples):
                gams[si] = (2 + 2 * purs[si]) / rdr.loc[s0, s]
            return gams
        for si, s in enumerate(samples):
            dom = (a + b - 2) * rdr.loc[s0, s] - 2 * rdr.loc[z, s]
            if dom == 0.0:
                return None
            gams[si] = (2 * (a + b) - 8) / dom
        return gams

    def _build_clonal(s0_cn, z, z_cn, fix_cn):
        d = {s0: s0_cn, z: z_cn}
        d.update(fix_cn)
        return d

    def _cn_candidates(maxcn):
        return [(1, 0)] + [
            (c - b, b)
            for c in range(2, maxcn + 1)
            for b in range(c // 2 + 1)
            if c - b != b
        ]

    def _is_concordant(z, a, b, purs, gammas_arr, base_ploidy):
        """Check BAF-derived and RDR-derived purities agree."""
        if (a + b) == base_ploidy or (a + b) == 2:
            return True
        for si, s in enumerate(samples):
            prdr = (gammas_arr[si] * rdr.loc[z, s] - 2) / (a + b - 2)
            if prdr <= 0.0 or prdr > 1.0:
                return False
            dom = (b - 1) - baf.loc[z, s] * (a + b - 2)
            var_p = ((b - a) / dom**2) ** 2 * baf_se.loc[z, s] ** 2
            var_r = (gammas_arr[si] / (a + b - 2)) ** 2 * rd_var.loc[z, s]
            se = np.sqrt(var_p + var_r)
            if se > 0 and 2 * norm.sf(abs(purs[si] - prdr) / se) < 0.05:
                return False
        return True

    def _fits_grid(purs, gammas_arr, maxcn_grid):
        """Check all clusters fit within CN grid bounds."""
        for si, s in enumerate(samples):
            p = purs[si]
            d = 2 * (1 - p) + maxcn_grid * p
            lo_baf, hi_baf = (1 - p) / d, (1 - p + maxcn_grid * p) / d
            lo_rdr, hi_rdr = 2 * (1 - p) / gammas_arr[si], d / gammas_arr[si]
            for c in clusters:
                if not (
                    lo_baf <= baf.loc[c, s] <= hi_baf
                    and lo_rdr <= rdr.loc[c, s] <= hi_rdr
                ):
                    return False
        return True

    def _rdr_mse(z, a, b, purs, gammas_arr):
        """Mean squared RDR residual across all bins and samples."""
        ss, n = 0.0, 0
        for si, s in enumerate(samples):
            denom = 2 * (1 - purs[si]) + (a + b) * purs[si]
            if denom <= 0:
                return None
            exp_rdr = denom / gammas_arr[si]
            rdrs = bin_rdrs[z][s]
            ss += np.sum((rdrs - exp_rdr) ** 2)
            n += len(rdrs)
        return ss / n if n >= 2 else None

    def _score_candidates(cn_list, base_ploidy, maxcn_grid, is_wgd, gammas_default):
        """Find best (z, cn, purs, gammas) by lowest RDR MSE."""
        rdr_s0 = rdr.loc[s0, :]
        results = {}
        for z in imbalanced_z:
            baf_z = baf.loc[z, :]
            rdr_z = rdr.loc[z, :]
            if not (np.all(baf_z <= 0.5) or np.all(baf_z >= 0.5)):
                continue
            if not (np.all(rdr_z <= rdr_s0) or np.all(rdr_z >= rdr_s0)):
                continue

            for a_orig, b_orig in cn_list:
                a, b = (b_orig, a_orig) if np.all(baf_z > 0.5) else (a_orig, b_orig)
                c = a + b
                if np.all(rdr_z >= rdr_s0) and c < base_ploidy:
                    continue
                if np.all(rdr_z <= rdr_s0) and c > base_ploidy:
                    continue

                purs = _purity_from_baf(z, a, b)
                if purs is None:
                    continue
                gammas_arr = _wgd_gamma(z, a, b, purs) if is_wgd else gammas_default
                if gammas_arr is None:
                    continue
                if not _is_concordant(z, a, b, purs, gammas_arr, base_ploidy):
                    continue
                if not _fits_grid(purs, gammas_arr, maxcn_grid):
                    continue

                mse = _rdr_mse(z, a, b, purs, gammas_arr)
                if mse is not None:
                    results[(z, a, b)] = (mse, purs, gammas_arr)

        if not results:
            return None

        # Pick best per z, then log all, then return overall best
        best_per_z = {}
        for (z, a, b), (mse, purs, gammas_arr) in results.items():
            if z not in best_per_z or mse < best_per_z[z][0]:
                best_per_z[z] = (mse, (a, b), purs, gammas_arr)

        for z in sorted(best_per_z):
            mse, cn, _, _ = best_per_z[z]
            logging.debug(f"  z={z} cn={cn} mse={mse:.6f}")

        best_z = min(best_per_z, key=lambda z: best_per_z[z][0])
        mse, cn, purs, gammas_arr = best_per_z[best_z]
        return (best_z, cn, purs, gammas_arr)

    # ----------------------------------------------------------------
    # user-specified imbalanced clusters → direct purity (skip search)
    # ----------------------------------------------------------------
    user_imb_dip = {cid: cn for cid, cn in fix_cn_dip.items() if cn != (1, 1)}
    user_imb_tet = {cid: cn for cid, cn in fix_cn_tet.items() if cn != (2, 2)}

    clonal_dip, purities_dip = None, None
    for z, (a, b) in user_imb_dip.items():
        purs = _purity_from_baf(z, a, b)
        if purs is not None:
            clonal_dip = _build_clonal((1, 1), z, (a, b), fix_cn_dip)
            purities_dip = dict(zip(samples, purs))
            logging.info(
                f"user-specified noWGD: z={z} cn=({a},{b}) purity={purities_dip}"
            )
            break

    clonal_tet, purities_tet, gammas_wgd = None, None, None
    for z, (a, b) in user_imb_tet.items():
        purs = _purity_from_baf(z, a, b)
        if purs is None:
            continue
        gams = _wgd_gamma(z, a, b, purs)
        if gams is not None:
            clonal_tet = _build_clonal((2, 2), z, (a, b), fix_cn_tet)
            purities_tet = dict(zip(samples, purs))
            gammas_wgd = dict(zip(samples, gams))
            logging.info(
                f"user-specified WGD: z={z} cn=({a},{b}) purity={purities_tet}"
            )
            break

    # early exit if nothing to search
    if not imbalanced_z or (clonal_dip is not None and clonal_tet is not None):
        if clonal_dip is None:
            clonal_dip = {s0: (1, 1)}
            clonal_dip.update(fix_cn_dip)
        return (
            clonal_dip,
            gammas_noWGD,
            purities_dip,
            clonal_tet,
            gammas_wgd,
            purities_tet,
            balanced_s,
        )

    # ----------------------------------------------------------------
    # search for best clonal pair
    # ----------------------------------------------------------------
    gammas_nowgd_arr = np.array([gammas_noWGD[s] for s in samples])

    if clonal_dip is None:
        result = _score_candidates(
            _cn_candidates(maxcn_z), 2, maxcn, False, gammas_nowgd_arr
        )
        if result is not None:
            z, cn, purs, _ = result
            clonal_dip = _build_clonal((1, 1), z, cn, fix_cn_dip)
            purities_dip = dict(zip(samples, purs))
            logging.info(f"best noWGD: z={z} cn={cn}")
        else:
            clonal_dip = {s0: (1, 1)}
            clonal_dip.update(fix_cn_dip)
            purities_dip = {s: 0.0 for s in samples}
            logging.warning("no valid noWGD pair, using s0 only")

    if clonal_tet is None:
        result = _score_candidates(
            _cn_candidates(maxcn_wgd_z), 4, maxcn_wgd, True, None
        )
        if result is not None:
            z, cn, purs, gams = result
            clonal_tet = _build_clonal((2, 2), z, cn, fix_cn_tet)
            purities_tet = dict(zip(samples, purs))
            gammas_wgd = dict(zip(samples, gams))
            logging.info(f"best WGD: z={z} cn={cn}")

    logging.info(f"clonal_dip={clonal_dip}")
    if clonal_tet is not None:
        logging.info(f"clonal_tet={clonal_tet}")

    return (
        clonal_dip,
        gammas_noWGD,
        purities_dip,
        clonal_tet,
        gammas_wgd,
        purities_tet,
        balanced_s,
    )
