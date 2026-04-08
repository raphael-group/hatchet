import itertools
import numpy as np
from scipy.stats import beta


def is_explained_mut(nref, nmut, est_vaf, gamma=0.05):
    """Check if predicted VAF falls within a binomial CI of observed counts.

    Uses a Beta(nref+1, nmut+1) posterior to build a (1-gamma) credible
    interval and tests whether ``est_vaf`` (or its complement) lies inside.
    """
    if est_vaf <= 1e-8:
        return False
    lo, hi = beta.ppf([gamma / 2, 1 - gamma / 2], nref + 1, nmut + 1)
    return (lo <= est_vaf <= hi) or (lo <= 1.0 - est_vaf <= hi)


def relative_error(est_vaf, obs_vaf):
    if obs_vaf == 0:
        return np.inf
    return abs(est_vaf - obs_vaf) / obs_vaf


def estimate_vaf(obs_vaf, clones, cns, props, mode="range"):
    """Estimate expected VAF for a somatic SNV given CN states and clone proportions.

    For each allele (A and B), enumerate possible mutated-copy assignments
    across clones and pick the one closest to observed VAF.

    Parameters
    ----------
    obs_vaf : float
        Observed variant allele frequency.
    clones : list of str
        Clone names (e.g., ["normal", "clone1", "clone2"]).
    cns : list of str
        Per-clone CN states as "a|b" strings (same order as clones).
    props : array-like
        Clone proportions (same order as clones, sum to 1).
    mode : str
        "range" — mutated copies in [0, allele_cn] per clone.
        "all_or_none" — mutated copies are 0 or allele_cn.

    Returns
    -------
    (best_assignment, predicted_vaf, ccf, allele) or (None, nan, nan, None) if infeasible.
    """
    props = np.asarray(props, dtype=float)
    cna = np.array([int(cn.split("|")[0]) for cn in cns])
    cnb = np.array([int(cn.split("|")[1]) for cn in cns])

    baf_denom = float(np.sum(props * (cna + cnb)))
    purity = float(np.sum(props[1:]))
    if not np.isfinite(obs_vaf) or baf_denom <= 0:
        return None, np.nan, np.nan, None

    def _estimate_allele(cap_cn):
        if mode == "range":
            avail = [[0]] + [range(cap_cn[i] + 1) for i in range(1, len(clones))]
        else:
            avail = [[0]] + [[0, cap_cn[i]] for i in range(1, len(clones))]

        def baf_num(y):
            return sum(float(e) * props[i] for i, e in enumerate(y))

        ests = {x: baf_num(x) / baf_denom for x in itertools.product(*avail)}
        ests = {k: v for k, v in ests.items() if 0.0 <= v <= 1.0}
        best = min(ests, key=lambda x: abs(ests[x] - obs_vaf))
        ccf = (
            float(np.sum(props * (np.array(best) > 0))) / purity if purity > 0 else 0.0
        )
        return best, ests[best], ccf

    bestA, vafA, ccfA = _estimate_allele(cna)
    bestB, vafB, ccfB = _estimate_allele(cnb)
    if abs(vafA - obs_vaf) < abs(vafB - obs_vaf):
        return bestA, vafA, ccfA, "A"
    return bestB, vafB, ccfB, "B"
