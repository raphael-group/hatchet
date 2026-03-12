import os
import logging

import kneed
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.special import betaln, gammaln


def _ll_gauss(obs_rdrs, exp_rdr, rdr_var, floor_var=1e-12):
    """Gaussian log-likelihood in RDR space with fixed variance."""
    valid = np.isfinite(obs_rdrs)
    x = obs_rdrs[valid]
    n = x.size
    if n == 0:
        return 0.0
    v = max(float(rdr_var), floor_var)
    return float(-0.5 * np.sum((x - exp_rdr) ** 2) / v
                 - 0.5 * n * np.log(2 * np.pi * v))


def _ll_betabinom(b_counts, total_counts, p, tau, eps=1e-9):
    """Beta-Binomial log-likelihood with fixed dispersion tau."""
    valid = np.isfinite(b_counts) & np.isfinite(total_counts) & (total_counts > 0)
    b = b_counts[valid]
    t = total_counts[valid]
    if b.size == 0:
        return 0.0
    a_counts = t - b
    p = float(np.clip(p, eps, 1.0 - eps))
    bb_alpha = tau * p
    bb_beta = tau * (1.0 - p)
    ll = (betaln(b + bb_alpha, a_counts + bb_beta) - betaln(bb_alpha, bb_beta)
          + gammaln(t + 1) - gammaln(b + 1) - gammaln(a_counts + 1))
    return float(np.sum(ll))


def _compute_loglik_from_ucn(
    ucn_file: str, n: int, gammas: dict, segs: pd.DataFrame
) -> float:
    """Compute conditioned log-likelihood for a single (ploidy, n) solution.

    Uses fixed RDR variance and BAF dispersion from the HMM (seg file) rather
    than concentrating out variance, so that a wrong exp_FCN is always penalised.

    RDR component: Gaussian with fixed per-(cluster, sample) variance from seg.
    BAF component: Beta-Binomial with fixed per-sample dispersion tau from seg.

    Args:
        ucn_file: path to the results.{ploidy}.n{n}.bbc.ucn.tsv file.
        n:        total number of clones (including normal).
        gammas:   dict {sample_id: gamma} for RDR-to-FCN scaling.
        segs:     DataFrame with columns #ID, SAMPLE, RD-var, BAF-tau.

    Returns:
        Scalar log-likelihood (higher is better fit).
    """
    bbcs = pd.read_table(ucn_file, sep="\t")

    # Build lookup for HMM-estimated variance/dispersion
    seg_lookup = {}
    for _, row in segs.iterrows():
        seg_lookup[(row["#ID"], row["SAMPLE"])] = (
            float(row["RD-var"]),
            float(row["BAF-tau"]),
        )

    ll = 0.0
    for (cluster_id, sample_id), bbc_sub in bbcs.groupby(
        ["CLUSTER", "SAMPLE"], sort=False
    ):
        exp_fcn = 2.0
        exp_fcn_b = 1.0
        if n > 1:
            row = bbc_sub.iloc[0]
            exp_fcn = float(row["u_normal"]) * 2
            exp_fcn_b = float(row["u_normal"])
            for nn in range(1, n):
                a_str, b_str = str(row[f"cn_clone{nn}"]).split("|")
                a, b = int(a_str), int(b_str)
                exp_fcn += float(row[f"u_clone{nn}"]) * (a + b)
                exp_fcn_b += float(row[f"u_clone{nn}"]) * b
        p = exp_fcn_b / exp_fcn if exp_fcn > 0 else 0.5

        # Look up HMM-estimated parameters
        rdr_var, baf_tau = seg_lookup[(cluster_id, sample_id)]

        # Gaussian RDR log-likelihood (in RDR space, fixed variance)
        obs_rdrs = bbc_sub["RD"].to_numpy().astype(float)
        exp_rdr = exp_fcn / gammas[sample_id]
        ll += _ll_gauss(obs_rdrs, exp_rdr, rdr_var)

        # Beta-Binomial BAF log-likelihood (fixed dispersion)
        b_counts = bbc_sub["BETA"].to_numpy().astype(float)
        total_counts = (bbc_sub["ALPHA"] + bbc_sub["BETA"]).to_numpy().astype(float)
        ll += _ll_betabinom(b_counts, total_counts, p, baf_tau)

    return ll


def model_selection(
    diploid_sols: dict,
    tetraploid_sols: dict,
    out_dir: str,
    gammas_noWGD: dict,
    gammas_WGD: dict,
    segs: pd.DataFrame,
    v=1,
):
    """Select best n and ploidy using conditioned Gaussian RDR + Beta-Binomial BAF log-likelihood.

    Uses fixed RDR variance and BAF dispersion from the HMM (seg file) so that
    the n=1 all-diploid baseline cannot absorb prediction error into variance.

    For each ploidy with at least one solution:
      1. Compute loglik for every n by reading the corresponding bbc.ucn file.
      2. Select best n via kneed elbow on the neg-loglik vs n curve.
         With < 3 solutions the smallest n is chosen (most parsimonious).
      3. Final ploidy chosen by parsimony (fewest clones; diploid preferred on tie).

    Args:
        diploid_sols:    {n: (obj, imf_obj)} for diploid runs.
        tetraploid_sols: {n: (obj, imf_obj)} for tetraploid runs.
        out_dir:         directory containing results.{ploidy}.n{n}.bbc.ucn.tsv files.
        gammas_noWGD:    {sample: gamma} RDR scaling for diploid solutions.
        gammas_WGD:      {sample: gamma} RDR scaling for tetraploid solutions.
        segs:            DataFrame with HMM-estimated RD-var and BAF-tau per (cluster, sample).
        v:               verbosity level.

    Returns:
        (n_dip, n_tet, final_selection) where final_selection is "diploid" or
        "tetraploid" (or None if no solutions exist for that ploidy).
    """
    if len(diploid_sols) == 0 and len(tetraploid_sols) == 0:
        logging.info(
            "ERROR! no solution found for either diploid or tetraploid setting!"
        )
        raise ValueError("final model selection error")

    def _select_by_loglik(sols: dict, ploidy: str, gammas: dict) -> int:
        ns_sorted = sorted(sols.keys())

        # Prepend fake n=1 (all-diploid baseline) using first solution's UCN file
        first_ucn = os.path.join(
            out_dir, f"results.{ploidy}.n{ns_sorted[0]}.bbc.ucn.tsv"
        )
        ll_n1 = _compute_loglik_from_ucn(first_ucn, 1, gammas, segs)
        logging.info(f"{ploidy}: n=1 (baseline), loglik={ll_n1:.4f}")

        ns_all = [1] + ns_sorted
        lls = [ll_n1]
        for n in ns_sorted:
            ucn_file = os.path.join(out_dir, f"results.{ploidy}.n{n}.bbc.ucn.tsv")
            ll = _compute_loglik_from_ucn(ucn_file, n, gammas, segs)
            lls.append(ll)
            logging.info(f"{ploidy}: n={n}, loglik={ll:.4f}")

        ns = np.array(ns_all, dtype=np.int32)
        neg_lls = -1.0 * np.array(lls)

        chosen_n = int(ns_sorted[0])  # default: fewest real clones
        if len(ns) >= 3:
            kl = kneed.KneeLocator(
                x=ns, y=neg_lls, curve="convex", direction="decreasing"
            )
            if kl.elbow is not None:
                chosen_n = max(int(ns_sorted[0]), int(kl.elbow))

        logging.info(f"{ploidy}: chosen n={chosen_n} via loglik elbow")
        return chosen_n, ns, neg_lls

    ploidy2scores = {}

    n2 = 0
    if len(diploid_sols) > 0:
        n2, ns_dip, neg_lls_dip = _select_by_loglik(
            diploid_sols, "diploid", gammas_noWGD
        )
        ploidy2scores["diploid"] = (n2, ns_dip, neg_lls_dip)
        logging.info(f"best diploid solution: n={n2}")

    n4 = 0
    if len(tetraploid_sols) > 0:
        n4, ns_tet, neg_lls_tet = _select_by_loglik(
            tetraploid_sols, "tetraploid", gammas_WGD
        )
        ploidy2scores["tetraploid"] = (n4, ns_tet, neg_lls_tet)
        logging.info(f"best tetraploid solution: n={n4}")

    # parsimony: prefer fewer clones; diploid preferred on tie
    if len(tetraploid_sols) == 0:
        final_selection = "diploid"
    elif len(diploid_sols) == 0:
        final_selection = "tetraploid"
    else:
        final_selection = "diploid" if n2 <= n4 else "tetraploid"

    logging.info(f"final selection: {final_selection}")

    # plot elbow curve
    fig, ax = plt.subplots(figsize=(4.5, 3.5))
    all_ns = []
    colors = {"diploid": "#1f77b4", "tetraploid": "#d62728"}
    for ploidy_key in ["diploid", "tetraploid"]:
        if ploidy_key not in ploidy2scores:
            continue
        chosen_n, ns, neg_lls = ploidy2scores[ploidy_key]
        c = colors[ploidy_key]
        ax.plot(ns, neg_lls, "-o", color=c, markersize=5, label=ploidy_key.capitalize())
        idx = list(ns).index(chosen_n)
        ax.plot(chosen_n, neg_lls[idx], "*", color=c, markersize=14, zorder=5)
        all_ns.extend(ns.tolist())
    ax.set_xticks(sorted(set(int(x) for x in all_ns)))
    ax.set_xlabel("Number of clones")
    ax.set_ylabel("Negative log-likelihood")
    best_n = n2 if final_selection == "diploid" else n4
    ax.set_title(f"Model selection (best: {final_selection}, n={best_n})")
    ax.legend(framealpha=0.9)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "elbow_curve.png"), dpi=150)
    plt.close(fig)

    return n2, n4, final_selection
