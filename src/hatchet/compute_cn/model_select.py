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
    return float(
        -0.5 * np.sum((x - exp_rdr) ** 2) / v - 0.5 * n * np.log(2 * np.pi * v)
    )


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
    ll = (
        betaln(b + bb_alpha, a_counts + bb_beta)
        - betaln(bb_alpha, bb_beta)
        + gammaln(t + 1)
        - gammaln(b + 1)
        - gammaln(a_counts + 1)
    )
    return float(np.sum(ll))


def _compute_loglik_from_ucn(ucn_file: str, n: int, gammas: dict, segs: pd.DataFrame):
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
        (ll, n_obs, n_clusters, n_samples): log-likelihood, observation count,
        number of clusters, and number of samples.
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
    n_obs = 0
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
        valid_rdr = np.isfinite(obs_rdrs)
        exp_rdr = exp_fcn / gammas[sample_id]
        ll += _ll_gauss(obs_rdrs, exp_rdr, rdr_var)
        n_obs += int(valid_rdr.sum())

        # Beta-Binomial BAF log-likelihood (fixed dispersion)
        b_counts = bbc_sub["BETA"].to_numpy().astype(float)
        total_counts = (bbc_sub["ALPHA"] + bbc_sub["BETA"]).to_numpy().astype(float)
        ll += _ll_betabinom(b_counts, total_counts, p, baf_tau)
        valid_baf = (
            np.isfinite(b_counts) & np.isfinite(total_counts) & (total_counts > 0)
        )
        n_obs += int(valid_baf.sum())

    n_clusters = bbcs["CLUSTER"].nunique()
    n_samples = bbcs["SAMPLE"].nunique()
    return ll, n_obs, n_clusters, n_samples


def _count_params(n, n_clusters, n_samples):
    """Count free parameters for the CN deconvolution model.

    Integer CN: 2 * n_clusters * (n-1)  (allele-specific per tumor clone)
    Proportions: (n-1) * n_samples      (sum-to-1 constraint removes 1 per sample)
    """
    if n <= 1:
        return 0
    n_tumor = n - 1
    return 2 * n_clusters * n_tumor + n_tumor * n_samples


def model_selection(
    diploid_sols: dict,
    tetraploid_sols: dict,
    out_dir: str,
    gammas_noWGD: dict,
    gammas_WGD: dict,
    segs: pd.DataFrame,
    method="elbow",
    plot_dir=None,
):
    """Select best n and ploidy.

    Two methods are supported (controlled by *method*):

    ``"elbow"`` (default)
        Compute conditioned log-likelihood for every n, select the elbow/knee
        on the neg-loglik curve (kneed).  Ploidy chosen by parsimony (fewer
        clones; diploid preferred on tie).

    ``"bic"``
        Compute BIC = -2·LL + k·ln(N) for every n, select the n with
        minimum BIC.  Ploidy chosen by minimum BIC across ploidies.

    Returns:
        (n_dip, n_tet, final_selection) where final_selection is "diploid" or
        "tetraploid" (or None if no solutions exist for that ploidy).
    """
    if len(diploid_sols) == 0 and len(tetraploid_sols) == 0:
        logging.info(
            "ERROR! no solution found for either diploid or tetraploid setting!"
        )
        raise ValueError("final model selection error")

    def _compute_scores(sols: dict, ploidy: str, gammas: dict):
        """Compute loglik (and BIC) for all n values including n=1 baseline.

        n_clusters and n_samples are constant across all n for a given ploidy
        (same clusters/samples in every UCN file).  The loop overwrites these
        variables on each iteration; using the final values is correct.
        """
        ns_sorted = sorted(sols.keys())
        first_ucn = os.path.join(
            out_dir, f"results.{ploidy}.n{ns_sorted[0]}.bbc.ucn.tsv"
        )
        ll_n1, nobs_n1, n_clusters, n_samples = _compute_loglik_from_ucn(
            first_ucn, 1, gammas, segs
        )
        logging.info(f"{ploidy}: n=1 (baseline), loglik={ll_n1:.4f}")

        ns_all = [1] + ns_sorted
        lls = [ll_n1]
        n_obs_list = [nobs_n1]
        for clone_n in ns_sorted:
            ucn_file = os.path.join(out_dir, f"results.{ploidy}.n{clone_n}.bbc.ucn.tsv")
            ll, nobs, n_clusters, n_samples = _compute_loglik_from_ucn(
                ucn_file, clone_n, gammas, segs
            )
            lls.append(ll)
            n_obs_list.append(nobs)
            logging.info(f"{ploidy}: n={clone_n}, loglik={ll:.4f}")

        ns = np.array(ns_all, dtype=np.int32)
        lls = np.array(lls)
        neg_lls = -lls
        n_obs = n_obs_list[0]  # same across all n for a given ploidy

        bics = np.array(
            [
                -2 * lls[i] + _count_params(ns[i], n_clusters, n_samples) * np.log(n_obs)
                for i in range(len(ns))
            ]
        )

        return ns, lls, neg_lls, bics, ns_sorted

    def _select_elbow(ns, neg_lls, ns_sorted, ploidy):
        chosen_n = int(ns_sorted[0])
        if len(ns) >= 3:
            kl = kneed.KneeLocator(
                x=ns, y=neg_lls, curve="convex", direction="decreasing"
            )
            if kl.elbow is not None:
                chosen_n = max(int(ns_sorted[0]), int(kl.elbow))
        logging.info(f"{ploidy}: chosen n={chosen_n} via loglik elbow")
        return chosen_n

    def _select_bic(ns, bics, ns_sorted, ploidy):
        best_idx = int(np.argmin(bics))
        chosen_n = int(ns[best_idx])
        # n=1 is synthetic; if BIC picks it, fall back to smallest real n
        if chosen_n < ns_sorted[0]:
            chosen_n = int(ns_sorted[0])
        for i, clone_n in enumerate(ns):
            logging.info(f"{ploidy}: n={clone_n}, BIC={bics[i]:.2f}")
        logging.info(f"{ploidy}: chosen n={chosen_n} via BIC")
        return chosen_n

    ploidy2scores = {}

    n_dip = 0
    bic_dip = np.inf
    if len(diploid_sols) > 0:
        ns, lls, neg_lls, bics, ns_sorted = _compute_scores(
            diploid_sols, "diploid", gammas_noWGD
        )
        if method == "bic":
            n_dip = _select_bic(ns, bics, ns_sorted, "diploid")
        else:
            n_dip = _select_elbow(ns, neg_lls, ns_sorted, "diploid")
        ploidy2scores["diploid"] = (n_dip, ns, neg_lls, bics)
        bic_dip = float(bics[list(ns).index(n_dip)])
        logging.info(f"best diploid solution: n={n_dip}")

    n_tet = 0
    bic_tet = np.inf
    if len(tetraploid_sols) > 0:
        ns, lls, neg_lls, bics, ns_sorted = _compute_scores(
            tetraploid_sols, "tetraploid", gammas_WGD
        )
        if method == "bic":
            n_tet = _select_bic(ns, bics, ns_sorted, "tetraploid")
        else:
            n_tet = _select_elbow(ns, neg_lls, ns_sorted, "tetraploid")
        ploidy2scores["tetraploid"] = (n_tet, ns, neg_lls, bics)
        bic_tet = float(bics[list(ns).index(n_tet)])
        logging.info(f"best tetraploid solution: n={n_tet}")

    # Ploidy selection
    if len(tetraploid_sols) == 0:
        final_selection = "diploid"
    elif len(diploid_sols) == 0:
        final_selection = "tetraploid"
    elif method == "bic":
        # pick ploidy by minimum BIC
        final_selection = "diploid" if bic_dip <= bic_tet else "tetraploid"
        logging.info(
            f"ploidy selection by BIC: diploid={bic_dip:.2f}, tetraploid={bic_tet:.2f}"
        )
    else:
        # parsimony: prefer fewer clones; diploid preferred on tie
        final_selection = "diploid" if n_dip <= n_tet else "tetraploid"

    logging.info(f"final selection ({method}): {final_selection}")

    # plot
    fig, axes = plt.subplots(
        1, 2 if method == "bic" else 1, figsize=(9 if method == "bic" else 4.5, 3.5)
    )
    if method != "bic":
        axes = [axes]
    all_ns = []
    colors = {"diploid": "#1f77b4", "tetraploid": "#d62728"}
    for ploidy_key in ["diploid", "tetraploid"]:
        if ploidy_key not in ploidy2scores:
            continue
        chosen_n, ns, neg_lls, bics = ploidy2scores[ploidy_key]
        c = colors[ploidy_key]

        # neg-loglik panel (always shown)
        ax = axes[0]
        ax.plot(ns, neg_lls, "-o", color=c, markersize=5, label=ploidy_key.capitalize())
        if method == "elbow":
            idx = list(ns).index(chosen_n)
            ax.plot(chosen_n, neg_lls[idx], "*", color=c, markersize=14, zorder=5)
        all_ns.extend(ns.tolist())

        # BIC panel
        if method == "bic":
            ax2 = axes[1]
            ax2.plot(
                ns, bics, "-s", color=c, markersize=5, label=ploidy_key.capitalize()
            )
            idx = list(ns).index(chosen_n)
            ax2.plot(chosen_n, bics[idx], "*", color=c, markersize=14, zorder=5)

    ax0 = axes[0]
    ax0.set_xticks(sorted(set(int(x) for x in all_ns)))
    ax0.set_xlabel("Number of clones")
    ax0.set_ylabel("Negative log-likelihood")
    ax0.legend(framealpha=0.9)
    ax0.grid(True, alpha=0.3)

    best_n = n_dip if final_selection == "diploid" else n_tet
    if method == "bic":
        ax0.set_title("Log-likelihood")
        ax2 = axes[1]
        ax2.set_xticks(sorted(set(int(x) for x in all_ns)))
        ax2.set_xlabel("Number of clones")
        ax2.set_ylabel("BIC")
        ax2.set_title(f"BIC (best: {final_selection}, n={best_n})")
        ax2.legend(framealpha=0.9)
        ax2.grid(True, alpha=0.3)
    else:
        ax0.set_title(f"Elbow (best: {final_selection}, n={best_n})")

    fig.tight_layout()
    elbow_dir = plot_dir if plot_dir is not None else out_dir
    fig.savefig(os.path.join(elbow_dir, "elbow_curve.png"), dpi=150)
    plt.close(fig)

    return n_dip, n_tet, final_selection
