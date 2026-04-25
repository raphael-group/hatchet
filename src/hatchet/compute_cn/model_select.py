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


def model_selection_ploidy(
    chosen_sols: dict,
    out_dir: str,
    scaling: dict,
    segs: pd.DataFrame,
    method="elbow",
):
    """Select best n and ploidy across ploidies.

    Args:
        chosen_sols: {ploidy: {n: best_sol_dict}} from the solve loop.
        out_dir: output directory (for reading UCN files).
        scaling: dict from get_scaling_factor with per-ploidy gammas.
        segs: cluster-level SEG DataFrame.
        method: "elbow" or "bic".

    Returns:
        (best_ploidy, best_n, fig) where fig is the elbow/BIC plot.
    """

    def _compute_scores(ploidy):
        gammas = scaling[ploidy]["gammas"]
        ns_sorted = sorted(chosen_sols[ploidy].keys())
        first_ucn = os.path.join(
            out_dir, f"results.{ploidy}.n{ns_sorted[0]}.bbc.ucn.tsv"
        )
        ll_n1, nobs_n1, n_clusters, n_samples = _compute_loglik_from_ucn(
            first_ucn, 1, gammas, segs
        )
        logging.info(f"{ploidy}: n=1 (baseline), loglik={ll_n1:.4f}")

        ns_all = [1] + ns_sorted
        lls = [ll_n1]
        for clone_n in ns_sorted:
            ucn_file = os.path.join(out_dir, f"results.{ploidy}.n{clone_n}.bbc.ucn.tsv")
            ll, nobs, n_clusters, n_samples = _compute_loglik_from_ucn(
                ucn_file, clone_n, gammas, segs
            )
            lls.append(ll)
            logging.info(f"{ploidy}: n={clone_n}, loglik={ll:.4f}")

        ns = np.array(ns_all, dtype=np.int32)
        lls = np.array(lls)
        bics = np.array(
            [
                -2 * lls[i]
                + _count_params(ns[i], n_clusters, n_samples) * np.log(nobs_n1)
                for i in range(len(ns))
            ]
        )
        return ns, lls, -lls, bics, ns_sorted

    def _select(ns, neg_lls, bics, ns_sorted, ploidy):
        if method == "bic":
            idx = int(np.argmin(bics))
            chosen = max(int(ns[idx]), int(ns_sorted[0]))
            logging.info(f"{ploidy}: chosen n={chosen} via BIC")
        else:
            chosen = int(ns_sorted[0])
            if len(ns) >= 3:
                kl = kneed.KneeLocator(
                    x=ns, y=neg_lls, curve="convex", direction="decreasing"
                )
                if kl.elbow is not None:
                    chosen = max(int(ns_sorted[0]), int(kl.elbow))
            logging.info(f"{ploidy}: chosen n={chosen} via loglik elbow")
        return chosen

    results = {}  # ploidy -> (chosen_n, ns, neg_lls, bics)
    for ploidy in chosen_sols:
        ns, lls, neg_lls, bics, ns_sorted = _compute_scores(ploidy)
        chosen_n = _select(ns, neg_lls, bics, ns_sorted, ploidy)
        results[ploidy] = (chosen_n, ns, neg_lls, bics)
        logging.info(f"best {ploidy} solution: n={chosen_n}")

    # Ploidy selection
    ploidies = list(results.keys())
    if len(ploidies) == 1:
        best_ploidy = ploidies[0]
    elif method == "bic":
        bic_vals = {
            p: float(bics[list(ns).index(n)]) for p, (n, ns, _, bics) in results.items()
        }
        best_ploidy = min(bic_vals, key=bic_vals.get)
        logging.info(f"ploidy selection by BIC: {bic_vals}")
    else:
        n_vals = {p: n for p, (n, _, _, _) in results.items()}
        best_ploidy = (
            "diploid"
            if n_vals.get("diploid", 999) <= n_vals.get("tetraploid", 999)
            else "tetraploid"
        )

    best_n = results[best_ploidy][0]
    logging.info(f"final selection ({method}): {best_ploidy}, n={best_n}")

    # Plot
    fig, axes = plt.subplots(
        1, 2 if method == "bic" else 1, figsize=(9 if method == "bic" else 4.5, 3.5)
    )
    if method != "bic":
        axes = [axes]
    colors = {"diploid": "#1f77b4", "tetraploid": "#d62728"}
    all_ns = []
    for ploidy, (chosen_n, ns, neg_lls, bics) in results.items():
        c = colors.get(ploidy, "gray")
        axes[0].plot(
            ns, neg_lls, "-o", color=c, markersize=5, label=ploidy.capitalize()
        )
        if method == "elbow":
            idx = list(ns).index(chosen_n)
            axes[0].plot(chosen_n, neg_lls[idx], "*", color=c, markersize=14, zorder=5)
        all_ns.extend(ns.tolist())
        if method == "bic":
            axes[1].plot(
                ns, bics, "-s", color=c, markersize=5, label=ploidy.capitalize()
            )
            idx = list(ns).index(chosen_n)
            axes[1].plot(chosen_n, bics[idx], "*", color=c, markersize=14, zorder=5)

    axes[0].set_xticks(sorted(set(int(x) for x in all_ns)))
    axes[0].set_xlabel("Number of clones")
    axes[0].set_ylabel("Negative log-likelihood")
    axes[0].legend(framealpha=0.9)
    axes[0].grid(True, alpha=0.3)
    if method == "bic":
        axes[0].set_title("Log-likelihood")
        axes[1].set_xticks(sorted(set(int(x) for x in all_ns)))
        axes[1].set_xlabel("Number of clones")
        axes[1].set_ylabel("BIC")
        axes[1].set_title(f"BIC (best: {best_ploidy}, n={best_n})")
        axes[1].legend(framealpha=0.9)
        axes[1].grid(True, alpha=0.3)
    else:
        axes[0].set_title(f"Elbow (best: {best_ploidy}, n={best_n})")
    fig.tight_layout()

    chosen_n = {p: n for p, (n, _, _, _) in results.items()}
    return best_ploidy, best_n, chosen_n, fig


def model_select_elbow_from_regularization(pool_instances):
    """Select the best solution from a regularization-path pool using elbow criterion.

    Each solution in pool_instances must have imf_obj and reg_obj.

    Args:
        pool_instances: {sol_id: {"imf_obj", "reg_obj", "cA", "cB", "u", ...}}.

    Returns (best_sol_id, df) where df has Pareto/selected annotations.
    """
    assert len(pool_instances) > 0, "no solutions to select from"

    reg_col = "REG"
    imf_col = "IMF"

    sol_ids = sorted(pool_instances)
    df = pd.DataFrame(
        {
            "instance_id": sol_ids,
            imf_col: [pool_instances[s]["imf_obj"] for s in sol_ids],
            reg_col: [pool_instances[s]["reg_obj"] for s in sol_ids],
        }
    )

    df["is_pareto"] = True
    df["selected"] = ""
    for i in range(len(df)):
        for j in range(len(df)):
            if (
                i != j
                and df[imf_col].iloc[j] <= df[imf_col].iloc[i]
                and df[reg_col].iloc[j] <= df[reg_col].iloc[i]
                and (
                    df[imf_col].iloc[j] < df[imf_col].iloc[i]
                    or df[reg_col].iloc[j] < df[reg_col].iloc[i]
                )
            ):
                df.iloc[i, df.columns.get_loc("is_pareto")] = False
                break
    pids = df.loc[df["is_pareto"]].index.to_numpy()
    if len(pids) <= 1:
        df["is_pareto"] = True
        df["selected"] = "*"
        return sol_ids[0], df
    logging.info(f"model selection, #pareto={len(pids)}/{len(df)}")

    pareto_df = df.loc[pids].sort_values(reg_col)
    pids_sorted = pareto_df.index.to_numpy()
    xs = pareto_df[reg_col].to_numpy()
    ys = pareto_df[imf_col].to_numpy()

    best_idx = pids_sorted[-1]
    elbow_x, elbow_y = None, None
    if len(pids_sorted) >= 3:
        kl = kneed.KneeLocator(x=xs, y=ys, curve="convex", direction="decreasing")
        elbow_x, elbow_y = kl.elbow, kl.elbow_y
        if elbow_x is not None and elbow_x != xs[0]:
            sol_indices = np.where(ys <= elbow_y)[0]
            if len(sol_indices) > 0:
                best_idx = pids_sorted[sol_indices[0]]
                logging.info(f"Model selection elbow at index={best_idx}")
    df.loc[best_idx, "selected"] = "*"
    best_sol_id = df.loc[best_idx, "instance_id"]
    return best_sol_id, df
