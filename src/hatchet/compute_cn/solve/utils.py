import os
import logging
import numpy as np
import pandas as pd
import kneed
import matplotlib.pyplot as plt

# A list of random states, used as a stack
random_states = []


class Random:
    """
    A context manager that pushes a random seed to the stack for reproducible results,
    and pops it on exit.
    """

    def __init__(self, seed=None):
        self.seed = seed

    def __enter__(self):
        if self.seed is not None:
            # Push current state on stack
            random_states.append(np.random.get_state())
            new_state = np.random.RandomState(self.seed)
            np.random.set_state(new_state.get_state())

    def __exit__(self, *args):
        if self.seed is not None:
            np.random.set_state(random_states.pop())


def store_solve_input(
    out_file: str, fcn_data: dict, weights: pd.Series, nbins: pd.DataFrame
):
    """Write fractional CN data to TSV."""
    fa = fcn_data["fa"]
    cluster_ids = fa.index.tolist()
    sample_ids = fa.columns.tolist()
    cols = ["fcn", "fa", "fb", "fa_lo", "fa_hi", "fb_lo", "fb_hi"]
    header = "CLUSTER\tSAMPLE\t#BINS\t" + "\t".join(cols) + "\tweight"
    with open(out_file, "w") as fd:
        fd.write(header + "\n")
        for sample in sample_ids:
            for cid in cluster_ids:
                nb = int(nbins.loc[cid, sample])
                vals = [str(fcn_data[c].loc[cid, sample]) for c in cols]
                fd.write(
                    f"{cid}\t{sample}\t{nb}\t" + "\t".join(vals) + f"\t{weights[cid]}\n"
                )


def _write_solution_tsv(
    fd, f_a, f_b, cA, cB, u, n, cluster_ids, sample_ids, header, fcn_data, nbins
):
    """Write a single solution's per-cluster/sample details to an open file."""
    cA_ = np.array(cA)
    cB_ = np.array(cB)
    u_ = np.array(u)
    exp_a = cA_ @ u_
    exp_b = cB_ @ u_

    fd.write(header + "\n")
    for ci, cid in enumerate(cluster_ids):
        for si, sample in enumerate(sample_ids):
            fa_lo = fcn_data["fa_lo"].iloc[ci, si]
            fa_hi = fcn_data["fa_hi"].iloc[ci, si]
            fb_lo = fcn_data["fb_lo"].iloc[ci, si]
            fb_hi = fcn_data["fb_hi"].iloc[ci, si]
            ea, eb = exp_a[ci, si], exp_b[ci, si]
            accepted = ea >= fa_lo and ea <= fa_hi and eb >= fb_lo and eb <= fb_hi
            fields = [
                cid,
                sample,
                int(nbins.loc[cid, sample]),
                f_a.loc[cid, sample],
                f_b.loc[cid, sample],
                ea,
                eb,
                fa_lo,
                fa_hi,
                fb_lo,
                fb_hi,
            ]
            for oi in range(n):
                fields.extend([f"{cA[ci][oi]}|{cB[ci][oi]}", u[oi][si]])
            fields.append(accepted)
            fd.write("\t".join(str(v) for v in fields) + "\n")


def dedup_solutions(solutions, u_atol=1e-3):
    """Remove duplicate (obj, cA, cB, u) tuples that differ only by clone ordering.

    Two solutions are duplicates if their CN profiles (cA, cB) are identical
    after sorting clones into a canonical order and their clone proportions (u)
    agree within ``u_atol``.

    Args:
        solutions: List of ``(obj, cA, cB, u)`` tuples.
        u_atol: Absolute tolerance for comparing clone proportions.

    Returns:
        Deduplicated list (order preserved, first occurrence kept).
    """
    if len(solutions) <= 1:
        return solutions

    def _canon(cA, cB, u):
        arr = np.array(cA + cB)
        u_arr = np.array(u)
        order = np.lexsort(arr[::-1])
        return arr[:, order], u_arr[order]

    keep = []
    canon_cache = []
    for sol in solutions:
        _, cA, cB, u = sol
        cn_c, u_c = _canon(cA, cB, u)
        is_dup = False
        for cn_k, u_k in canon_cache:
            if np.array_equal(cn_c, cn_k) and np.allclose(u_c, u_k, atol=u_atol):
                is_dup = True
                break
        if not is_dup:
            keep.append(sol)
            canon_cache.append((cn_c, u_c))
    n_removed = len(solutions) - len(keep)
    if n_removed > 0:
        logging.info(
            f"removed {n_removed} duplicate solutions (same CN up to clone reordering)"
        )
    return keep


def store_instance_tofile(
    pool_instances: dict,
    f_a: pd.DataFrame,
    f_b: pd.DataFrame,
    outdir: str,
    solve_mode: str,
    n: int,
    fcn_data: dict,
    nbins: pd.DataFrame,
):
    """Store all solution detail TSVs with naming: sol{pparam}_pool{idx}.tsv."""
    cluster_ids = f_a.index.tolist()
    sample_ids = f_a.columns.tolist()
    clone_cols = ["cn_normal\tu_normal"] + [
        f"cn_clone{i}\tu_clone{i}" for i in range(1, n)
    ]
    cols = (
        [
            "CLUSTER",
            "SAMPLE",
            "#BINS",
            "f_a",
            "f_b",
            "exp_f_a",
            "exp_f_b",
            "fa_lo",
            "fa_hi",
            "fb_lo",
            "fb_hi",
        ]
        + clone_cols
        + ["ci_accepted"]
    )
    header = "\t".join(cols)

    for pparam, solutions in pool_instances.items():
        for pool_idx, (obj, cA, cB, u) in enumerate(solutions):
            path = os.path.join(outdir, f"{solve_mode}_sol{pparam}_pool{pool_idx}.tsv")
            with open(path, "w") as fd:
                _write_solution_tsv(
                    fd,
                    f_a,
                    f_b,
                    cA,
                    cB,
                    u,
                    n,
                    cluster_ids,
                    sample_ids,
                    header,
                    fcn_data,
                    nbins,
                )


def compute_pairwise_cnt(cA_list, cB_list, bbcs, cluster_ids):
    """Compute per-pair CNT and weighted-CNT distances at cluster level.

    Orders clusters by genomic position using bin-level BBC data, then
    computes pairwise distances among **tumor clones** (excluding normal).

    Parameters
    ----------
    cA_list, cB_list : list of list
        Cluster-level allele-specific CN, shape (n_clusters, n_clones).
    bbcs : pd.DataFrame
        Bin-level BBC with ``#CHR``, ``START``, ``END``, ``SAMPLE``, ``CLUSTER``.
    cluster_ids : list
        Ordered cluster identifiers matching rows of cA_list/cB_list.

    Returns
    -------
    dict
        ``{"CNT_ci_cj": float, "WCNT_ci_cj": float, ...}`` for each
        tumor clone pair (i < j), using 1-based clone indices.
        Returns empty dict when fewer than 2 tumor clones.
    """
    from hatchet.compute_cn.solve.cnt_distance import (
        compute_cnt_distances,
        compute_weighted_cnt_distances,
    )

    cA_arr = np.array(cA_list)
    cB_arr = np.array(cB_list)
    n_clones = cA_arr.shape[1]
    n_tumor = n_clones - 1  # exclude normal (col 0)
    if n_tumor < 2:
        return {}

    # Order clusters by genomic position; use sorted-first sample for determinism
    first_sample = sorted(bbcs["SAMPLE"].unique())[0]
    bbc_s = bbcs[bbcs["SAMPLE"] == first_sample]
    spans = (
        bbc_s.groupby("CLUSTER")
        .agg(
            chr_first=("#CHR", "first"),
            start=("START", "min"),
            end=("END", "max"),
        )
        .reset_index()
    )
    spans = spans[spans["CLUSTER"].isin(cluster_ids)]

    # Map cluster_ids to cA/cB row indices
    cid_to_row = {cid: i for i, cid in enumerate(cluster_ids)}
    spans["_row"] = spans["CLUSTER"].map(cid_to_row)
    spans = spans.sort_values(["chr_first", "start"]).reset_index(drop=True)

    order = spans["_row"].to_numpy()
    cA_ordered = cA_arr[order]
    cB_ordered = cB_arr[order]
    seg_lengths = (spans["end"] - spans["start"]).to_numpy(dtype=float)

    m = len(order)
    chr_b = np.zeros(m, dtype=bool)
    chr_b[0] = True
    chr_b[1:] = (
        spans["chr_first"].iloc[1:].values != spans["chr_first"].iloc[:-1].values
    )

    # Exclude normal clone (col 0)
    cA_t = cA_ordered[:, 1:]
    cB_t = cB_ordered[:, 1:]

    cnt_dist = compute_cnt_distances(cA_t, cB_t, chr_b)
    wcnt_dist = compute_weighted_cnt_distances(
        cA_t, cB_t, chr_b, seg_lengths=seg_lengths
    )

    result = {}
    for i in range(n_tumor):
        for j in range(i + 1, n_tumor):
            ci, cj = i + 1, j + 1  # 1-based clone indices
            cnt_val = cnt_dist[i, j]
            wcnt_val = wcnt_dist[i, j]
            result[f"CNT_c{ci}_c{cj}"] = cnt_val if np.isfinite(cnt_val) else "inf"
            result[f"WCNT_c{ci}_c{cj}"] = wcnt_val if np.isfinite(wcnt_val) else "inf"

    # Sum of CNT distances from clone 1 (MRCA) to all other tumor clones
    cnt_from_c1 = cnt_dist[0, 1:]  # row 0 = clone1, cols 1+ = clone2,3,...
    wcnt_from_c1 = wcnt_dist[0, 1:]
    result["CNT_from_c1"] = (
        float(cnt_from_c1.sum()) if np.all(np.isfinite(cnt_from_c1)) else "inf"
    )
    result["WCNT_from_c1"] = (
        float(wcnt_from_c1.sum()) if np.all(np.isfinite(wcnt_from_c1)) else "inf"
    )
    return result


def compute_individual_objs(
    pname: str,
    weights: pd.Series,
    fA: pd.DataFrame,
    fB: pd.DataFrame,
    cA: list,
    cB: list,
    u: list,
):
    """Compute the IMF and regularisation objective values for one solution.

    Args:
        pname: Regularisation objective name (e.g. ``"DROOT_SUM"``).
            If not in the known set, the regularisation objective is 0.
        weights: Per-cluster weights, shape (n_clusters,).
        fA: Observed fractional A copy numbers (clusters x samples).
        fB: Observed fractional B copy numbers (clusters x samples).
        cA: Integer allele-A copy numbers, shape (n_clusters, n_clones).
        cB: Integer allele-B copy numbers, shape (n_clusters, n_clones).
        u: Clone proportions, shape (n_clones, n_samples).

    Returns:
        List ``[imf_obj, reg_obj]`` of float values.
    """
    w_ = weights.to_numpy().reshape((len(weights), 1))
    fA_ = fA.to_numpy()
    fB_ = fB.to_numpy()
    cA_ = np.array(cA)
    cB_ = np.array(cB)
    u_ = np.array(u)

    imf_obj = compute_obj_IMF(w_, fA_, fB_, cA_, cB_, u_)
    reg_objs = {
        "MAXCN": compute_obj_MAXCN,
        "DROOT_SUM": compute_obj_DROOT_SUM,
        "DADJ_SUM": compute_obj_DADJ_SUM,
        "DMRCA_SUM": compute_obj_DMRCA_SUM,
    }
    sub_obj = reg_objs[pname](w_, fA_, fB_, cA_, cB_, u_) if pname in reg_objs else 0.0
    return [imf_obj, sub_obj]


def compute_obj_IMF(weights, fA, fB, cA, cB, u):
    """Compute the weighted IMF (integer matrix factorisation) objective."""
    leftA_w = weights * np.abs(fA - cA @ u)
    leftB_w = weights * np.abs(fB - cB @ u)
    obj = np.sum(leftA_w) + np.sum(leftB_w)
    return obj


def compute_obj_DROOT_SUM(weights, _fA, _fB, cA, cB, _u):
    """Compute weighted sum of CN distance from the diploid root (1,1) for tumor clones.

    All regularisation objective functions share the same call signature
    ``(weights, fA, fB, cA, cB, u)`` so they can be dispatched uniformly via a
    dict in ``compute_individual_objs``.  Parameters not needed by this objective
    are prefixed with ``_``.
    """
    distA = weights * np.abs(cA[:, 1:] - cA[:, :1])
    distB = weights * np.abs(cB[:, 1:] - cB[:, :1])
    obj = np.sum(distA) + np.sum(distB)
    return obj


def compute_obj_DMRCA_SUM(weights, _fA, _fB, cA, cB, _u):
    """Compute weighted sum of Manhattan distance from the MRCA clone to all subclonal clones.

    Clone layout: col 0 = normal, col 1 = MRCA, cols 2+ = subclonal tumor clones.
    Returns 0.0 when n < 3 (no subclonal clones exist beyond the MRCA).

    All regularisation objective functions share the same call signature
    ``(weights, fA, fB, cA, cB, u)`` so they can be dispatched uniformly via a
    dict in ``compute_individual_objs``.  Parameters not needed by this objective
    are prefixed with ``_``.
    """
    if cA.shape[1] < 3:
        return 0.0
    distA = weights * np.abs(cA[:, 2:] - cA[:, 1:2])
    distB = weights * np.abs(cB[:, 2:] - cB[:, 1:2])
    return np.sum(distA) + np.sum(distB)


def compute_obj_DADJ_SUM(weights, _fA, _fB, cA, cB, _u):
    """Compute weighted pairwise Hamming distance between clone CN states.

    All regularisation objective functions share the same call signature
    ``(weights, fA, fB, cA, cB, u)`` so they can be dispatched uniformly via a
    dict in ``compute_individual_objs``.  Parameters not needed by this objective
    are prefixed with ``_``.
    """
    obj = 0
    num_clusters, num_clones = cA.shape
    for cluster_idx in range(num_clusters):
        obj_cluster = 0.0
        for clone_i in range(num_clones - 1):
            for clone_j in range(clone_i + 1, num_clones):
                obj_cluster += abs(cA[cluster_idx, clone_i] - cA[cluster_idx, clone_j])
                obj_cluster += abs(cB[cluster_idx, clone_i] - cB[cluster_idx, clone_j])
        obj += weights[cluster_idx, 0] * obj_cluster
    return obj


def compute_obj_MAXCN(weights, _fA, _fB, cA, cB, _u):
    """Compute weighted sum of maximum CN per cluster across tumor clones.

    All regularisation objective functions share the same call signature
    ``(weights, fA, fB, cA, cB, u)`` so they can be dispatched uniformly via a
    dict in ``compute_individual_objs``.  Parameters not needed by this objective
    are prefixed with ``_``.
    """
    maxA_w = np.dot(np.max(cA[:, 1:], axis=1), weights)[0]
    maxB_w = np.dot(np.max(cB[:, 1:], axis=1), weights)[0]
    return maxA_w + maxB_w


def filter_non_pareto(points: np.ndarray):
    """Return a boolean mask of Pareto-optimal points.

    A point is Pareto-optimal if no other point dominates it.
    Point j dominates point i if j is no worse on all objectives and
    strictly better on at least one.

    Args:
        points: Array of shape (n_points, n_objectives) with minimisation objectives.

    Returns:
        Boolean array of length n_points; True where the point is Pareto-optimal.
    """
    is_pareto = np.ones(len(points), dtype=bool)
    for i in range(len(points)):
        for j in range(len(points)):
            if i == j:
                continue
            if (
                points[j, 0] <= points[i, 0]
                and points[j, 1] <= points[i, 1]
                and (points[j, 0] < points[i, 0] or points[j, 1] < points[i, 1])
            ):
                is_pareto[i] = False
                break
    return is_pareto


def model_select_elbow(df: pd.DataFrame, xid: str, yid: str, pareto_img: str):
    """Select the best solution from a two-objective Pareto front using the elbow criterion.

    Given a set of solutions with two minimising objectives ``xid`` and ``yid``:
    1. Select the Pareto-optimal subset.
    2. Pick the best solution using the elbow/knee criterion on the Pareto curve.
       If elbow detection fails, fall back to the Pareto point with the lowest
       ``yid`` (best IMF fit).

    Pareto points are sorted by ``xid`` ascending before elbow detection so that
    KneeLocator always receives a monotone-increasing x sequence.
    """
    df.loc[:, "is_pareto"] = filter_non_pareto(df[[xid, yid]].to_numpy())
    df.loc[:, "selected"] = ""

    pids = df.loc[df["is_pareto"]].index.to_numpy()

    logging.info(f"model selection, #pareto={len(pids)}/{len(df)}")

    if len(pids) == 0:
        logging.info("WARN! no Pareto points found; falling back to row 0")
        sol_index = df.index[0]
        df.loc[sol_index, "selected"] = "*"
        return df, sol_index

    pareto_df = df.loc[pids].sort_values(xid)
    pids_sorted = pareto_df.index.to_numpy()
    xs = pareto_df[xid].to_numpy()
    ys = pareto_df[yid].to_numpy()

    sol_index = pids_sorted[-1]

    # Run elbow detection only when enough Pareto points exist for a meaningful curve.
    elbow_x, elbow_y = None, None
    if len(pids_sorted) >= 3:
        kl = kneed.KneeLocator(x=xs, y=ys, curve="convex", direction="decreasing")
        elbow_x, elbow_y = kl.elbow, kl.elbow_y
        if elbow_x is not None and elbow_x != xs[0]:
            sol_indices = np.where(ys <= elbow_y)[0]
            if len(sol_indices) > 0:
                sol_index = pids_sorted[sol_indices[0]]
                logging.info(f"Model selection elbow at index={sol_index}")

    # Always write the Pareto plot when a path is given, regardless of point count.
    if pareto_img is not None:
        fig, ax = plt.subplots()
        if (~df["is_pareto"]).any():
            ax.scatter(
                df.loc[~df["is_pareto"], xid].to_numpy(),
                df.loc[~df["is_pareto"], yid].to_numpy(),
                c="gray",
                marker="x",
                alpha=0.6,
                label="non-Pareto",
            )
        ax.plot(xs, ys, c="green", linewidth=1, zorder=2)
        ax.scatter(xs, ys, c="green", marker="o", zorder=3, label="Pareto")
        if elbow_x is not None:
            ax.axvline(elbow_x, linestyle="--", color="steelblue", label="knee/elbow")
        ax.set_xlabel(xid)
        ax.set_ylabel(yid)
        ax.set_title("Model Selection Pareto Curve")
        ax.legend()
        plt.savefig(pareto_img, dpi=150)
        plt.close()

    df.loc[sol_index, "selected"] = "*"
    return df, sol_index


def _count_ci_violations(cA, cB, u, fcn_data, nbins):
    """Count bins where expected FCN falls outside the CI.

    nbins: DataFrame (clusters x samples) of bin counts per cluster/sample.
    """
    cA_ = np.array(cA)
    cB_ = np.array(cB)
    u_ = np.array(u)
    exp_a = cA_ @ u_  # (n_clusters, n_samples)
    exp_b = cB_ @ u_

    fa_lo = fcn_data["fa_lo"].to_numpy()
    fa_hi = fcn_data["fa_hi"].to_numpy()
    fb_lo = fcn_data["fb_lo"].to_numpy()
    fb_hi = fcn_data["fb_hi"].to_numpy()
    nb = nbins.to_numpy()

    violations = (exp_a < fa_lo) | (exp_a > fa_hi) | (exp_b < fb_lo) | (exp_b > fb_hi)
    n_violated_bins = int(np.sum(violations * nb))
    n_total_bins = int(np.sum(nb))
    ratio = n_violated_bins / n_total_bins if n_total_bins > 0 else float("nan")
    return n_violated_bins, ratio


def model_selection_instance(
    f_a: pd.DataFrame,
    f_b: pd.DataFrame,
    weights: pd.Series,
    pool_instances: dict,
    pname: str,
    solve_mode: str,
    outdir: str,
    fcn_data: dict = None,
    nbins: pd.DataFrame = None,
    pareto_img: str = None,
):
    """Select the best solution from a regularisation path using the elbow criterion.

    Args:
        f_a: DataFrame of observed fractional A copy numbers (clusters x samples).
        f_b: DataFrame of observed fractional B copy numbers (clusters x samples).
        weights: Series of per-cluster weights.
        pool_instances: Mapping ``{pparam: [(obj, cA, cB, u), ...]}``, where index 0
            is the primary solution and index 1+ are pool alternatives.
        pname: Regularisation objective name (e.g. ``"DROOT_SUM"``), or None for
            raw (unregularised) selection.
        solve_mode: Label used in output filenames (e.g. ``"cd"`` or ``"ilp"``).
        outdir: Directory for TSV and PNG output, or None to suppress all file I/O.
            Ignored when ``pareto_img`` is supplied explicitly.
        fcn_data: Dict of fractional-CN DataFrames used for CI-violation counting,
            or None to skip CI accounting.
        nbins: DataFrame of bin counts per cluster/sample used for CI-violation
            counting, or None to skip CI accounting.
        pareto_img: Explicit path for the Pareto-curve PNG.  Takes priority over
            the ``outdir``-based auto-constructed path.  Pass None (default) to
            let the function derive the path from ``outdir``.

    Returns:
        Tuple ``(best_solution, imf_obj, key)`` where ``best_solution`` is the
        ``(obj, cA, cB, u)`` tuple for the selected solution, ``imf_obj`` is its
        unregularised IMF objective value, and ``key`` is ``(pparam, pool_idx)``.
    """
    assert len(pool_instances) > 0, "no solutions to select from"

    def _build_row(pparam, pool_idx, tobj, cA, cB, u):
        [imf_obj, reg_obj] = compute_individual_objs(
            pname, weights, f_a, f_b, cA, cB, u
        )
        lambda_val = float(pparam)
        if fcn_data is not None and nbins is not None:
            n_viol, viol_ratio = _count_ci_violations(cA, cB, u, fcn_data, nbins)
        else:
            n_viol, viol_ratio = 0, 0.0
        return [
            pparam,
            pool_idx,
            lambda_val,
            tobj,
            imf_obj,
            reg_obj,
            n_viol,
            round(viol_ratio, 4),
        ]

    data = []
    all_solutions = {}
    for pparam, solutions in sorted(pool_instances.items(), key=lambda tp: tp[0]):
        for pool_idx, (tobj, cA, cB, u) in enumerate(solutions):
            data.append(_build_row(pparam, pool_idx, tobj, cA, cB, u))
            all_solutions[(pparam, pool_idx)] = (tobj, cA, cB, u)

    if pname is None or len(data) == 1:
        key = next(iter(all_solutions))
        return all_solutions[key], all_solutions[key][0], key

    columns = [
        "instance_id",
        "pool_idx",
        "Lambda",
        "Objective",
        "IMF-objective",
        f"{pname}-objective",
        "ci_violations",
        "ci_violation_ratio",
    ]

    df = pd.DataFrame(data=data, columns=columns)

    # Deduplicate by actual CN states (up to clone reordering) + purity tolerance
    keys = list(all_solutions.keys())
    all_sols_list = [all_solutions[key] for key in keys]
    deduped = dedup_solutions(all_sols_list)
    if len(deduped) < len(all_sols_list):
        deduped_set = {id(s) for s in deduped}
        keep_mask = [id(all_sols_list[i]) in deduped_set for i in range(len(keys))]
        df = df.loc[keep_mask].reset_index(drop=True)
        all_solutions = {
            keys[i]: all_sols_list[i] for i, should_keep in enumerate(keep_mask) if should_keep
        }

    if pareto_img is None and outdir is not None:
        pareto_img = os.path.join(outdir, f"pareto_curve.{solve_mode}.{pname}.png")

    df, sol_index = model_select_elbow(
        df,
        f"{pname}-objective",
        "IMF-objective",
        pareto_img,
    )

    df = df.sort_values(["instance_id", "pool_idx"]).reset_index(drop=True)

    if outdir is not None:
        sols_parent = os.path.dirname(outdir)
        subdir_name = os.path.basename(outdir)
        df.to_csv(
            os.path.join(
                sols_parent, f"{subdir_name}.solutions.{solve_mode}.{pname}.tsv"
            ),
            sep="\t",
            header=True,
            index=False,
        )

    sel = df.loc[sol_index]
    key = (sel["instance_id"], int(sel["pool_idx"]))
    return all_solutions[key], sel["IMF-objective"], key
