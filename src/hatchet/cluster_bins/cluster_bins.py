import os
import time
import logging
import shutil

import numpy as np
import pandas as pd
from hatchet.utils import *
from hatchet.cluster_bins.cluster_utils import *
from hatchet.cluster_bins.hmm.hmm_init import *
from hatchet.cluster_bins.hmm.hmm_transitions import *
from hatchet.cluster_bins.hmm.hmm_model import *
from hatchet.plot.plot_1d2d import plot_rdr_baf


def run(args=None):
    """Main entry point for the cluster-bins step.
    1. run BAF+RDR factorial HMM with K cluster states for K in [minK, maxK].
    2. select best K by BIC or ICL.

    Input data (read from bb_dir/):
        bb.tsv.gz         — bin metadata (chr, start, end, region_id, switchprobs)
        bb.rdr.npz        — (N, M+1) RDR matrix
        bb.depth.npz      - (N, M+1) Read-depth matrix
        bb.{A,B,T}allele.npz — (N, M+1) allele count matrices

    Outputs (written to bbc_dir/):
        bulk.bbc / bulk.seg     — BBC and SEG files for the optimal K
        labels/bulkK.bbc|seg    — per-K results
        cluster_infos/          — per-K cluster-label TSVs
        plots/                  — ELBO traces, RDR-BAF scatter, model score
        model_scores.tsv / .png — BIC or ICL scores across K

    Args:
        args: dict or Namespace of CLI arguments (see hatchet_parser.py).
    """
    args = normalize_args(args)
    setup_logging(args)
    log_arguments(args)
    logging.info("cluster bins")
    _log_done = log_step_start()

    bb_dir = args["bb_dir"]
    bb_file = os.path.join(bb_dir, "bb.tsv.gz")
    sample_file = os.path.join(bb_dir, "sample_ids.tsv")
    rdr_mfile = os.path.join(bb_dir, "bb.rdr.npz")
    depth_mfile = os.path.join(bb_dir, "bb.depth.npz")
    a_mfile = os.path.join(bb_dir, "bb.Aallele.npz")
    b_mfile = os.path.join(bb_dir, "bb.Ballele.npz")
    t_mfile = os.path.join(bb_dir, "bb.Tallele.npz")
    genome_size = args["genome_size"]
    out_dir = args["bbc_dir"]

    out_bbc = os.path.join(out_dir, "bulk.bbc")
    out_seg = os.path.join(out_dir, "bulk.seg")
    if not args["force"] and os.path.exists(out_bbc) and os.path.exists(out_seg):
        logging.info(
            f"skip cluster-bins: {out_bbc} and {out_seg} already exist (use --force to re-run)"
        )
        _log_done("cluster-bins")
        return

    min_tau = args["min_tau"]
    max_tau = args["max_tau"]
    baf_eps = args["baf_eps"]
    min_covar = args["min_covar"]
    ig_alpha = args["ig_alpha"]
    tau_iters = args["tau_iters"]
    diag_t = args["t"]

    minK = args["minK"]
    maxK = args["maxK"]
    restarts = args["restarts"]
    top_restarts = (
        args["top_restarts"] if args["top_restarts"] is not None else restarts
    )
    n_local_trials = args["n_local_trials"]
    n_iter = args["niters"]

    seed = args["seed"]
    decode_method = args["decode_method"]
    score_method = args["score_method"]
    log_rdr = args["log_rdr"]
    init_method = args["init_method"]
    training_method = args["training_method"]
    baf_k_start = 0 if args["free_baf_c0"] else 1

    os.makedirs(out_dir, exist_ok=True)
    add_file_logging(out_dir, "cluster-bins")
    label_dir = os.path.join(out_dir, "labels")
    plot_dir = os.path.join(out_dir, "plots")
    os.makedirs(label_dir, exist_ok=True)
    os.makedirs(plot_dir, exist_ok=True)

    ##################################################
    logging.info("load arguments")
    _, samples, no_normal = read_sample_file(sample_file)
    tumor_sidx = 0 if no_normal else 1
    tumor_samples = samples[tumor_sidx:]

    bbs = pd.read_table(bb_file, sep="\t")

    X_depths = np.load(depth_mfile)["mat"].astype(np.float32)
    X_depths_tumor = X_depths[:, tumor_sidx:]

    X_rdrs = np.load(rdr_mfile)["mat"].astype(np.float32)
    X_alphas_all = np.load(a_mfile)["mat"].astype(np.int32)
    X_betas_all = np.load(b_mfile)["mat"].astype(np.int32)
    X_totals_all = np.load(t_mfile)["mat"].astype(np.int32)
    if not no_normal:
        X_alphas_normal = X_alphas_all[:, 0]
        X_betas_normal = X_betas_all[:, 0]
    X_alphas = X_alphas_all[:, tumor_sidx:]
    X_betas = X_betas_all[:, tumor_sidx:]
    X_totals = X_totals_all[:, tumor_sidx:]
    nbbs, ntumor_samples = X_rdrs.shape
    assert len(bbs) == nbbs, f"unmatched {len(bbs)} and {nbbs}"

    X_log_rdrs = np.log(np.clip(X_rdrs, 1e-6, None)).astype(np.float32)
    X_hmm_rdrs = X_log_rdrs if log_rdr else X_rdrs

    ##################################################
    logging.info("prepare HMM inputs")
    X_lengths = (
        bbs.groupby(by="region_id", sort=False).agg("size").to_numpy(dtype=np.int64)
    )
    nsegments = len(X_lengths)
    assert np.sum(X_lengths) == nbbs, f"unmatched {np.sum(X_lengths)} and {nbbs}"
    logging.info(f"#bbs={nbbs}, #segments={nsegments}")

    X_bafs = np.clip(X_betas / X_totals, baf_eps, 1 - baf_eps)
    switchprobs = bbs["switchprobs"].to_numpy()

    X_hmm_rdrs = np.ascontiguousarray(X_hmm_rdrs, dtype=np.float64)
    X_alphas = np.ascontiguousarray(X_alphas, dtype=np.float64)  # (N, M)
    X_betas = np.ascontiguousarray(X_betas, dtype=np.float64)  # (N, M)
    X_totals = np.ascontiguousarray(X_totals, dtype=np.float64)  # (N, M)
    X_bafs = np.ascontiguousarray(X_bafs, dtype=np.float64)  # (N, M)
    log_switchprobs = np.ascontiguousarray(np.log(switchprobs), dtype=np.float64)
    log_stayprobs = np.ascontiguousarray(np.log(1 - switchprobs), dtype=np.float64)

    ##################################################
    logging.debug(
        f"nbbs={nbbs}, ntumor_samples={ntumor_samples}, bbs.columns={list(bbs.columns)}"
    )
    bbcs = pd.DataFrame(
        {
            "#CHR": np.repeat(bbs["#CHR"].to_numpy(), ntumor_samples),
            "START": np.repeat(bbs["START"].to_numpy(), ntumor_samples),
            "END": np.repeat(bbs["END"].to_numpy(), ntumor_samples),
            "SAMPLE": np.tile(tumor_samples, nbbs),
            "#SNPS": np.repeat(bbs["#SNPS"].to_numpy(), ntumor_samples),
        }
    )
    bbcs["CLUSTER"] = 0
    bbcs["RD"] = X_rdrs.ravel()
    bbcs["COV"] = X_depths_tumor.ravel()

    bbs["CLUSTER"] = 0
    bbs["PHASE"] = 0.0
    bbs["PHASE_POSTS"] = 0.0

    ##################################################
    plot_rdr_baf(
        tumor_samples,
        bbs,
        X_bafs,
        X_rdrs,
        genome_size,
        xlab="BAF",
        ylab="RDR",
        out_dir=plot_dir,
        out_prefix="raw_",
        dpi=100,
    )

    if not no_normal:
        baf_taus0 = estimate_BB_dispersion_normal(
            X_alphas_normal,
            X_betas_normal,
            ntumor_samples,
            min_tau=min_tau,
            max_tau=max_tau,
        )
    else:
        baf_taus0 = estimate_BB_dispersion_segment(
            X_alphas,
            X_betas,
            X_bafs,
            X_lengths,
            ntumor_samples,
            min_tau=min_tau,
            max_tau=max_tau,
        )
    logging.info("estimated BAF per-sample dispersion:      %s", np.round(baf_taus0, 3))
    rdr_vars0 = estimate_rdr_vars(X_hmm_rdrs, X_lengths, min_var=min_covar)
    logging.info("estimated RDR per-sample variance:      %s", np.round(rdr_vars0, 3))

    ig_beta = rdr_vars0 * (ig_alpha + 1)
    logging.info("IG prior: alpha=%.2f, beta=%s", ig_alpha, np.round(ig_beta, 6))

    chrom_sizes = read_genome_sizes(genome_size)
    DEBUG = logging.getLogger().isEnabledFor(logging.DEBUG)

    logging.info("HMM init method: %s", init_method)
    if init_method == "kmeans_plus_plus":
        X_mhbafs = np.minimum(X_bafs, 1.0 - X_bafs)
        inits_maxK, inits_diag = init_hmm_kmeans_plus_plus(
            X_mhbafs=X_mhbafs,
            X_rdrs=X_hmm_rdrs,
            rdr_vars=rdr_vars0,
            K=maxK,
            random_state=seed,
            restarts=restarts,
            n_local_trials=n_local_trials,
            baf_eps=baf_eps,
        )
    else:
        inits_maxK, inits_diag = init_hmm_cna_plus_plus(
            X_hmm_rdrs,
            X_bafs,
            X_alphas,
            X_betas,
            X_totals,
            baf_taus0,
            rdr_vars=rdr_vars0,
            K=maxK,
            random_state=seed,
            restarts=restarts,
            n_local_trials=n_local_trials,
            log_rdr=log_rdr,
            baf_eps=baf_eps,
            bal_margin=args["bal_margin"],
            collect_diag=DEBUG,
        )
    plot_2d_inits(
        X_rdrs,
        X_bafs,
        inits_maxK,
        ntumor_samples,
        maxK,
        os.path.join(plot_dir, "hmm_init.pdf"),
        baf_taus=baf_taus0,
        log_rdr=log_rdr,
        bbs=bbs,
        chrom_sizes=chrom_sizes,
        init_method=init_method,
        sample_names=tumor_samples,
    )

    ##################################################
    score_records = []
    elbo_data = []  # list of (K, all_elbo_traces, best_it)

    sorted_inits = sorted(inits_maxK.items(), key=lambda x: x[1][-1], reverse=False)
    inits_run = dict(sorted_inits[:top_restarts])

    if DEBUG and inits_diag:
        init_diag_dir = os.path.join(plot_dir, "init_diag")
        os.makedirs(init_diag_dir, exist_ok=True)
        for it in inits_run:
            diag = inits_diag[it]
            plot_init_sampling_probs(
                X_hmm_rdrs,
                X_bafs,
                diag["probs_history"],
                diag["centroids_history"],
                init_diag_dir,
                name=f"restart{it}",
                log_rdr=log_rdr,
                bin_info=bbs,
                chrom_sizes=chrom_sizes,
                selected_bins_history=diag["selected_bins_history"],
                candidates_history=diag.get("candidates_history"),
                final_baf_means=diag["final_baf_means"],
                final_rdr_means=diag["final_rdr_means"],
            )
    logging.info(
        "use top %d/%d restarts to run HMM",
        top_restarts,
        len(inits_maxK),
    )

    for K in range(minK, maxK + 1):
        logging.info("==================================================")
        logging.info(
            f"running HMM on K={K}, {len(inits_run)} restarts ({training_method})"
        )
        log_transmat0 = np.log(make_transmat(1 - diag_t, K))
        run_fn = (
            run_baum_welch if training_method == "baum_welch" else run_viterbi_training
        )

        t0 = time.perf_counter()
        best_ll = -np.inf
        best_it = 0
        best_sol = None
        all_elbo_traces = {}
        for it, (baf_means_it, rdr_means_it, rdr_vars_it, _) in inits_run.items():
            sol = run_fn(
                K,
                X_hmm_rdrs,
                X_alphas,
                X_betas,
                X_totals,
                X_lengths,
                log_switchprobs,
                log_stayprobs,
                log_transmat0,
                np.ascontiguousarray(rdr_means_it[:K], dtype=np.float64),
                np.ascontiguousarray(rdr_vars_it[:K], dtype=np.float64),
                np.ascontiguousarray(baf_means_it[:K], dtype=np.float64),
                baf_taus0.copy(),
                X_rdrs_orig=X_rdrs,
                X_totals_orig=X_totals,
                n_iter=n_iter,
                min_covar=min_covar,
                tau_iters=tau_iters,
                min_tau=min_tau,
                max_tau=max_tau,
                baf_eps=baf_eps,
                log_rdr=log_rdr,
                restart_id=it,
                ig_alpha=ig_alpha,
                ig_beta=ig_beta,
                baf_k_start=baf_k_start,
            )
            model_ll = sol["model_ll"]
            obj_ll = sol["obj_ll"]
            all_elbo_traces[it] = sol["elbo_trace"]
            logging.info(
                f"K={K} restart {it}: model_ll={model_ll:.6f} obj_ll={obj_ll:.6f}"
            )
            if score_method == "bic":
                score = score_BIC(model_ll, K, ntumor_samples, nbbs)
            else:
                score = score_ICL(
                    sol["cluster_posts"],
                    model_ll,
                    K,
                    ntumor_samples,
                    nbbs,
                )
            score_records.append(
                {"K": K, "restart_it": it, "ll": model_ll, score_method: score}
            )
            if model_ll > best_ll:
                best_ll = model_ll
                best_it = it
                best_sol = sol
        t_elapsed = time.perf_counter() - t0
        logging.info(
            f"K={K} best restart: it={best_it}, model_ll={best_ll:.6f}, time={t_elapsed:.1f}s"
        )
        elbo_data.append((K, all_elbo_traces, best_it))

        # Save EM parameter trace for best restart
        trace_dir = os.path.join(out_dir, "traces")
        os.makedirs(trace_dir, exist_ok=True)
        np.savez_compressed(
            os.path.join(trace_dir, f"K{K}.em_trace.npz"),
            elbo_trace=np.array(best_sol["elbo_trace"]),
            rdr_means=best_sol["trace_rdr_means"],
            rdr_vars=best_sol["trace_rdr_vars"],
            baf_means=best_sol["trace_baf_means"],
            baf_taus=best_sol["trace_baf_taus"],
        )

        k_labels, k_phases, decode_ll = decode_hmm(
            best_sol,
            decode_method,
            X_lengths,
            log_switchprobs,
            log_stayprobs,
            best_sol.get("log_transmat", log_transmat0),
        )
        logging.info(
            f"K={K} decode path loglik={decode_ll:.6f} "
            f"(model_ll={best_ll:.6f}, diff={best_ll - decode_ll:.6f})"
        )
        k_betas_phased = (
            X_alphas * (1 - k_phases[:, None]) + X_betas * k_phases[:, None]
        )
        k_bafs = k_betas_phased / X_totals
        k_cids = np.unique(k_labels)
        k_rdr_means = best_sol["RDR_means"][k_cids]
        k_rdr_vars = best_sol["RDR_vars"][k_cids]
        k_baf_means = best_sol["BAF_means"][k_cids]
        k_baf_taus = best_sol["BAF_taus"]

        # mhBAF fold: flip BAF means and phases for clusters with BAF > 0.5
        if not args["skip_mhbafs"]:
            for ci, c in enumerate(k_cids):
                if np.mean(k_baf_means[ci]) > 0.5:
                    k_baf_means[ci] = 1.0 - k_baf_means[ci]
                    mask = k_labels == c
                    k_phases[mask] = 1 - k_phases[mask]
                    k_betas_phased[mask] = X_totals[mask] - k_betas_phased[mask]
                    k_bafs[mask] = k_betas_phased[mask] / X_totals[mask]

        balanced_ids = label_balanced_clusters(
            k_cids,
            k_labels,
            X_betas,
            X_totals,
            k_baf_means,
            k_baf_taus,
            alpha=args["bal_lrt_alpha"],
            margin=args["bal_margin"],
            baf_eps=baf_eps,
        )
        for ci, c in enumerate(k_cids):
            if c in balanced_ids:
                mask = k_labels == c
                k_baf_means[ci] = 0.5
                k_bafs[mask] = X_bafs[mask]
                k_betas_phased[mask] = X_betas[mask]
        logging.info(f"K={K} balanced clusters: {sorted(int(x) for x in balanced_ids)}")

        filtered_ids = filter_clusters(
            k_cids,
            k_labels,
            X_rdrs,
            k_bafs,
            k_rdr_means if not log_rdr else np.exp(k_rdr_means),
            k_baf_means,
            fstd=args["filter_std"],
            min_nbins=args["min_nbins"],
            ub_nbins=args["ub_nbins"],
        )
        if filtered_ids:
            logging.info(
                f"K={K} filtered clusters: {sorted(int(x) for x in filtered_ids)}"
            )

        if log_rdr:
            k_rdr_means_nat = np.exp(k_rdr_means)
            k_rdr_vars_nat = np.exp(2.0 * k_rdr_means) * k_rdr_vars
        else:
            k_rdr_means_nat = k_rdr_means
            k_rdr_vars_nat = k_rdr_vars

        for ci, c in enumerate(k_cids):
            mask = k_labels == c
            inf_rdr = k_rdr_means_nat[ci]
            inf_mhbaf = np.minimum(k_baf_means[ci], 1.0 - k_baf_means[ci])
            logging.info(
                f"K={K} cluster {c:2d} (n={mask.sum():5d}): "
                f"RDR={np.round(inf_rdr, 3)} | "
                f"mhBAF={np.round(inf_mhbaf, 3)}"
            )

        plot_rdr_baf(
            tumor_samples,
            bbs,
            k_bafs,
            X_rdrs,
            genome_size,
            cluster_labels=k_labels,
            expected_rdrs=k_rdr_means_nat,
            expected_bafs=k_baf_means,
            unique_labels=k_cids,
            label_clone=False,
            xlab="mhBAF",
            ylab="RDR",
            out_dir=plot_dir,
            out_prefix=f"K{K}_",
            dpi=100,
            rdr_means=k_rdr_means,
            rdr_vars=k_rdr_vars,
            baf_taus=k_baf_taus,
            log_rdr=log_rdr,
            filtered_ids=filtered_ids,
        )

        bbs["PHASE"] = k_phases
        bbs["PHASE_POSTS"] = best_sol["phase_posts"][:, 1]
        bbs[["#CHR", "START", "END", "PHASE", "PHASE_POSTS", "switchprobs"]].to_csv(
            os.path.join(label_dir, f"bulk{K}.bb.phased.tsv.gz"),
            sep="\t",
            header=True,
            index=False,
        )

        bbcs["CLUSTER"] = np.repeat(k_labels, ntumor_samples)
        bbcs["BAF"] = k_bafs.ravel()
        bbcs["BETA"] = k_betas_phased.ravel().astype(int)
        bbcs["ALPHA"] = (X_totals - k_betas_phased).ravel().astype(int)
        k_baf_ses = compute_baf_se(k_labels, k_bafs, k_cids)
        k_rdr_ses = compute_rdr_se(k_labels, k_rdr_vars_nat, k_cids)
        k_segs = mat2segs(
            bbcs,
            tumor_samples,
            k_baf_means,
            k_baf_taus,
            k_baf_ses,
            k_rdr_means_nat,
            k_rdr_vars_nat,
            k_rdr_ses,
            k_cids,
        )
        k_segs["is_balanced"] = k_segs["#ID"].isin(balanced_ids)
        k_segs["is_filtered"] = k_segs["#ID"].isin(filtered_ids)
        bbcs.to_csv(
            os.path.join(label_dir, f"bulk{K}.bbc"),
            sep="\t",
            header=True,
            index=False,
        )
        k_segs.to_csv(
            os.path.join(label_dir, f"bulk{K}.seg"),
            sep="\t",
            header=True,
            index=False,
        )

    plot_elbo_traces(elbo_data, os.path.join(plot_dir, "elbo_traces.pdf"))

    scores_df = pd.DataFrame(score_records)
    best_idx = scores_df[score_method].idxmin()
    best_K = int(scores_df.loc[best_idx, "K"])
    best_score = scores_df.loc[best_idx, score_method]
    logging.info(f"model selection: best K={best_K} {score_method}={best_score:.4f}")
    scores_df.to_csv(os.path.join(out_dir, "model_scores.tsv"), sep="\t", index=False)
    plot_score(scores_df, score_method, os.path.join(plot_dir, "model_scores.png"))

    ##################################################
    # copy best-K results to top-level output
    for suffix in ["bbc", "seg"]:
        shutil.copy2(
            os.path.join(label_dir, f"bulk{best_K}.{suffix}"),
            os.path.join(out_dir, f"bulk.{suffix}"),
        )
    shutil.copy2(
        os.path.join(label_dir, f"bulk{best_K}.bb.phased.tsv.gz"),
        os.path.join(out_dir, "bb.phased.tsv.gz"),
    )
    shutil.copy2(
        os.path.join(plot_dir, f"K{best_K}.pdf"),
        os.path.join(out_dir, f"bulk.K{best_K}.pdf"),
    )

    _log_done("cluster-bins")
    return
