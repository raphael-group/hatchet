import os
import logging
import shutil
import argparse

import pandas as pd

from hatchet.utils import (
    add_file_logging,
    log_arguments,
    read_bbc_file,
    setup_logging,
)
from hatchet.compute_cn.compute_cn_utils import (
    store_gammas,
    store_solve_input,
    store_instance_tofile,
    build_data,
    compute_fractional_cn,
    load_pool_from_disk,
    plot_pareto_pdf,
    run_plot_cn,
    segmentation,
)
from hatchet.compute_cn.model_select import (
    model_selection_ploidy,
    model_select_elbow_from_regularization,
)
from hatchet.compute_cn.scaling import get_scaling_factor
from hatchet.hatchet_parser import parse_arguments_compute_cn, add_arguments_compute_cn
from hatchet.compute_cn.solve.datatypes import SolverParams, SolverInputs
from hatchet.compute_cn.solve.inference import (
    run_full_ilp,
    run_coordinate_descent,
    run_coordinate_descent_cnt,
)
from hatchet.plot.plot_pool import plot_pool_cnp


def run(args=None):
    logging.info("run hatchet compute cn")
    if isinstance(args, argparse.Namespace):
        args = vars(args)
    args = parse_arguments_compute_cn(args)

    bbc_file = args["bbc"]
    seg_file = args["seg"]
    out_dir = args["result_dir"]
    os.makedirs(out_dir, exist_ok=True)
    add_file_logging(out_dir, "compute-cn")
    log_arguments(args)
    plot_dir = os.path.join(out_dir, "plots")
    sols_dir = os.path.join(out_dir, "sols")
    os.makedirs(plot_dir, exist_ok=True)
    os.makedirs(sols_dir, exist_ok=True)

    bbcs = read_bbc_file(bbc_file)
    segs = pd.read_table(seg_file, sep="\t")

    samples = sorted(bbcs["SAMPLE"].unique().tolist())

    # Remove clusters marked as filtered by cluster-bins
    filtered_ids = segs.loc[segs["is_filtered"], "#ID"].unique().tolist()
    if filtered_ids:
        logging.info(f"Excluding filtered clusters from seg: {filtered_ids}")
        segs = segs[~segs["is_filtered"]].reset_index(drop=True)
        bbcs = bbcs[~bbcs["CLUSTER"].isin(filtered_ids)].reset_index(drop=True)

    scaling, balanced_clusters = get_scaling_factor(
        samples,
        segs,
        bbcs,
        fix_cn_dip=args["fix_cn_dip"],
        fix_cn_tet=args["fix_cn_tet"],
        maxcn=args["diploidcmax"],
        maxcn_wgd=args["tetraploidcmax"],
    )
    gamma_outfile = os.path.join(out_dir, "gammas.tsv")
    store_gammas(gamma_outfile, scaling, samples)

    solve_mode = args["mode"]
    input_data = build_data(bbcs, segs, segment=(solve_mode == "cnt_cd"))
    minClone = args["minClone"]
    maxClone = args["maxClone"] + 1

    run_ploidy = {"diploid": args["diploid"], "tetraploid": args["tetraploid"]}
    if not run_ploidy["diploid"] and not run_ploidy["tetraploid"]:
        run_ploidy["diploid"] = True
        run_ploidy["tetraploid"] = True

    # 3. Run solves
    whole_pool = {}
    chosen_sols = {}
    model_selection_df = []
    for ploidy, run_it in run_ploidy.items():
        if not run_it or scaling[ploidy] is None:
            continue
        logging.info(f"running {ploidy} with n={minClone}..{maxClone}")
        whole_pool[ploidy] = {}
        chosen_sols[ploidy] = {}
        gammas = scaling[ploidy]["gammas"]
        clonals = scaling[ploidy]["clonal"]
        purities = scaling[ploidy]["purities"]
        logging.info(f"{ploidy} clonal CN: {clonals}")
        for sample, gamma in gammas.items():
            logging.info(f"  {sample}\tgamma={gamma}")

        fcn_data = compute_fractional_cn(
            input_data,
            gammas,
            alpha=args["fcn_ci_alpha"],
            min_ci_margin=args["min_ci_margin"],
        )
        store_solve_input(
            os.path.join(out_dir, "sols", f"solver_input.{ploidy}.tsv"),
            fcn_data,
        )

        for n in range(minClone, maxClone):
            out_bbc = os.path.join(out_dir, f"results.{ploidy}.n{n}.bbc.ucn.tsv")
            out_seg = os.path.join(out_dir, f"results.{ploidy}.n{n}.seg.ucn.tsv")
            sol_dir = os.path.join(out_dir, f"sols/{ploidy}_n{n}")
            if (
                not args["force"]
                and os.path.exists(out_bbc)
                and os.path.exists(out_seg)
                and os.path.isdir(sol_dir)
            ):
                logging.info(
                    f"skip {ploidy} n={n}: results already exist (use --force to re-solve)"
                )
                pool_instances = load_pool_from_disk(
                    sol_dir, fcn_data["cluster_ids"], fcn_data["sample_ids"]
                )
            else:
                os.makedirs(sol_dir, exist_ok=True)
                pool_instances = solve(
                    n,
                    clonals,
                    args,
                    ploidy,
                    bbcs,
                    fcn_data,
                    purities,
                    sol_dir,
                    solve_mode=solve_mode,
                    balanced_clusters=balanced_clusters,
                )
            whole_pool[ploidy][n] = pool_instances

            selected_id, sel_df = model_select_elbow_from_regularization(pool_instances)
            best_sol = pool_instances[selected_id]
            chosen_sols[ploidy][n] = best_sol
            logging.info(
                f"{ploidy} n={n} selected={selected_id} "
                f"fit_loss={best_sol['fit_loss']:.4f} imf={best_sol['imf_obj']:.4f}"
            )
            sel_df["ploidy"] = ploidy
            sel_df["n_clones"] = n
            model_selection_df.append(sel_df)

            cn_segs = {}
            for sol_id, sol in pool_instances.items():
                bbc_out = out_bbc if sol_id == selected_id else None
                seg_out = out_seg if sol_id == selected_id else None
                cn_segs[sol_id] = segmentation(
                    sol["cA"],
                    sol["cB"],
                    sol["u"],
                    fcn_data,
                    bbcs=bbcs,
                    region_file=args["region_bed"],
                    bbc_out_file=bbc_out,
                    seg_out_file=seg_out,
                )

            run_plot_cn(
                args,
                out_bbc,
                out_seg,
                gamma_outfile,
                os.path.join(plot_dir, f"{ploidy}_n{n}"),
                ploidy,
            )

            plot_pool_cnp(
                pool_instances,
                args["region_bed"],
                os.path.join(plot_dir, f"{ploidy}_n{n}_pool"),
                sel_df=sel_df,
                segs=cn_segs,
                title=f"{ploidy} n={n}",
                solve_mode=solve_mode,
                sample_names=fcn_data["sample_ids"],
            )

    summary_df = (
        pd.concat(model_selection_df, ignore_index=True)
        if model_selection_df
        else pd.DataFrame()
    )
    summary_path = os.path.join(out_dir, "summary.tsv")
    summary_df.to_csv(summary_path, sep="\t", index=False)
    logging.info(f"wrote {summary_path} ({len(summary_df)} solutions)")

    best_ploidy, best_n, chosen_n, elbow_fig = model_selection_ploidy(
        chosen_sols,
        out_dir,
        scaling,
        segs,
        method=args["model_select"],
    )
    plot_pareto_pdf(summary_df, plot_dir, args["reg_term"], elbow_fig)

    # Write chosen per-ploidy
    for ploidy, n in chosen_n.items():
        shutil.copy2(
            os.path.join(out_dir, f"results.{ploidy}.n{n}.bbc.ucn.tsv"),
            os.path.join(out_dir, f"chosen.{ploidy}.bbc.ucn"),
        )
        shutil.copy2(
            os.path.join(out_dir, f"results.{ploidy}.n{n}.seg.ucn.tsv"),
            os.path.join(out_dir, f"chosen.{ploidy}.seg.ucn"),
        )
        logging.info(
            f"chosen {ploidy} n={n}: {os.path.join(out_dir, f'chosen.{ploidy}.bbc.ucn')}"
        )

    # Write best (across ploidies)
    shutil.copy2(
        os.path.join(out_dir, f"chosen.{best_ploidy}.bbc.ucn"),
        os.path.join(out_dir, "best.bbc.ucn"),
    )
    shutil.copy2(
        os.path.join(out_dir, f"chosen.{best_ploidy}.seg.ucn"),
        os.path.join(out_dir, "best.seg.ucn"),
    )
    logging.info(f"model-selected: {best_ploidy} n={best_n}")


def solve(
    n: int,
    clonal: dict,
    args: dict,
    ploidy: str,
    bbcs: pd.DataFrame,
    input_data: dict,
    purities: dict,
    sol_dir: str,
    solve_mode="ilp",
    balanced_clusters=None,
):
    """Solve for allele-specific integer copy numbers and clone proportions.

    Runs coordinate descent (CD), integer linear programming (ILP), or both
    (CD warm-starting ILP) over a regularization path, selects the best
    instance via Pareto-elbow model selection, and writes BBC/SEG UCN output.

    Returns (objective, IMF-objective, pool_dict) of the selected solution.
    """
    logging.info(f"running {ploidy} with n={n}")
    cn_max = {"diploid": args["diploidcmax"], "tetraploid": args["tetraploidcmax"]}[
        ploidy
    ]
    base = {"diploid": 1, "tetraploid": 2}[ploidy]
    if args["purities"] is not None:
        purities = args["purities"]
        logging.info(f"purities overridden by user: {purities}")
    elif purities is not None:
        logging.info(f"purities: {purities}")

    reg_term = args["reg_term"]
    reg_steps = args["reg_steps"]
    solver_type = args["solver"]
    timelimit = args["timelimit"]
    ampdel = not args["no_ampdel"]

    cd_instances = None
    pool_instances = {}
    u0_tsv_path = os.path.join(sol_dir, "u0_seeds.tsv") if sol_dir is not None else None
    cd_run_kwargs = dict(
        solver_type=solver_type,
        max_iters=args["cd_niters"],
        max_convergence_iters=args["cd_convergence_iters"],
        n_seed=args["cd_nseeds"],
        j=args["cd_njobs"],
        random_seed=args["cd_seed"],
        timelimit=timelimit,
        u0_tsv_path=u0_tsv_path,
    )

    # Build SolverParams shared by both CD and ILP
    weights = input_data["weights"]
    cluster_ids = input_data["cluster_ids"]
    sample_ids = input_data["sample_ids"]
    nbins = input_data["nbins"]

    fixed_rows = set()
    for _m, cid in enumerate(cluster_ids):
        if cid in clonal:
            fixed_rows.add(_m)
    free_rows_list = [_m for _m in range(len(cluster_ids)) if _m not in fixed_rows]

    inputs = SolverInputs(
        f_a=input_data["fa"],
        f_b=input_data["fb"],
        w=weights,
        cluster_ids=cluster_ids,
        sample_ids=sample_ids,
        copy_numbers=clonal,
        free_rows=free_rows_list,
        fixed_rows=fixed_rows,
        purities=purities,
        balanced_clusters=balanced_clusters,
        fa_lo=input_data["fa_lo"],
        fa_hi=input_data["fa_hi"],
        fb_lo=input_data["fb_lo"],
        fb_hi=input_data["fb_hi"],
        nbins=nbins,
        chr_boundaries=input_data.get("chr_boundaries"),
    )
    params = SolverParams(
        n=n,
        cn_max=cn_max,
        base=base,
        ampdel=ampdel,
        minprop=args["min_prop"],
        max_ncns_seg=args["num_cnstates"],
        tol=args["tol"],
        zero_cn_thres=args["zero_cn_thres"],
        reg_name=reg_term if reg_term is not None else "RAW",
        obj_type=args["obj_type"],
        eps_fit=args["eps_fit"],
    )

    if solve_mode == "cnt_cd":
        pool_instances = run_coordinate_descent_cnt(
            params=params,
            inputs=inputs,
            solver_type=solver_type,
            max_iters=args["cd_niters"],
            max_convergence_iters=args["cd_convergence_iters"],
            n_seed=args["cd_nseeds"],
            j=args["cd_njobs"],
            cd_tol=args["cd_tol"],
            random_seed=args["cd_seed"],
            timelimit=timelimit,
            tree_file=args["tree_file"],
            u_dir_alpha=args["u_dir_alpha"],
        )

    elif solve_mode in ("cd", "both"):
        cd_instances = run_coordinate_descent(
            params=params,
            inputs=inputs,
            reg_steps=reg_steps,
            reg_bound=args["reg_bound"],
            u_init_method=args["u_init"],
            u_dir_alpha=args["u_dir_alpha"],
            solver_threads=args["solver_threads"],
            cd_tol=args["cd_tol"],
            **cd_run_kwargs,
        )
        pool_instances = cd_instances

    if solve_mode in ("ilp", "both"):
        warm_cA = warm_cB = None
        if solve_mode == "both":
            best_cd = min(cd_instances.values(), key=lambda s: s["fit_loss"])
            logging.info(
                f"use CD local opt with obj={best_cd['fit_loss']:.4f} to initialize ILP model"
            )
            warm_cA, warm_cB = best_cd["cA"], best_cd["cB"]

        pool_instances = run_full_ilp(
            params=params,
            inputs=inputs,
            reg_steps=reg_steps,
            reg_bound=args["reg_bound"],
            solver_type=solver_type,
            timelimit=timelimit,
            warm_start_cA=warm_cA,
            warm_start_cB=warm_cB,
        )

    store_instance_tofile(pool_instances, input_data, sol_dir, solve_mode)
    return pool_instances


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="HATCHet compute-cn",
        description="solve copy-numbers",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    add_arguments_compute_cn(parser)
    args = parser.parse_args()
    setup_logging(args)
    run(args)
