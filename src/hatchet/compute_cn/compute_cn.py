import os
import logging
import shutil
import argparse

import numpy as np
import pandas as pd

from hatchet.utils import *
from hatchet.compute_cn.compute_cn_utils import *
from hatchet.compute_cn.scaling import get_scaling_factor
from hatchet.compute_cn.model_select import *
from hatchet.hatchet_parser import parse_arguments_compute_cn
from hatchet.compute_cn.solve.utils import (
    store_solve_input,
    store_instance_tofile,
    model_selection_instance,
    compute_individual_objs,
    filter_non_pareto,
)
from hatchet.compute_cn.solve.ilp_subset import ILPSubset
from hatchet.compute_cn.solve.cd import CoordinateDescent
from hatchet.plot.plot_cnp_panel import plot_pool_cnp


def run(args=None):
    logging.info("run hatchet compute cn")
    if isinstance(args, argparse.Namespace):
        args = vars(args)

    bbc_file = args["bbc"]
    seg_file = args["seg"]

    out_dir = args["result_dir"]
    os.makedirs(out_dir, exist_ok=True)
    add_file_logging(out_dir, "compute-cn")
    plot_dir = os.path.join(out_dir, "plots")
    sols_dir = os.path.join(out_dir, "sols")
    os.makedirs(plot_dir, exist_ok=True)
    os.makedirs(sols_dir, exist_ok=True)

    logging.info("load arguments")
    bbcs = read_bbc_file(bbc_file)
    segs = pd.read_table(seg_file, sep="\t")

    samples = sorted(bbcs["SAMPLE"].unique().tolist())
    clusters = sorted(segs["#ID"].unique().tolist())
    # TODO: move cluster filtering into cluster-bins step
    if args["filter_cluster"]:
        good_clusters, bad_clusters = filtering(
            bbc=bbcs,
            seg=segs,
            samples=samples,
            clusters=clusters,
            fstd=args["filter_std"],
            min_nbins=args["min_nbins"],
            ub_nbins=args["ub_nbins"],
        )
        if len(bad_clusters) > 0:
            segs = segs[segs["#ID"].isin(good_clusters)].reset_index(drop=True)
            bbcs = bbcs[bbcs["CLUSTER"].isin(good_clusters)].reset_index(drop=True)
            fseg_path = os.path.join(out_dir, "bulk.good.seg")
            segs.to_csv(fseg_path, header=True, index=False, sep="\t")
            args["seg"] = fseg_path
            fbbc_path = os.path.join(out_dir, "bulk.good.bbc")
            bbcs.to_csv(fbbc_path, header=True, index=False, sep="\t")
            args["bbc"] = fbbc_path

    (
        s0,
        pair_noWGD,
        gammas_noWGD,
        purities_noWGD,
        pair_WGD,
        gammas_WGD,
        purities_WGD,
    ) = get_scaling_factor(
        samples,
        segs,
        bal_tost_alpha=args["bal_tost_alpha"],
        bal_tost_margin=args["bal_tost_margin"],
        tol_nstd=args["tol_nstd"],
        maxcn=args["diploidcmax"],
        maxcn_wgd=args["tetraploidcmax"],
    )
    gamma_outfile = os.path.join(out_dir, "gammas.tsv")
    with open(gamma_outfile, "w") as fd:
        for sample in samples:
            gamma_noWGD = gammas_noWGD.get(sample, 0)
            gamma_WGD = gammas_WGD.get(sample, 0) if gammas_WGD is not None else 0
            fd.write(f"{sample}\t{gamma_noWGD}\t{gamma_WGD}\n")

    if args.get("segment", False):
        seg_data = build_segment_data(bbcs, segs)
    else:
        seg_data = build_cluster_data(segs)
    rdr = seg_data["rdr"]
    baf = seg_data["baf"]
    rdr_se = seg_data["rdr_se"]
    baf_se = seg_data["baf_se"]
    nbins = seg_data["nbins"]
    weights = seg_data["weights"]
    cluster_ids = rdr.index.tolist()
    sample_ids = rdr.columns.tolist()

    minClone = args["minClone"]
    maxClone = args["maxClone"] + 1

    run_diploid = args["diploid"]
    run_tetraploid = args["tetraploid"]
    if not run_diploid and not run_tetraploid:
        run_diploid = run_tetraploid = True

    diploid_sols = {}
    if run_diploid:
        if pair_noWGD is not None:
            logging.info(f"Inferred (1,1) balanced cluster={s0}")
            (s, z, (sa, sb), (za, zb)) = pair_noWGD
            logging.info(f"Inferred clonal pair: {s}:({sa},{sb}), {z}:({za},{zb})")
            clonal_dip = {s: (sa, sb), z: (za, zb)}
            logging.info("Inferred diploid RD scaling factor gamma per sample:")
            for sample, gamma in gammas_noWGD.items():
                logging.info(f"{sample}\tgamma={gamma}")
            gammas_dip = pd.Series(gammas_noWGD).sort_index()
            fcn_dip = compute_fractional_cn(
                rdr,
                baf,
                rdr_se,
                baf_se,
                gammas_dip,
                alpha=args.get("fcn_ci_alpha", 0.05),
            )
            store_solve_input(
                os.path.join(out_dir, "sols", "diploid_input.tsv"),
                fcn_dip,
                weights,
                nbins,
            )
            for n in range(minClone, maxClone):
                logging.info(f"running diploid with n={n}")
                obj, imf_obj, pool = solve(
                    n,
                    clonal_dip,
                    args,
                    "diploid",
                    out_dir,
                    plot_dir,
                    bbcs,
                    fcn_dip,
                    weights,
                    cluster_ids,
                    sample_ids,
                    purities_noWGD,
                    nbins,
                    args["mode"],
                    args["verbosity"],
                )
                diploid_sols[n] = (obj, imf_obj)
                logging.info(f"diploid n={n} objective={obj} imf-objective={imf_obj}")
                out_bbc = os.path.join(out_dir, f"results.diploid.n{n}.bbc.ucn.tsv")
                out_seg = os.path.join(out_dir, f"results.diploid.n{n}.seg.ucn.tsv")
                run_plot_cn(
                    args,
                    out_bbc,
                    out_seg,
                    gamma_outfile,
                    os.path.join(plot_dir, f"diploid_n{n}"),
                    "diploid",
                )
                if pool:
                    pool_entries = [
                        (tag, seg_df, obj_, pareto, selected)
                        for tag, (seg_df, obj_, pareto, selected) in pool.items()
                    ]
                    plot_pool_cnp(
                        pool_entries,
                        args["region_bed"],
                        os.path.join(plot_dir, f"diploid_n{n}_pool_pareto.pdf"),
                        title=f"diploid n={n} pool solutions",
                    )
        else:
            logging.warning("run_diploid=True, but failed to infer clonal pair")

    tetraploid_sols = {}
    if run_tetraploid and pair_WGD is not None:
        (s, z, (sa, sb), (za, zb)) = pair_WGD
        logging.info(f"Inferred clonal pair: {s}:({sa},{sb}), {z}:({za},{zb})")
        logging.info("Inferred tetraploid RD scaling factor gamma per sample:")
        for sample, gamma in gammas_WGD.items():
            logging.info(f"{sample}\tgamma={gamma}")
        clonal_tet = {s: (sa, sb), z: (za, zb)}
        gammas_tet = pd.Series(gammas_WGD).sort_index()
        fcn_tet = compute_fractional_cn(
            rdr,
            baf,
            rdr_se,
            baf_se,
            gammas_tet,
            alpha=args.get("fcn_ci_alpha", 0.05),
        )
        store_solve_input(
            os.path.join(out_dir, "sols", "tetraploid_input.tsv"),
            fcn_tet,
            weights,
            nbins,
        )
        for n in range(minClone, maxClone):
            logging.info(f"running tetraploid with n={n}")
            obj, imf_obj, pool = solve(
                n,
                clonal_tet,
                args,
                "tetraploid",
                out_dir,
                plot_dir,
                bbcs,
                fcn_tet,
                weights,
                cluster_ids,
                sample_ids,
                purities_WGD,
                nbins,
                args["mode"],
                args["verbosity"],
            )
            tetraploid_sols[n] = (obj, imf_obj)
            logging.info(f"tetraploid n={n} objective={obj} imf-objective={imf_obj}")
            out_bbc = os.path.join(out_dir, f"results.tetraploid.n{n}.bbc.ucn.tsv")
            out_seg = os.path.join(out_dir, f"results.tetraploid.n{n}.seg.ucn.tsv")
            run_plot_cn(
                args,
                out_bbc,
                out_seg,
                gamma_outfile,
                os.path.join(plot_dir, f"tetraploid_n{n}"),
                "tetraploid",
            )
            if pool:
                pool_entries = [
                    (tag, seg_df, obj_, pareto)
                    for tag, (seg_df, obj_, pareto) in pool.items()
                ]
                plot_pool_cnp(
                    pool_entries,
                    args["region_bed"],
                    os.path.join(plot_dir, f"tetraploid_n{n}_pool_pareto.pdf"),
                    title=f"tetraploid n={n} pool solutions",
                )
    elif run_tetraploid:
        logging.warning("run_tetraploid=True, but failed to infer clonal pair")

    if len(diploid_sols) == 0 and len(tetraploid_sols) == 0:
        raise ValueError("No solutions found for either noWGD or WGD case, exit..")

    n_dip, n_tet, best_type = model_selection(
        diploid_sols,
        tetraploid_sols,
        out_dir,
        gammas_noWGD,
        gammas_WGD,
        segs,
        args["verbosity"],
    )

    if n_dip > 0:
        shutil.copy2(
            os.path.join(out_dir, f"results.diploid.n{n_dip}.bbc.ucn.tsv"),
            os.path.join(out_dir, "chosen.diploid.bbc.ucn"),
        )
        shutil.copy2(
            os.path.join(out_dir, f"results.diploid.n{n_dip}.seg.ucn.tsv"),
            os.path.join(out_dir, "chosen.diploid.seg.ucn"),
        )
        logging.info(
            f"chosen diploid n={n_dip}: {os.path.join(out_dir, 'chosen.diploid.bbc.ucn')}"
        )
        logging.info(
            f"chosen diploid plots: {os.path.join(plot_dir, f'diploid_n{n_dip}')}"
        )

    if n_tet > 0:
        shutil.copy2(
            os.path.join(out_dir, f"results.tetraploid.n{n_tet}.bbc.ucn.tsv"),
            os.path.join(out_dir, "chosen.tetraploid.bbc.ucn"),
        )
        shutil.copy2(
            os.path.join(out_dir, f"results.tetraploid.n{n_tet}.seg.ucn.tsv"),
            os.path.join(out_dir, "chosen.tetraploid.seg.ucn"),
        )
        logging.info(
            f"chosen tetraploid n={n_tet}: {os.path.join(out_dir, 'chosen.tetraploid.bbc.ucn')}"
        )
        logging.info(
            f"chosen tetraploid plots: {os.path.join(plot_dir, f'tetraploid_n{n_tet}')}"
        )

    if best_type is not None:
        shutil.copy2(
            os.path.join(out_dir, f"chosen.{best_type}.bbc.ucn"),
            os.path.join(out_dir, "best.bbc.ucn"),
        )
        shutil.copy2(
            os.path.join(out_dir, f"chosen.{best_type}.seg.ucn"),
            os.path.join(out_dir, "best.seg.ucn"),
        )
        best_n = {"diploid": n_dip, "tetraploid": n_tet}[best_type]
        logging.info(
            f"model-selected result ({best_type}): {os.path.join(out_dir, 'best.bbc.ucn')}"
        )
        logging.info(f"best plots: {os.path.join(plot_dir, f'{best_type}_n{best_n}')}")

    return


def solve(
    n: int,
    clonal: dict,
    args: dict,
    ploidy: str,
    out_dir: str,
    plot_dir: str,
    bbcs: pd.DataFrame,
    fcn_data: dict,
    weights: pd.Series,
    cluster_ids: list,
    sample_ids: list,
    purities: dict,
    nbins: pd.DataFrame,
    solve_mode="ilp",
    verbosity=0,
):
    """Solve for allele-specific integer copy numbers and clone proportions.

    Runs coordinate descent (CD), integer linear programming (ILP), or both
    (CD warm-starting ILP) over a regularization path, selects the best
    instance via Pareto-elbow model selection, and writes BBC/SEG UCN output.

    Returns (objective, IMF-objective, pool_dict) of the selected solution.
    """
    f_a = fcn_data["fa"]
    f_b = fcn_data["fb"]

    sol_dir = os.path.join(out_dir, f"sols/{ploidy}_n{n}")
    os.makedirs(sol_dir, exist_ok=True)

    out_bbc = os.path.join(out_dir, f"results.{ploidy}.n{n}.bbc.ucn.tsv")
    out_seg = os.path.join(out_dir, f"results.{ploidy}.n{n}.seg.ucn.tsv")

    # TODO: support user-defined fixed clonal states & fixed clone proportions
    copy_number_fixed = None

    cn_max = {"diploid": args["diploidcmax"], "tetraploid": args["tetraploidcmax"]}[
        ploidy
    ]
    ampdel = not args["no_ampdel"]
    base = {"diploid": 1, "tetraploid": 2}[ploidy]
    if args["purities"] is not None:
        purities = args["purities"]
        logging.info(f"purities overridden by user: {purities}")
    elif purities is not None:
        logging.info(f"purities: {purities}")

    reg_term = args["reg_term"]
    reg_steps = args["reg_steps"]
    reg_stepsize = args["reg_stepsize"]
    solver_type = args["solver"]
    verbose = verbosity >= 1
    timelimit = args["timelimit"]
    pool_size = args.get("pool_size", 1)
    pool_gap = args.get("pool_gap", None)

    cd_instances = None
    pool_instances = {}
    if solve_mode in ("cd", "both"):
        cd = CoordinateDescent(
            f_a=f_a,
            f_b=f_b,
            n=n,
            minprop=args["min_prop"],
            max_ncns_seg=args["num_cnstates"],
            cn_max=cn_max,
            w=weights,
            ampdel=ampdel,
            cn=clonal,
            purities=purities,
            copy_numbers_fixed=copy_number_fixed,
            reg_term=reg_term,
            reg_steps=reg_steps,
            reg_stepsize=reg_stepsize,
            base=base,
        )

        u0_tsv_path = (
            os.path.join(sol_dir, "u0_seeds.tsv") if sol_dir is not None else None
        )
        cd_instances = cd.run(
            solver_type=solver_type,
            max_iters=args["cd_niters"],
            max_convergence_iters=args["cd_convergence_iters"],
            n_seed=args["cd_nseeds"],
            j=args["cd_njobs"],
            random_seed=args["cd_seed"],
            timelimit=timelimit,
            u0_tsv_path=u0_tsv_path,
        )
        pool_instances = {k: [v] for k, v in cd_instances.items()}
        store_instance_tofile(
            pool_instances,
            f_a,
            f_b,
            sol_dir,
            "cd",
            n,
            fcn_data=fcn_data,
            nbins=nbins,
        )

    sol_instances = None
    if solve_mode in ("ilp", "both"):
        sol_instances = {}
        solver = ILPSubset(
            n,
            cn_max,
            max_ncns_seg=args["num_cnstates"],
            minprop=args["min_prop"],
            ampdel=ampdel,
            copy_numbers=clonal,
            f_a=f_a,
            f_b=f_b,
            w=weights,
            purities=purities,
            copy_numbers_fixed=copy_number_fixed,
            penalty_param=[reg_term if reg_term is not None else "RAW", 0.0],
            base=base,
        )
        solver.create_model(pprint=verbose)
        if solve_mode == "both":
            _, [obj_, cA_, cB_, _] = min(cd_instances.items(), key=lambda tp: tp[1][0])
            logging.info(f"use CD local opt with obj={obj_} to initialize ILP model")
            solver.hot_start(cA_, cB_)

        # DMRCA_SUM only penalises clones at index >= 2; with n <= 2 there are no
        # subclonal clones beyond the MRCA, so the regularisation path has no effect
        # and a single unregularised solve (i0=0, lambda=0) is sufficient.
        dmrca_no_effect = reg_term == "DMRCA_SUM" and n <= 2
        effective_reg_steps = 0 if dmrca_no_effect else reg_steps

        pool_instances = {}
        for i0 in range(0, effective_reg_steps + 1):
            logging.debug(f"running instance {i0}/{effective_reg_steps}")
            pparam = reg_stepsize * i0
            solver.model.pparam = pparam
            if i0 > 0:
                cA_, cB_ = sol_instances[0][1:3]
                solver.hot_start(cA_, cB_)
            sol_instances[pparam] = solver.run(
                solver_type=solver_type,
                timelimit=timelimit,
                pool_size=pool_size,
                pool_gap=pool_gap,
            )
            assert sol_instances[pparam] is not None, "optimization failed"

            pool_instances[pparam] = [sol_instances[pparam]]
            if pool_size > 1 and solver_type in ("gurobi", "gurobipy"):
                pool_sols = solver.get_pool_solutions(pool_size=pool_size)
                if pool_sols:
                    pool_instances[pparam].extend(pool_sols)

        store_instance_tofile(
            pool_instances,
            f_a,
            f_b,
            sol_dir,
            solve_mode,
            n,
            fcn_data=fcn_data,
            nbins=nbins,
        )

    pareto_img = os.path.join(
        plot_dir, f"{ploidy}_n{n}_pareto_curve.{solve_mode}.{reg_term}.png"
    )
    best_instance, imf_obj, selected_key = model_selection_instance(
        f_a,
        f_b,
        weights,
        pool_instances,
        reg_term,
        solve_mode,
        sol_dir,
        fcn_data,
        nbins,
        pareto_img=pareto_img,
    )

    assert best_instance is not None, f"no solution for {ploidy} and n={n}"
    [obj, cA, cB, u] = best_instance
    segmentation(
        cA,
        cB,
        u,
        cluster_ids,
        sample_ids,
        bbcs=bbcs,
        region_file=args["region_bed"],
        bbc_out_file=out_bbc,
        seg_out_file=out_seg,
    )

    all_pool = {}
    pool_objs = []
    pool_tags = []
    pool_keys = []
    for pparam, sols in pool_instances.items():
        for pidx, (pobj, pcA, pcB, pu) in enumerate(sols):
            tag = f"pool_p{pparam}_s{pidx}"
            seg_df = segmentation(
                pcA,
                pcB,
                pu,
                cluster_ids,
                sample_ids,
                bbcs=bbcs,
                region_file=args["region_bed"],
            )
            p_imf, p_reg = compute_individual_objs(
                reg_term, weights, f_a, f_b, pcA, pcB, pu
            )
            pool_objs.append([p_imf, p_reg])
            pool_tags.append(tag)
            pool_keys.append((pparam, pidx))
            all_pool[tag] = (seg_df, p_imf, False, False)

    is_pareto = filter_non_pareto(np.array(pool_objs))
    for tag, key, pareto in zip(pool_tags, pool_keys, is_pareto):
        seg_df_, obj_, _, _ = all_pool[tag]
        is_selected = key == selected_key
        all_pool[tag] = (seg_df_, obj_, bool(pareto), is_selected)

    return obj, imf_obj, all_pool


if __name__ == "__main__":
    args = parse_arguments_compute_cn()
    setup_logging(args)
    run(args)
