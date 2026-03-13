import os
import sys
import time
import logging
import shutil
import argparse
import numpy as np
import pandas as pd


from hatchet.utils import *
from hatchet.compute_cn.compute_cn_utils import *
from hatchet.compute_cn.compute_cn_utils import run_plot_cn
from hatchet.compute_cn.model_select import *
from hatchet.hatchet_parser import parse_arguments_compute_cn

from hatchet.compute_cn.solve import *


def run(args=None):
    logging.info("run hatchet compute cn")
    if isinstance(args, argparse.Namespace):
        args = vars(args)

    bbc_file = args["bbc"]
    seg_file = args["seg"]

    out_dir = args["result_dir"]
    # output files
    os.makedirs(out_dir, exist_ok=True)
    add_file_logging(out_dir, "compute-cn")
    plot_dir = os.path.join(out_dir, "plots")
    sols_dir = os.path.join(out_dir, "sols")
    os.makedirs(plot_dir, exist_ok=True)
    os.makedirs(sols_dir, exist_ok=True)

    bbcs = read_bbc_file(bbc_file)
    segs = pd.read_table(seg_file, sep="\t")

    samples = sorted(bbcs["SAMPLE"].unique().tolist())
    clusters = sorted(segs["#ID"].unique().tolist())
    # filter outlier clusters
    # TODO, move cluster filtering into cluster-bins step, account for usage as well.
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

    # infer balanced clusters
    segs = annotate_balanced_clusters(
        segs, args["balanced_baf_tol"], colname="balanced"
    )
    balanced_s = segs.loc[segs["balanced"], "#ID"].unique().tolist()
    unbalanced_z = segs.loc[~segs["balanced"], "#ID"].unique().tolist()

    # estimate RDR scaling factor
    s0, pair_noWGD, gammas_noWGD, purities_noWGD, pair_WGD, gammas_WGD, purities_WGD = (
        get_scaling_factor(
            samples,
            segs,
            balanced_s,
            unbalanced_z,
            args["toleranceRDR"],
            args["toleranceBAF"],
            args["diploidcmax"],
            args["tetraploidcmax"],
        )
    )
    gamma_outfile = os.path.join(out_dir, "gammas.tsv")
    with open(gamma_outfile, "w") as fd:
        for sample in samples:
            gamma_noWGD = gammas_noWGD.get(sample, 0)
            gamma_WGD = gammas_WGD.get(sample, 0)
            fd.write(f"{sample}\t{gamma_noWGD}\t{gamma_WGD}\n")

    segs_sorted = segs.sort_values(["#ID", "SAMPLE"])
    rdr = segs_sorted.pivot(index="#ID", columns="SAMPLE", values="RD")
    baf = segs_sorted.pivot(index="#ID", columns="SAMPLE", values="BAF")
    samples_sorted = sorted(segs["SAMPLE"].unique().tolist())
    bins = (
        segs.loc[segs["SAMPLE"] == samples_sorted[0]]
        .set_index("#ID")["LENGTH"]
        .sort_index()
    )
    for cid, length in bins.items():
        logging.info(f"#{cid}\t{length}")
    weights = 100 * bins / sum(bins)
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
            logging.info(f"Inferred diploid RD scaling factor gamma per sample:")
            for sample, gamma in gammas_noWGD.items():
                logging.info(f"{sample}\tgamma={gamma}")
            for n in range(minClone, maxClone):
                logging.info(f"running diploid with n={n}")
                (obj, imf_obj) = solve_wrapper(
                    n,
                    clonal_dip,
                    gammas_noWGD,
                    args,
                    "diploid",
                    out_dir,
                    bbcs,
                    rdr,
                    baf,
                    weights,
                    cluster_ids,
                    sample_ids,
                    purities_noWGD,
                    args["mode"],
                    args["verbosity"],
                )
                diploid_sols[n] = (obj, imf_obj)
                logging.info(f"diploid n={n} objective={obj} imf-objective={imf_obj}")
                run_plot_cn(args, out_dir, plot_dir, gamma_outfile, "diploid", n)
        else:
            logging.warn(f"run_diploid=True, but failed to infer clonal pair")

    tetraploid_sols = {}
    if run_tetraploid and pair_WGD != None:
        if pair_WGD is not None:
            (s, z, (sa, sb), (za, zb)) = pair_WGD
            logging.info(f"Inferred clonal pair: {s}:({sa},{sb}), {z}:({za},{zb})")
            logging.info("Inferred tetraploid RD scaling factor gamma per sample:")
            for sample, gamma in gammas_WGD.items():
                logging.info(f"{sample}\tgamma={gamma}")
            clonal_tet = {s: (sa, sb), z: (za, zb)}
            for n in range(minClone, maxClone):
                logging.info(f"running tetraploid with n={n}")
                (obj, imf_obj) = solve_wrapper(
                    n,
                    clonal_tet,
                    gammas_WGD,
                    args,
                    "tetraploid",
                    out_dir,
                    bbcs,
                    rdr,
                    baf,
                    weights,
                    cluster_ids,
                    sample_ids,
                    purities_WGD,
                    args["mode"],
                    args["verbosity"],
                )
                tetraploid_sols[n] = (obj, imf_obj)
                logging.info(f"tetraploid n={n} objective={obj} imf-objective={imf_obj}")
                run_plot_cn(args, out_dir, plot_dir, gamma_outfile, "tetraploid", n)
        else:
            logging.warn(f"run_tetraploid=True, but failed to infer clonal pair")
    if len(diploid_sols) == 0 and len(tetraploid_sols) == 0:
        raise ValueError("No solutions found for either noWGD or WGD case, exit..")

    # final model selection between diploid and tetraploid with varying n.
    n_dip, n_tet, best_type = model_selection(
        diploid_sols,
        tetraploid_sols,
        out_dir,
        gammas_noWGD,
        gammas_WGD,
        segs,
        args["verbosity"],
    )

    # save model selected result here
    if n_dip > 0:
        shutil.copy2(
            os.path.join(out_dir, f"results.diploid.n{n_dip}.bbc.ucn.tsv"),
            os.path.join(out_dir, "chosen.diploid.bbc.ucn"),
        )
        shutil.copy2(
            os.path.join(out_dir, f"results.diploid.n{n_dip}.seg.ucn.tsv"),
            os.path.join(out_dir, "chosen.diploid.seg.ucn"),
        )
        logging.info(f"chosen diploid n={n_dip}: {os.path.join(out_dir, 'chosen.diploid.bbc.ucn')}")
        logging.info(f"chosen diploid plots: {os.path.join(plot_dir, f'diploid_n{n_dip}')}")

    if n_tet > 0:
        shutil.copy2(
            os.path.join(out_dir, f"results.tetraploid.n{n_tet}.bbc.ucn.tsv"),
            os.path.join(out_dir, "chosen.tetraploid.bbc.ucn"),
        )
        shutil.copy2(
            os.path.join(out_dir, f"results.tetraploid.n{n_tet}.seg.ucn.tsv"),
            os.path.join(out_dir, "chosen.tetraploid.seg.ucn"),
        )
        logging.info(f"chosen tetraploid n={n_tet}: {os.path.join(out_dir, 'chosen.tetraploid.bbc.ucn')}")
        logging.info(f"chosen tetraploid plots: {os.path.join(plot_dir, f'tetraploid_n{n_tet}')}")

    if best_type != None:
        shutil.copy2(
            os.path.join(out_dir, f"chosen.{best_type}.bbc.ucn"),
            os.path.join(out_dir, "best.bbc.ucn"),
        )
        shutil.copy2(
            os.path.join(out_dir, f"chosen.{best_type}.seg.ucn"),
            os.path.join(out_dir, "best.seg.ucn"),
        )
        best_n = {"diploid": n_dip, "tetraploid": n_tet}[best_type]
        logging.info(f"model-selected result ({best_type}): {os.path.join(out_dir, 'best.bbc.ucn')}")
        logging.info(f"best plots: {os.path.join(plot_dir, f'{best_type}_n{best_n}')}")

    return


def solve_wrapper(
    n: int,
    clonal: dict,
    gammas: dict,
    args: dict,
    ploidy: str,
    out_dir: str,
    bbcs: pd.DataFrame,
    rdr: pd.DataFrame,
    baf: pd.DataFrame,
    weights: pd.Series,
    cluster_ids: list,
    sample_ids: list,
    purities: dict,
    solve_mode="ilp",
    verbosity=0,
):
    """
    execute optimization prog
    return:
    obj
    """
    sol_dir = os.path.join(out_dir, f"sols/{ploidy}_n{n}")
    instances_dir = os.path.join(sol_dir, "instances")
    os.makedirs(sol_dir, exist_ok=True)
    os.makedirs(instances_dir, exist_ok=True)

    out_bbc = os.path.join(out_dir, f"results.{ploidy}.n{n}.bbc.ucn.tsv")
    out_seg = os.path.join(out_dir, f"results.{ploidy}.n{n}.seg.ucn.tsv")

    gammas_ = pd.Series(gammas).sort_index()
    fcn = rdr * gammas_
    f_b = fcn * baf
    f_a = fcn - f_b

    # check all user-defined fixed clonal states & fixed clone proportions. TODO
    copy_number_fixed = None

    store_solve_input(
        os.path.join(sol_dir, "input.tsv"),
        baf,
        rdr,
        fcn,
        f_a,
        f_b,
        weights,
    )

    # pre-process some args
    cn_max = {"diploid": args["diploidcmax"], "tetraploid": args["tetraploidcmax"]}[
        ploidy
    ]
    ampdel = not args["no_ampdel"]
    if args["purities"] is not None:
        purities = args["purities"]
        logging.info(f"purities overridden by user: {purities}")
    elif purities:
        logging.info(f"purities: {purities}")

    instances = solve(
        f_a=f_a,
        f_b=f_b,
        n=n,
        minprop=args["min_prop"],
        max_ncns_seg=args["num_cnstates"],
        cn_max=cn_max,
        weights=weights,
        ampdel=ampdel,
        clonal=clonal,
        purities=purities,
        baf=baf,
        copy_numbers_fixed=copy_number_fixed,
        reg_term=args["reg_term"],
        reg_steps=args["reg_steps"],
        reg_stepsize=args["reg_stepsize"],
        solver_type=args["solver"],
        solve_mode=solve_mode,
        max_iters=args["cd_niters"],
        max_convergence_iters=args["cd_convergence_iters"],
        n_seed=args["cd_nseeds"],
        n_worker=args["cd_njobs"],
        random_seed=args["cd_seed"],
        timelimit=args["timelimit"],
        instances_dir=instances_dir,
        verbose=verbosity >= 1,
        pool_size=args.get("pool_size", 1),
        pool_gap=args.get("pool_gap", None),
    )

    best_instance, imf_obj = model_selection_instance(
        f_a,
        f_b,
        weights,
        instances,
        args["reg_term"],
        solve_mode,
        sol_dir,
        verbose=verbosity >= 1,
    )

    assert best_instance != None, f"no solution for {ploidy} and n={n}"
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
    return obj, imf_obj


if __name__ == "__main__":
    args = parse_arguments_compute_cn()
    setup_logging(args)
    run(args)
