import os
import sys
import shlex
import pandas as pd

# parse inputs
bb_dir = config["bb_dir"]
bbc_dir = config["bbc_dir"]
result_dir = config["result_dir"]
plot_dir = config["plot_dir"]

# additional output-suffix
manual_k = config["compute_cn"]["k"]
if manual_k is not None:
    result_dir = f"{result_dir}_K{manual_k}"
    plot_dir = f"{plot_dir}_K{manual_k}"

sample_file = os.path.join(bb_dir, "sample_ids.tsv")
samples_df = pd.read_table(sample_file, sep="\t")
samples = samples_df["SAMPLE"].tolist()
tumor_samples = samples_df.loc[samples_df["sample_type"] == "tumor", "SAMPLE"].tolist()

##################################################
xargs = ["--verbosity", int(config["verbosity"])]

##################################################
cluster_bins_extra = []
if bool(config["cluster_bins"].get("log_rdr", False)):
    cluster_bins_extra.append("--log_rdr")
if not bool(config["cluster_bins"].get("share_tau", True)):
    cluster_bins_extra.append("--no-share_tau")
if config["cluster_bins"].get("top_restarts") is not None:
    cluster_bins_extra.extend(["--top_restarts", int(config["cluster_bins"]["top_restarts"])])
cluster_bins_args = " ".join(shlex.quote(str(x)) for x in xargs + cluster_bins_extra)

##################################################
compute_cn_args_list = []
if config["compute_cn"]["timelimit"] is not None:
    compute_cn_args_list.extend(["--timelimit", int(config["compute_cn"]["timelimit"])])
if bool(config["compute_cn"]["diploid"]):
    compute_cn_args_list.append("--diploid")
if bool(config["compute_cn"]["tetraploid"]):
    compute_cn_args_list.append("--tetraploid")
if bool(config["compute_cn"]["no_ampdel"]):
    compute_cn_args_list.append("--no_ampdel")
if config["compute_cn"].get("purities"):
    compute_cn_args_list.extend(["--purities", str(config["compute_cn"]["purities"])])
if config["compute_cn"].get("fix_cn_dip"):
    compute_cn_args_list.extend(["--fix_cn_dip", str(config["compute_cn"]["fix_cn_dip"])])
if config["compute_cn"].get("fix_cn_tet"):
    compute_cn_args_list.extend(["--fix_cn_tet", str(config["compute_cn"]["fix_cn_tet"])])
if config["compute_cn"].get("solver_threads") is not None:
    compute_cn_args_list.extend(["--solver_threads", int(config["compute_cn"]["solver_threads"])])
if config["compute_cn"].get("tree_file"):
    compute_cn_args_list.extend(["--tree_file", str(config["compute_cn"]["tree_file"])])
if config["compute_cn"].get("eps_fit") is not None:
    compute_cn_args_list.extend(["--eps_fit", float(config["compute_cn"]["eps_fit"])])
compute_cn_args = " ".join(shlex.quote(str(x)) for x in compute_cn_args_list + xargs)

##################################################
plot_cn_args_list = []
if bool(config["plot_cn"]["transparent"]):
    plot_cn_args_list.append("--transparent")
if bool(config["plot_cn"]["keep_gap"]):
    plot_cn_args_list.append("--keep_gap")
plot_cn_args = " ".join(shlex.quote(str(x)) for x in plot_cn_args_list)

##################################################
final_targets = [
    os.path.join(result_dir, "best.bbc.ucn"),
    os.path.join(result_dir, "best.seg.ucn"),
]

summary_targets = []
img_type = str(config["plot_cn"]["img_type"])
min_clone = int(config["compute_cn"]["minClone"])
max_clone = int(config["compute_cn"]["maxClone"])
if config["compute_cn"]["diploid"]:
    for n in range(min_clone, max_clone + 1):
        seg_ucn = os.path.join(result_dir, f"results.diploid.n{n}.seg.ucn.tsv")
        bbc_ucn = os.path.join(result_dir, f"results.diploid.n{n}.bbc.ucn.tsv")
        if os.path.exists(seg_ucn) and os.path.exists(bbc_ucn):
            for tumor_sample in tumor_samples:
                summary_targets.append(
                    os.path.join(
                        config["plot_dir"], f"diploid-n{n}/{tumor_sample}.1D.{img_type}"
                    )
                )
                summary_targets.append(
                    os.path.join(
                        config["plot_dir"], f"diploid-n{n}/{tumor_sample}.2D.{img_type}"
                    )
                )

if config["compute_cn"]["tetraploid"]:
    for n in range(min_clone, max_clone + 1):
        seg_ucn = os.path.join(result_dir, f"results.tetraploid.n{n}.seg.ucn.tsv")
        bbc_ucn = os.path.join(result_dir, f"results.tetraploid.n{n}.bbc.ucn.tsv")
        if os.path.exists(seg_ucn) and os.path.exists(bbc_ucn):
            for tumor_sample in tumor_samples:
                summary_targets.append(
                    os.path.join(
                        config["plot_dir"],
                        f"tetraploid-n{n}/{tumor_sample}.1D.{img_type}",
                    )
                )
                summary_targets.append(
                    os.path.join(
                        config["plot_dir"],
                        f"tetraploid-n{n}/{tumor_sample}.2D.{img_type}",
                    )
                )


# entry-point
rule all:
    input:
        final_targets,
    default_target: True


rule summary:
    input:
        summary_targets,


##################################################
rule run_cluster_bins:
    input:
        bb_dir=bb_dir,
        bb=os.path.join(bb_dir, "bb.tsv.gz"),
        sample_file=os.path.join(bb_dir, "sample_ids.tsv"),
        rdr_mfile=os.path.join(bb_dir, "bb.rdr.npz"),
        a_mfile=os.path.join(bb_dir, "bb.Aallele.npz"),
        b_mfile=os.path.join(bb_dir, "bb.Ballele.npz"),
        t_mfile=os.path.join(bb_dir, "bb.Tallele.npz"),
        genome_size=config["genome_size"],
    output:
        bbc_dir=directory(bbc_dir),
        bbc=os.path.join(bbc_dir, "bulk.bbc"),
        seg=os.path.join(bbc_dir, "bulk.seg"),
    threads: config["threads"]
    params:
        minK=int(config["cluster_bins"]["minK"]),
        maxK=int(config["cluster_bins"]["maxK"]),
        t=float(config["cluster_bins"]["t"]),
        restarts=int(config["cluster_bins"]["restarts"]),
        niters=int(config["cluster_bins"]["niters"]),
        decode_method=str(config["cluster_bins"]["decode_method"]),
        score_method=str(config["cluster_bins"]["score_method"]),
        init_method=str(config["cluster_bins"]["init_method"]),
        min_tau=float(config["cluster_bins"]["min_tau"]),
        max_tau=float(config["cluster_bins"]["max_tau"]),
        baf_eps=float(config["cluster_bins"]["baf_eps"]),
        min_covar=float(config["cluster_bins"]["min_covar"]),
        tau_iters=int(config["cluster_bins"]["tau_iters"]),
        bal_lrt_alpha=float(config["cluster_bins"]["bal_lrt_alpha"]),
        bal_margin=float(config["cluster_bins"]["bal_margin"]),
        filter_std=float(config["cluster_bins"]["filter_std"]),
        min_nbins=int(config["cluster_bins"]["min_nbins"]),
        ub_nbins=int(config["cluster_bins"]["ub_nbins"]),
        optional_args=cluster_bins_args,
    log:
        os.path.join(config["log_dir"], "cluster_bins.log"),
    shell:
        r"""
        hatchet cluster-bins \
            --bb_dir {input.bb_dir} \
            --bbc_dir {output.bbc_dir} \
            --genome_size {input.genome_size} \
            --minK {params.minK} \
            --maxK {params.maxK} \
            -t {params.t} \
            --restarts {params.restarts} \
            --niters {params.niters} \
            --decode_method {params.decode_method} \
            --score_method {params.score_method} \
            --init_method {params.init_method} \
            --min_tau {params.min_tau} \
            --max_tau {params.max_tau} \
            --baf_eps {params.baf_eps} \
            --min_covar {params.min_covar} \
            --tau_iters {params.tau_iters} \
            --bal_lrt_alpha {params.bal_lrt_alpha} \
            --bal_margin {params.bal_margin} \
            --filter_std {params.filter_std} \
            --min_nbins {params.min_nbins} \
            --ub_nbins {params.ub_nbins} \
            {params.optional_args} > {log} 2>&1
        """


##################################################
rule run_compute_cn:
    input:
        bbc_default=os.path.join(bbc_dir, "bulk.bbc"),
        seg_default=os.path.join(bbc_dir, "bulk.seg"),
        genome_size=config["genome_size"],
        region_bed=config["region_bed"],
    output:
        result_dir=directory(result_dir),
        gamma_file=os.path.join(result_dir, "gammas.tsv"),
        best_bbc_ucn=os.path.join(result_dir, "best.bbc.ucn"),
        best_seg_ucn=os.path.join(result_dir, "best.seg.ucn"),
    threads: config["threads"]
    params:
        bbc=(
            os.path.join(
                bbc_dir,
                "bulk.bbc" if manual_k is None else f"labels/bulk{manual_k}.bbc",
            )
        ),
        seg=(
            os.path.join(
                bbc_dir,
                "bulk.seg" if manual_k is None else f"labels/bulk{manual_k}.seg",
            )
        ),
        mode=str(config["compute_cn"]["mode"]),
        model_select=str(config["compute_cn"]["model_select"]),
        solver=str(config["compute_cn"]["solver"]),
        fcn_ci_alpha=float(config["compute_cn"]["fcn_ci_alpha"]),
        min_ci_margin=float(config["compute_cn"]["min_ci_margin"]),
        obj_type=str(config["compute_cn"]["obj_type"]),
        minClone=int(config["compute_cn"]["minClone"]),
        maxClone=int(config["compute_cn"]["maxClone"]),
        reg_term=str(config["compute_cn"]["reg_term"]),
        reg_steps=int(config["compute_cn"]["reg_steps"]),
        reg_bound=float(config["compute_cn"]["reg_bound"]),
        num_cnstates=int(config["compute_cn"]["num_cnstates"]),
        diploidcmax=int(config["compute_cn"]["diploidcmax"]),
        tetraploidcmax=int(config["compute_cn"]["tetraploidcmax"]),
        min_prop=float(config["compute_cn"]["min_prop"]),
        zero_cn_thres=float(config["compute_cn"]["zero_cn_thres"]),
        cd_niters=int(config["compute_cn"]["cd_niters"]),
        cd_convergence_iters=int(config["compute_cn"]["cd_convergence_iters"]),
        cd_nseeds=int(config["compute_cn"]["cd_nseeds"]),
        cd_njobs=int(config["compute_cn"]["cd_njobs"]),
        cd_seed=int(config["compute_cn"]["cd_seed"]),
        cd_tol=float(config["compute_cn"]["cd_tol"]),
        u_init=str(config["compute_cn"]["u_init"]),
        u_dir_alpha=float(config["compute_cn"]["u_dir_alpha"]),
        optional_args=compute_cn_args,
    log:
        os.path.join(config["log_dir"], "compute_cn.log"),
    shell:
        r"""
        hatchet compute-cn \
            --bbc {params.bbc} \
            --seg {params.seg} \
            --result_dir {output.result_dir} \
            --mode {params.mode} \
            --model_select {params.model_select} \
            --solver {params.solver} \
            --genome_size {input.genome_size} \
            --region_bed {input.region_bed} \
            --fcn_ci_alpha {params.fcn_ci_alpha} \
            --min_ci_margin {params.min_ci_margin} \
            --obj_type {params.obj_type} \
            --minClone {params.minClone} \
            --maxClone {params.maxClone} \
            --reg_term {params.reg_term} \
            --reg_steps {params.reg_steps} \
            --reg_bound {params.reg_bound} \
            --num_cnstates {params.num_cnstates} \
            --diploidcmax {params.diploidcmax} \
            --tetraploidcmax {params.tetraploidcmax} \
            --min_prop {params.min_prop} \
            --zero_cn_thres {params.zero_cn_thres} \
            --cd_niters {params.cd_niters} \
            --cd_convergence_iters {params.cd_convergence_iters} \
            --cd_nseeds {params.cd_nseeds} \
            --cd_njobs {params.cd_njobs} \
            --cd_seed {params.cd_seed} \
            --cd_tol {params.cd_tol} \
            --u_init {params.u_init} \
            --u_dir_alpha {params.u_dir_alpha} \
            {params.optional_args} > {log} 2>&1
        """


##################################################
rule run_plot_cn:
    input:
        bbc_ucn=lambda wc: os.path.join(
            result_dir, f"results.{wc.ploidy}.n{wc.n}.bbc.ucn.tsv"
        ),
        seg_ucn=lambda wc: os.path.join(
            result_dir, f"results.{wc.ploidy}.n{wc.n}.seg.ucn.tsv"
        ),
        gamma_file=lambda wc: os.path.join(result_dir, "gammas.tsv"),
        genome_size=config["genome_size"],
        region_bed=config["region_bed"],
    output:
        plot_dir=directory(os.path.join(config["plot_dir"], "{ploidy}-n{n}/")),
        plot1ds=[
            os.path.join(
                config["plot_dir"], "{ploidy}-n{n}/" + f"{sample_id}.1D.{img_type}"
            )
            for sample_id in tumor_samples
        ],
        plot2ds=[
            os.path.join(
                config["plot_dir"], "{ploidy}-n{n}/" + f"{sample_id}.2D.{img_type}"
            )
            for sample_id in tumor_samples
        ],
    wildcard_constraints:
        ploidy="(diploid|tetraploid)",
    threads: config["threads"]
    params:
        dpi=int(config["plot_cn"]["dpi"]),
        img_type=str(config["plot_cn"]["img_type"]),
        tail_alpha=float(config["plot_cn"]["tail_alpha"]),
        center_alpha=float(config["plot_cn"]["center_alpha"]),
        onetail_area=float(config["plot_cn"]["onetail_area"]),
        maxlim_fcn=int(config["plot_cn"]["maxlim_fcn"]),
        optional_args=plot_cn_args,
    log:
        os.path.join(config["log_dir"], "plot_cn.{ploidy}-n{n}.log"),
    shell:
        r"""
        hatchet plot-cn \
            --bbc {input.bbc_ucn} \
            --seg {input.seg_ucn} \
            --plot_dir {output.plot_dir} \
            --gamma_file {input.gamma_file} \
            --genome_size {input.genome_size} \
            --region_bed {input.region_bed} \
            --dpi {params.dpi} \
            --img_type {params.img_type} \
            --tail_alpha {params.tail_alpha} \
            --center_alpha {params.center_alpha} \
            --onetail_area {params.onetail_area} \
            --maxlim_fcn {params.maxlim_fcn} \
            --ploidy {wildcards.ploidy} \
            {params.optional_args} > {log} 2>&1
        """
