import os
import shlex
import argparse

from hatchet.hatchet_parser import (
    add_arguments_cluster_bins,
    add_arguments_compute_cn,
)


##################################################
def render_cli_args(add_args_fn, params):
    """Turn a section of config params into a CLI arg string for a subcommand.

    Values are introspected against the subcommand's parser so flags
    (store_true / BooleanOptionalAction) render correctly and any key that is
    not a real option for the subcommand (e.g. the Snakefile-only ``k``) is
    ignored, keeping the Snakefile in sync with the package's argument list.
    """
    parser = argparse.ArgumentParser(add_help=False)
    add_args_fn(parser)
    parts = []
    for action in parser._actions:
        dest = action.dest
        if dest not in params:
            continue
        val = params[dest]
        if val is None:
            continue
        long_opts = [o for o in action.option_strings if o.startswith("--")]
        opt = (
            long_opts[0]
            if long_opts
            else (action.option_strings[0] if action.option_strings else None)
        )
        if opt is None:
            continue  # positional; handled explicitly
        if isinstance(action, argparse.BooleanOptionalAction):
            parts.append(f"--{dest}" if bool(val) else f"--no-{dest}")
        elif action.nargs == 0:  # store_true / store_false flags
            if bool(val):
                parts.append(opt)
        else:
            parts.extend([opt, str(val)])
    return " ".join(shlex.quote(str(x)) for x in parts)


##################################################
# parse inputs
bb_dir = config["bb_dir"]
bbc_dir = config["bbc_dir"]
result_dir = config["result_dir"]

cluster_bins_params = config.get("cluster-bins", {})
compute_cn_params = config.get("compute-cn", {})

manual_k = compute_cn_params.get("k")

cluster_bins_args = render_cli_args(add_arguments_cluster_bins, cluster_bins_params)
compute_cn_args = render_cli_args(add_arguments_compute_cn, compute_cn_params)

##################################################
final_targets = [
    os.path.join(result_dir, "best.bbc.ucn"),
    os.path.join(result_dir, "best.seg.ucn"),
]


# entry-point
rule all:
    input:
        final_targets,
    default_target: True


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
        region_bed=config["region_bed"],
    output:
        bbc_dir=directory(bbc_dir),
        bbc=os.path.join(bbc_dir, "bulk.bbc"),
        seg=os.path.join(bbc_dir, "bulk.seg"),
    threads: config["threads"]
    params:
        args=cluster_bins_args,
    log:
        os.path.join(config["log_dir"], "cluster_bins.log"),
    benchmark:
        os.path.join(config["log_dir"], "cluster_bins.benchmark.tsv")
    shell:
        r"""
        hatchet cluster-bins \
            --bb_dir {input.bb_dir} \
            --bbc_dir {output.bbc_dir} \
            --genome_size {input.genome_size} \
            --region_bed {input.region_bed} \
            {params.args} > {log} 2>&1
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
        patient_id=config.get("patient_id", "panel"),
        args=compute_cn_args,
    log:
        os.path.join(config["log_dir"], "compute_cn.log"),
    benchmark:
        os.path.join(config["log_dir"], "compute_cn.benchmark.tsv")
    shell:
        r"""
        hatchet compute-cn \
            --bbc {params.bbc} \
            --seg {params.seg} \
            --result_dir {output.result_dir} \
            --genome_size {input.genome_size} \
            --region_bed {input.region_bed} \
            --patient_id {params.patient_id} \
            {params.args} > {log} 2>&1
        """
