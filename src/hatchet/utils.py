import os
import sys
import gzip
import time
import logging
import resource
import subprocess
from io import StringIO
from collections import OrderedDict
from contextlib import contextmanager

import pandas as pd
import numpy as np


@contextmanager
def log_step(name):
    """Context manager that logs wall time, CPU time, and peak RSS on exit.

    Usage::

        with log_step("cluster-bins"):
            ...  # work
        # logs: cluster-bins done: wall=12.3s, cpu=45.6s, peak_rss=1.23 GB
    """
    t_wall = time.perf_counter()
    t_cpu = time.process_time()
    yield
    wall = time.perf_counter() - t_wall
    cpu = time.process_time() - t_cpu
    # ru_maxrss is in bytes on macOS, kilobytes on Linux
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    rss_gb = rss / 1e9 if rss > 1e9 else rss / 1e6
    logging.info(
        f"{name} done: wall={wall:.1f}s, cpu={cpu:.1f}s, peak_rss={rss_gb:.2f} GB"
    )


def log_step_start():
    """Start a profiling step. Returns a callable that logs the summary.

    Usage::

        log_done = log_step_start()
        ...  # work
        log_done("cluster-bins")
        # logs: cluster-bins done: wall=12.3s, cpu=45.6s, peak_rss=1.23 GB
    """
    t_wall = time.perf_counter()
    t_cpu = time.process_time()

    def finish(name):
        wall = time.perf_counter() - t_wall
        cpu = time.process_time() - t_cpu
        rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        rss_gb = rss / 1e9 if rss > 1e9 else rss / 1e6
        logging.info(
            f"{name} done: wall={wall:.1f}s, cpu={cpu:.1f}s, peak_rss={rss_gb:.2f} GB"
        )

    return finish


def symlink_force(src, dst):
    try:
        os.remove(dst)
    except FileNotFoundError:
        pass
    os.symlink(os.path.abspath(src), os.path.abspath(dst))


def get_ord2chr(ch="chr"):
    return [f"{ch}{i}" for i in range(1, 23)] + [f"{ch}X", f"{ch}Y"]


def get_chr2ord(ch):
    chr2ord = {}
    for i in range(1, 23):
        chr2ord[f"{ch}{i}"] = i
    chr2ord[f"{ch}X"] = 23
    chr2ord[f"{ch}Y"] = 24
    return chr2ord


def sort_chroms(chromosomes: list):
    assert len(chromosomes) != 0
    chromosomes = [str(c) for c in chromosomes]
    ch = "chr" if str(chromosomes[0]).startswith("chr") else ""
    chr2ord = get_chr2ord(ch)
    return sorted(chromosomes, key=lambda x: chr2ord[x])


def sort_df_chr(df: pd.DataFrame, ch="#CHR", pos="POS"):
    chs = sort_chroms(df[ch].unique().tolist())
    df[ch] = pd.Categorical(df[ch], categories=chs, ordered=True)
    df.sort_values(by=[ch, pos], inplace=True, ignore_index=True)
    return df


def read_sample_file(sample_file: str):
    sample_df = pd.read_table(sample_file, sep="\t")
    samples = sample_df["SAMPLE"].tolist()
    sample_types = sample_df["sample_type"].tolist()
    no_normal = "normal" not in sample_types
    return sample_df, samples, no_normal


def read_genome_sizes(sz_file: str):
    chr_sizes = OrderedDict()
    with open(sz_file, "r") as rfd:
        for line in rfd.readlines():
            ch, sizes = line.strip().split()
            chr_sizes[ch] = int(sizes)
        rfd.close()
    return chr_sizes


def read_bbc_file(bbc_file: str):
    df = pd.read_table(bbc_file, sep="\t")
    df = sort_df_chr(df, pos="START")
    return df


def read_seg_ucn_file(seg_ucn_file: str):
    segs_df = pd.read_table(seg_ucn_file, sep="\t")
    segs_df = sort_df_chr(segs_df, pos="START")

    n_clones = len([cname for cname in segs_df.columns if cname.startswith("cn_")])
    clones = ["normal"] + [f"clone{c}" for c in range(1, n_clones)]
    segs_df.loc[:, "CNP"] = segs_df.apply(
        func=lambda r: ";".join(r[f"cn_{c}"] for c in clones), axis=1
    )
    segs_df["PROPS"] = segs_df.apply(
        func=lambda r: ";".join(str(r[f"u_{c}"]) for c in clones), axis=1
    )

    # TODO fix 1clone
    clone_props = segs_df[[f"u_{clone}" for clone in clones]].iloc[0].tolist()
    return segs_df, clones, clone_props


def read_region_bed(bed_file: str, names=["#CHR", "START", "END", "NAME"]):
    regions = pd.read_table(
        bed_file,
        sep="\t",
        header=None,
        names=names,
    )
    return regions


def compute_tumor_ploidy(segs: pd.DataFrame, clones: list, tumor_purity: float):
    if tumor_purity <= 1e-8:
        return 0.0
    segment_lengths = (segs["END"] - segs["START"]).to_numpy()
    rho = 0.0
    for clone in clones[1:]:
        u = segs[f"u_{clone}"].iloc[0]
        cn = (
            segs[f"cn_{clone}"]
            .apply(func=lambda v: int(v.split("|")[0]) + int(v.split("|")[1]))
            .to_numpy()
        )
        rho += u * np.sum(cn * segment_lengths)

    l = np.sum(segment_lengths)
    rho = (1 / tumor_purity) * (rho / l)
    return rho


def setup_logging(args) -> None:
    d = vars(args) if hasattr(args, "__dict__") else args
    verbosity = d.get("verbosity", 1)
    level = {0: logging.WARNING, 1: logging.INFO, 2: logging.DEBUG}.get(
        verbosity, logging.DEBUG
    )
    logging.basicConfig(
        level=level,
        format="%(asctime)s.%(msecs)03d %(levelname)s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        force=True,
    )
    logging.getLogger("adjustText").setLevel(logging.ERROR)
    logging.getLogger("fontTools").setLevel(logging.ERROR)
    logging.getLogger("jax").setLevel(logging.ERROR)
    logging.getLogger("matplotlib").setLevel(logging.ERROR)
    logging.getLogger("numba").setLevel(logging.ERROR)
    logging.getLogger("pyomo").setLevel(logging.WARNING)


def log_arguments(args) -> None:
    d = vars(args) if hasattr(args, "__dict__") else args
    lines = "\n".join(f"  {k}: {v}" for k, v in sorted(d.items()) if k != "func")
    logging.info(f"parsed arguments:\n{lines}")
