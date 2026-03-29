import os
import time
import logging
import resource
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
    return prepare_seg_ucn(segs_df)


def prepare_seg_ucn(segs_df: pd.DataFrame):
    """Add CNP/PROPS columns to a seg UCN DataFrame and return (df, clones, clone_props)."""
    segs_df = sort_df_chr(segs_df, pos="START")

    n_clones = len([cname for cname in segs_df.columns if cname.startswith("cn_")])
    clones = ["normal"] + [f"clone{c}" for c in range(1, n_clones)]
    segs_df.loc[:, "CNP"] = segs_df.apply(
        func=lambda r: ";".join(r[f"cn_{c}"] for c in clones), axis=1
    )
    segs_df["PROPS"] = segs_df.apply(
        func=lambda r: ";".join(str(r[f"u_{c}"]) for c in clones), axis=1
    )

    # Clone proportions are read from the first row; all rows share the same
    # per-sample proportions for a given clone, so any row gives the same result.
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


def build_seg_from_bbc(df: pd.DataFrame, regions: pd.DataFrame) -> pd.DataFrame:
    """Build a segment-level DataFrame from a bin-level BBC DataFrame with CN columns.

    Adjacent bins with the same copy-number state are merged into segments.
    Region boundaries (from a BED file) act as merge barriers so that segments
    never span across regions.

    Parameters
    ----------
    df : pd.DataFrame
        Bin-level DataFrame that must contain columns ``#CHR``, ``START``,
        ``END``, ``SAMPLE``, plus ``cn_*`` and ``u_*`` columns.
    regions : pd.DataFrame
        Region BED DataFrame with ``#CHR``, ``START``, ``END`` columns.

    Returns
    -------
    pd.DataFrame
        Segment-level DataFrame with the same CN/u columns.
    """
    cn_cols = [c for c in df.columns if c.startswith("cn_")]
    u_cols = [c for c in df.columns if c.startswith("u_")]
    extra_columns = [col for pair in zip(cn_cols, u_cols) for col in pair]

    df = df.sort_values(["#CHR", "START", "END", "SAMPLE"]).reset_index(drop=True)
    df["all_copy_numbers"] = df[cn_cols].apply(",".join, axis=1)
    first_sample = df["SAMPLE"].iloc[0]
    df["segment"] = (
        (df["SAMPLE"] == first_sample)
        & (
            (df["#CHR"] != df["#CHR"].shift())
            | (df["all_copy_numbers"] != df["all_copy_numbers"].shift())
            | (df["START"] != df["END"].shift())
        )
    ).cumsum()

    agg = {"#CHR": "first", "START": "min", "END": "max", "SAMPLE": "first"}
    if "CLUSTER" in df.columns:
        agg["CLUSTER"] = "first"
    agg.update({c: "first" for c in extra_columns})
    seg = df.groupby(["segment", "SAMPLE"]).agg(agg).reset_index(drop=True)

    base_cols = ["#CHR", "START", "END", "SAMPLE"]
    if "CLUSTER" in seg.columns:
        base_cols.append("CLUSTER")
    out_cols = base_cols + extra_columns

    # Assign each segment to a region index (-1 = outside all regions)
    seg["_region"] = -1
    for r_idx, region in regions.iterrows():
        mask = (
            (seg["#CHR"] == region["#CHR"])
            & (seg["START"] >= region["START"])
            & (seg["END"] <= region["END"])
        )
        seg.loc[mask, "_region"] = r_idx

    # Merge adjacent same-CN segments within each region
    merged_rows = []
    for _, grp in seg.groupby(["SAMPLE", "#CHR", "_region"], sort=False, observed=True):
        grp = grp.sort_values("START").reset_index(drop=True)
        state_key = grp[cn_cols].apply(tuple, axis=1)
        grp["_run"] = (state_key != state_key.shift()).cumsum()
        for _, run_grp in grp.groupby("_run"):
            row = run_grp.iloc[0].copy()
            row["START"] = run_grp["START"].min()
            row["END"] = run_grp["END"].max()
            merged_rows.append(row[out_cols])

    out = sort_df_chr(pd.DataFrame(merged_rows, columns=out_cols), pos="START")
    out = out.sort_values(["#CHR", "START", "SAMPLE"]).reset_index(drop=True)
    return out


def compute_clone_ploidies(segs: pd.DataFrame, clones: list):
    """Length-weighted average total CN per clone.

    Returns a dict mapping clone name to ploidy (float).
    """
    segment_lengths = (segs["END"] - segs["START"]).to_numpy()
    total_length = np.sum(segment_lengths)
    ploidies = {}
    for clone in clones:
        cn = (
            segs[f"cn_{clone}"]
            .apply(func=lambda v: int(v.split("|")[0]) + int(v.split("|")[1]))
            .to_numpy()
        )
        ploidies[clone] = float(np.sum(cn * segment_lengths) / total_length)
    return ploidies


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

    total_length = np.sum(segment_lengths)
    rho = (1 / tumor_purity) * (rho / total_length)
    return rho


_NOISY_LOGGERS = ["adjustText", "fontTools", "matplotlib", "numba", "pyomo"]


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
    for name in _NOISY_LOGGERS:
        logging.getLogger(name).setLevel(logging.WARNING)


def add_file_logging(out_dir: str, command: str = "hatchet") -> None:
    """Attach a FileHandler to the root logger so logs are also written to *out_dir/<command>.log*."""
    os.makedirs(out_dir, exist_ok=True)
    level = (
        logging.root.level if logging.root.level != logging.WARNING else logging.INFO
    )
    fh = logging.FileHandler(os.path.join(out_dir, f"{command}.log"), mode="w")
    fh.setLevel(level)
    fh.setFormatter(
        logging.Formatter(
            "%(asctime)s.%(msecs)03d %(levelname)s %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )
    )
    logging.root.addHandler(fh)
    if logging.root.level > level:
        logging.root.setLevel(level)
    for name in _NOISY_LOGGERS:
        logging.getLogger(name).setLevel(logging.WARNING)


def log_arguments(args) -> None:
    d = vars(args) if hasattr(args, "__dict__") else args
    lines = "\n".join(f"  {k}: {v}" for k, v in sorted(d.items()) if k != "func")
    logging.info(f"parsed arguments:\n{lines}")
