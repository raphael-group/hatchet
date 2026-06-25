import argparse
import os
import sys
import time
import logging
import resource
import threading
from collections import OrderedDict
from importlib.resources import files

import pandas as pd
import numpy as np
import yaml

try:
    import psutil
except ImportError:
    psutil = None


def load_defaults() -> dict:
    """Load packaged tuning defaults from src/hatchet/hatchet.yaml."""
    text = files("hatchet").joinpath("hatchet.yaml").read_text(encoding="utf-8")
    return yaml.safe_load(text)


def normalize_args(args) -> dict:
    """Convert Namespace→dict if needed and merge YAML defaults under it (args win)."""
    if isinstance(args, argparse.Namespace):
        args = vars(args)
    return {**load_defaults(), **args}


class _RSSSampler(threading.Thread):
    """Track the peak summed RSS of this process and all descendants over time."""

    def __init__(self, interval=0.2):
        super().__init__(daemon=True)
        self.interval = interval
        self._stop_event = threading.Event()
        self.peak_bytes = 0
        self._proc = psutil.Process()

    def run(self):
        while not self._stop_event.is_set():
            try:
                total = self._proc.memory_info().rss
                for c in self._proc.children(recursive=True):
                    try:
                        total += c.memory_info().rss
                    except (psutil.NoSuchProcess, psutil.AccessDenied):
                        pass
                self.peak_bytes = max(self.peak_bytes, total)
            except (psutil.NoSuchProcess, psutil.AccessDenied):
                pass
            self._stop_event.wait(self.interval)

    def stop(self):
        self._stop_event.set()
        self.join(timeout=2 * self.interval)


def log_step_start():
    """Start a profiling step; returns finish(name, out_file=None) that logs wall/cpu/peak_rss.

    peak_rss is the concurrent peak of the whole process tree (workers, CBC
    subprocess) sampled via psutil, falling back to rusage max(self, child).
    """
    t_wall = time.perf_counter()
    t_cpu = time.process_time()
    # ru_maxrss: bytes on macOS, kibibytes on Linux
    _rss_to_gb = (1 / 1e9) if sys.platform == "darwin" else (1024 / 1e9)

    sampler = None
    if psutil is not None:
        sampler = _RSSSampler()
        sampler.start()

    def finish(name, out_file=None):
        wall = time.perf_counter() - t_wall
        cpu = time.process_time() - t_cpu
        ru_self = resource.getrusage(resource.RUSAGE_SELF)
        ru_child = resource.getrusage(resource.RUSAGE_CHILDREN)
        peak_self = ru_self.ru_maxrss * _rss_to_gb
        peak_child = ru_child.ru_maxrss * _rss_to_gb  # largest single child
        child_cpu = ru_child.ru_utime + ru_child.ru_stime
        if sampler is not None:
            sampler.stop()
            peak = sampler.peak_bytes / 1e9  # psutil RSS is bytes
            peak = max(
                peak, peak_self
            )  # sampler may miss instantaneous high-water mark
        else:
            peak = max(peak_self, peak_child)
        rows = [
            ("step", name),
            ("wall_s", f"{wall:.3f}"),
            ("cpu_self_s", f"{cpu:.3f}"),
            ("cpu_children_s", f"{child_cpu:.3f}"),
            ("peak_rss_tree_gb", f"{peak:.3f}"),
            ("peak_rss_self_gb", f"{peak_self:.3f}"),
            ("peak_rss_largest_child_gb", f"{peak_child:.3f}"),
        ]
        logging.info(f"{name} runtime:\n" + "\n".join(f"  {k}: {v}" for k, v in rows))
        if out_file is not None:
            with open(out_file, "w") as fh:
                fh.write("".join(f"{k}\t{v}\n" for k, v in rows))

    return finish


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
    sample_types = sample_df["sample_type"].tolist()
    if "normal" in sample_types:
        normal_idx = [i for i, t in enumerate(sample_types) if t == "normal"]
        tumor_idx = [i for i, t in enumerate(sample_types) if t == "tumor"]
    else:
        normal_idx = []
        tumor_idx = list(range(len(sample_types)))

    assays = (
        sample_df["assay_type"].tolist()
        if "assay_type" in sample_df.columns
        else [None] * len(sample_df)
    )
    normal_set, tumor_set = set(normal_idx), set(tumor_idx)
    assay2samples = {}
    for i, a in enumerate(assays):
        grp = assay2samples.setdefault(a, {"normal": [], "tumor": []})
        if i in normal_set:
            grp["normal"].append(i)
        if i in tumor_set:
            grp["tumor"].append(i)
    return sample_df, normal_idx, tumor_idx, assay2samples


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
