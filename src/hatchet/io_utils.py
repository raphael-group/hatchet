"""File readers and table reshaping for HATCHet inputs/outputs.

Readers for sample sheets, genome sizes, region BEDs, BBC and seg.ucn tables and
gamma files, plus :func:`override_solution` which loads an alternative
copy-number solution onto BBC/seg tables. Chromosome ordering and segment
merging live in :mod:`hatchet.utils`.
"""

import logging
from collections import OrderedDict

import pandas as pd

from hatchet.utils import build_seg_from_bbc, sort_df_chr


# =============================================================================
# File readers
# =============================================================================


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


def read_region_bed(bed_file: str, names=["#CHR", "START", "END", "NAME"]):
    regions = pd.read_table(
        bed_file,
        sep="\t",
        header=None,
        names=names,
    )
    return regions


def read_bbc_file(bbc_file: str):
    df = pd.read_table(bbc_file, sep="\t")
    df = sort_df_chr(df, pos="START")
    return df


def read_seg_ucn_file(seg_ucn_file: str):
    """Read a seg.ucn table, sort by chromosome, and derive the clone list.

    Per-clone proportions are per sample (``u_<clone>`` columns), so read them at
    the call site from the relevant sample's rows rather than here.

    Args:
        seg_ucn_file: Path to a seg.ucn table with ``cn_<clone>`` / ``u_<clone>``
            columns.

    Returns:
        (df, clones): the chromosome-sorted table and the clone-name list
        ("normal", "clone1", ...).
    """
    segs_df = pd.read_table(seg_ucn_file, sep="\t")
    segs_df = sort_df_chr(segs_df, pos="START")
    n_clones = len([c for c in segs_df.columns if c.startswith("cn_")])
    clones = ["normal"] + [f"clone{c}" for c in range(1, n_clones)]
    return segs_df, clones


def read_gamma_file(gamma_file: str, is_diploid=True):
    """Read per-sample RDR scaling factors from a gammas.tsv file.

    Each line is ``sample\\tgamma_diploid\\tgamma_tetraploid``.

    Args:
        gamma_file: Path to the gammas TSV.
        is_diploid: Return the diploid gamma when True, else the tetraploid one.

    Returns:
        {sample: gamma} mapping.
    """
    gammas = {}
    with open(gamma_file, "r") as fd:
        for line in fd.readlines():
            sample, gamma_dip, gamma_tet = line.strip().split("\t")
            gammas[sample] = float(gamma_dip) if is_diploid else float(gamma_tet)
    return gammas


# =============================================================================
# Solution loading
# =============================================================================


def override_solution(
    bbcs: pd.DataFrame,
    samples: list,
    clusters: list,
    n_clones: int,
    solfile: str,
    regions: pd.DataFrame,
):
    """Overwrite BBC copy-number fields with an alternative solution TSV.

    Merges a per-cluster solution onto the bins and re-segments via
    :func:`hatchet.utils.build_seg_from_bbc`.

    Returns:
        (bbcs, segs, n_clones, n_tumors, solID).
    """
    solID = solfile[str.rindex(solfile, "/") + 1 : -len(".tsv")]
    logging.info(f"overwrite BBC fields with solution {solID}!")
    sol = pd.read_table(solfile)
    assert sorted(sol.CLUSTER.unique().tolist()) == clusters
    assert sorted(sol.SAMPLE.unique().tolist()) == samples

    clones = ["normal"] + [f"clone{i}" for i in range(1, n_clones)]
    for clone in clones:
        bbcs.drop(columns=[f"u_{clone}", f"cn_{clone}"], inplace=True)

    bbcs = pd.merge(
        left=bbcs,
        right=sol,
        on=["SAMPLE", "CLUSTER"],
        how="left",
        validate="m:1",
        sort=False,
    )

    n_tumors = len([c for c in bbcs.columns.tolist() if str.startswith(c, "cn_clone")])
    n_clones = n_tumors + 1
    segs = build_seg_from_bbc(bbcs, regions)
    return bbcs, segs, n_clones, n_tumors, solID
