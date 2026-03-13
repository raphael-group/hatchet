import os
import sys
import logging
from collections import OrderedDict
import numpy as np
import pandas as pd

import seaborn as sns

WHITE = (1, 1, 1, 1)
BLACK = (0, 0, 0, 1)
RED = (1, 0, 0, 1)
BLUE = (0, 0, 1, 1)


def load_gammas(gamma_file: str, is_diploid=True):
    gammas = {}
    with open(gamma_file, "r") as fd:
        for line in fd.readlines():
            sample, gamma_dip, gamma_tet = line.strip().split("\t")
            if is_diploid:
                gammas[sample] = float(gamma_dip)
            else:
                gammas[sample] = float(gamma_tet)
        fd.close()
    return gammas


def override_solution(
    bbcs: pd.DataFrame,
    samples: list,
    clusters: list,
    n_clones: int,
    solfile: str,
    regions: pd.DataFrame,
):
    from hatchet.utils import build_seg_from_bbc

    solID = solfile[str.rindex(solfile, "/") + 1 : -len(".tsv")]
    logging.info(f"overwrite BBC fields with solution {solID}!")
    sol = pd.read_table(solfile)
    assert sorted(sol.CLUSTER.unique().tolist()) == clusters
    assert sorted(sol.SAMPLE.unique().tolist()) == samples

    clones = ["normal"] + [f"clone{i}" for i in range(1, n_clones)]
    for clone in clones:
        bbcs.drop(columns=[f"u_{clone}", f"cn_{clone}"], inplace=True)
    bbcs.drop(columns=["CNP", "PROPS"], inplace=True, errors="ignore")

    bbcs = pd.merge(
        left=bbcs,
        right=sol,
        on=["SAMPLE", "CLUSTER"],
        how="left",
        validate="m:1",
        sort=False,
    )

    # Rebuild derived columns
    n_tumors = len([c for c in bbcs.columns.tolist() if str.startswith(c, "cn_clone")])
    n_clones = n_tumors + 1
    clones = ["normal"] + [f"clone{i}" for i in range(1, n_clones)]
    bbcs["CNP"] = bbcs.apply(
        func=lambda r: ";".join(r[f"cn_{c}"] for c in clones), axis=1
    )
    bbcs["PROPS"] = bbcs.apply(
        func=lambda r: ";".join(str(r[f"u_{c}"]) for c in clones), axis=1
    )

    segs = build_seg_from_bbc(bbcs, regions)
    segs["CNP"] = segs.apply(
        func=lambda r: ";".join(r[f"cn_{c}"] for c in clones), axis=1
    )
    segs["PROPS"] = segs.apply(
        func=lambda r: ";".join(str(r[f"u_{c}"]) for c in clones), axis=1
    )

    return bbcs, segs, n_clones, n_tumors, solID


def cn2total(s):
    tkns = s.split("|")
    assert len(tkns) == 2
    return int(tkns[0]) + int(tkns[1])


def get_expected_baf_fcn(cns, props):
    assert len(cns) == len(props)
    A = np.array([x[0] for x in cns])
    B = np.array([x[1] for x in cns])
    y_fcn_a = np.sum(A * props)
    y_fcn_b = np.sum(B * props)
    y_fcn = y_fcn_a + y_fcn_b
    y_baf = np.sum(B * props) / y_fcn

    return y_fcn_a, y_fcn_b, y_fcn, y_baf


def set_palette(num_colors=8, style="whitegrid"):
    sns.set_style(style)
    if num_colors > 8:
        palette = sns.color_palette("husl", n_colors=num_colors)
    else:
        palette = sns.color_palette("Set2", n_colors=num_colors)
    sns.set_palette(palette)
    return palette
