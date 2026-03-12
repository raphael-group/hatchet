import os
import sys

import pandas as pd
import numpy as np

from scripts_utils import *

# overwrite HATCHet solution with regularized alternatives
if __name__ == "__main__":
    _, fbbc, solfile, outprefix = sys.argv
    solID = solfile[str.rindex(solfile, "/") + 1 : -len(".tsv")]
    print(f"overwrite BBC fields with solution {solID}!")

    bbcs = pd.read_table(fbbc, sep="\t")
    sol = pd.read_table(solfile, sep="\t")
    print(f"bbc clusters: ", np.unique(bbcs["CLUSTER"]))
    print(f"sol clusters: ", np.unique(sol["CLUSTER"]))
    n_clones = (
        len([c for c in sol.columns.tolist() if str.startswith(c, "cn_clone")]) + 1
    )
    clones = [f"normal"] + [f"clone{i}" for i in range(1, n_clones)]

    cnp_cols = []
    for clone in clones:
        cnp_cols.extend([f"cn_{clone}", f"u_{clone}"])
    sol = sol[["CLUSTER", "SAMPLE"] + cnp_cols].reset_index(drop=True)

    bbcs = pd.merge(
        left=bbcs, right=sol, on=["SAMPLE", "CLUSTER"], how="left", sort=False
    )
    bbcs.dropna(subset=["cn_normal"], inplace=True)
    bbcs = bbcs.sort_values(by="SAMPLE")
    bbcs = bbcs.reset_index(drop=True)

    chs = sort_chroms(bbcs["#CHR"].unique().tolist())
    bbcs["#CHR"] = pd.Categorical(bbcs["#CHR"], categories=chs, ordered=True)
    bbcs_sorted = bbcs.sort_values(["SAMPLE", "#CHR", "START"]).copy()
    bbcs_sorted.to_csv(f"{outprefix}.bbc.ucn.tsv", sep="\t", header=True, index=False)

    bbcs_sorted = bbcs_sorted[
        ["#CHR", "START", "END", "SAMPLE"] + cnp_cols
    ].reset_index(drop=True)
    cn_cols = [f"cn_{clone}" for clone in clones]
    cn_changed = (bbcs_sorted[cn_cols] != bbcs_sorted[cn_cols].shift()).any(axis=1)
    bbcs_sorted["_gid"] = (
        (bbcs_sorted["SAMPLE"] != bbcs_sorted["SAMPLE"].shift())
        | (bbcs_sorted["#CHR"] != bbcs_sorted["#CHR"].shift())
        | (bbcs_sorted["START"] != bbcs_sorted["END"].shift())
        | cn_changed
    ).cumsum()

    agg = {"SAMPLE": "first", "#CHR": "first", "START": "min", "END": "max"}
    agg.update({c: "first" for c in cnp_cols})

    segs = bbcs_sorted.groupby("_gid", as_index=False).agg(agg).drop(columns="_gid")
    segs.to_csv(f"{outprefix}.seg.ucn.tsv", sep="\t", header=True, index=False)
