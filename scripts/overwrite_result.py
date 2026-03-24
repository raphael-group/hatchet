import os
import sys

import pandas as pd
import numpy as np

from scripts_utils import sort_chroms

# Overwrite HATCHet BBC with a pool solution file, produce bbc.ucn.tsv + seg.ucn.tsv
if __name__ == "__main__":
    _, fbbc, solfile, outprefix = sys.argv
    solID = os.path.splitext(os.path.basename(solfile))[0]
    print(f"overwrite BBC with solution {solID}")

    bbcs = pd.read_table(fbbc, sep="\t")
    sol = pd.read_table(solfile, sep="\t")

    # Detect clone columns from solution file
    cn_col_names = [c for c in sol.columns if c.startswith("cn_")]
    u_col_names = [c for c in sol.columns if c.startswith("u_")]
    cnp_cols = [c for pair in zip(cn_col_names, u_col_names) for c in pair]

    # Drop any pre-existing CN/u columns from bbcs to avoid merge conflicts
    bbcs = bbcs.drop(columns=[c for c in cnp_cols if c in bbcs.columns], errors="ignore")

    sol = sol[["CLUSTER", "SAMPLE"] + cnp_cols].drop_duplicates()
    bbcs = bbcs.merge(sol, on=["CLUSTER", "SAMPLE"], how="left")
    bbcs = bbcs.dropna(subset=[cn_col_names[0]]).reset_index(drop=True)

    chs = sort_chroms(bbcs["#CHR"].unique().tolist())
    bbcs["#CHR"] = pd.Categorical(bbcs["#CHR"], categories=chs, ordered=True)
    bbcs = bbcs.sort_values(["SAMPLE", "#CHR", "START"]).reset_index(drop=True)
    bbcs.to_csv(f"{outprefix}.bbc.ucn.tsv", sep="\t", header=True, index=False)

    # Build seg: merge adjacent bins with same sample/chr/CN into segments
    sub = bbcs[["#CHR", "START", "END", "SAMPLE"] + cnp_cols].copy()
    sub["_gid"] = (
        (sub["SAMPLE"] != sub["SAMPLE"].shift())
        | (sub["#CHR"] != sub["#CHR"].shift())
        | (sub["START"] != sub["END"].shift())
        | (sub[cn_col_names] != sub[cn_col_names].shift()).any(axis=1)
    ).cumsum()

    agg = {"SAMPLE": "first", "#CHR": "first", "START": "min", "END": "max"}
    agg.update({c: "first" for c in cnp_cols})
    segs = sub.groupby("_gid", as_index=False).agg(agg).drop(columns="_gid")
    segs.to_csv(f"{outprefix}.seg.ucn.tsv", sep="\t", header=True, index=False)

    print(f"wrote {outprefix}.bbc.ucn.tsv and {outprefix}.seg.ucn.tsv")
