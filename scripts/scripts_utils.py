import gzip
import subprocess
from io import StringIO
from collections import OrderedDict

import pandas as pd
import numpy as np
import pysam


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


def read_genome_sizes(sz_file: str):
    chr_sizes = OrderedDict()
    with open(sz_file, "r") as rfd:
        for line in rfd.readlines():
            ch, sizes = line.strip().split()
            chr_sizes[ch] = int(sizes)
        rfd.close()
    return chr_sizes


def read_seg_ucn_file(seg_ucn_file: str):
    segs_df = pd.read_table(seg_ucn_file, sep="\t")
    segs_df = sort_df_chr(segs_df, pos="START")

    n_clones = len([cname for cname in segs_df.columns if cname.startswith("cn_")])
    clones = ["normal"] + [f"clone{c}" for c in range(1, n_clones)]
    clone_props = segs_df[[f"u_{clone}" for clone in clones]].iloc[0].tolist()
    segs_df.loc[:, "CNP"] = segs_df.apply(
        func=lambda r: ";".join(r[f"cn_{c}"] for c in clones), axis=1
    )
    segs_df["PROPS"] = ";".join([str(p) for p in clone_props])
    return segs_df, clones, clone_props
