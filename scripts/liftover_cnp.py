import os
import sys
import subprocess

import pandas as pd

from scripts_utils import read_seg_ucn_file

"""
Liftover HATCHet3 copy-number profile seg.ucn/bbc.ucn to different reference version
"""


def liftover_segs(
    sample: str,
    seg_ucn_file: str,
    chain_file: str,
    in_refvers: str,
    out_refvers: str,
    out_dir: str,
    max_gap=int(5e4),
):
    """
    Liftover HATCHet copy-number profiles to <out_refvers> reference version.
    UCSC chains file required
    """
    print(f"sample={sample}")
    print(f"liftover {seg_ucn} via chain file={chain_file}")
    print(f"refvers: {in_refvers}->{out_refvers}")
    print(f"max_gap: {max_gap}")
    print(f"out_dir={out_dir}")

    os.makedirs(out_dir, exist_ok=True)
    bed_in = os.path.join(out_dir, f"{sample}.{in_refvers}.bed")
    bed_out = os.path.join(out_dir, f"{sample}.{out_refvers}.bed")

    segs, clones, clone_props = read_seg_ucn_file(seg_ucn_file)
    print(f"#segments={len(segs)}")
    segs[["#CHR", "START", "END", "CNP", "PROPS"]].to_csv(
        bed_in, header=False, index=False, sep="\t"
    )

    crossmap_cmd = ["CrossMap", "bed", chain_file, bed_in, bed_out]
    try:
        subprocess.run(crossmap_cmd, check=True)
    except FileNotFoundError:
        raise RuntimeError(
            "CrossMap not found in PATH. Install with `pip install CrossMap` "
        )
    assert os.path.exists(bed_out), f"Liftover output BED not found: {bed_out}"

    lifted = pd.read_csv(
        bed_out,
        sep="\t",
        header=None,
        names=["#CHR", "START", "END", "CNP", "PROPS"],
    )
    print(f"#lifted segments={len(lifted)}")
    print(f"merge lifted segments with max-gap={max_gap}bp")
    gap = lifted["START"] - lifted["END"].shift()
    new_block = (
        (lifted["#CHR"] != lifted["#CHR"].shift())
        | (lifted["CNP"] != lifted["CNP"].shift())
        | (gap > max_gap)
    )

    block_id = new_block.cumsum()
    merged = lifted.groupby(block_id, as_index=False).agg(
        **{
            "#CHR": ("#CHR", "first"),
            "START": ("START", "first"),
            "END": ("END", "last"),
            "CNP": ("CNP", "first"),
            "PROPS": ("PROPS", "first"),
        }
    )
    print(f"#merged-segments={len(merged)}")

    out_file = os.path.join(out_dir, f"{sample}.{out_refvers}.seg.ucn.tsv")
    merged.to_csv(out_file, header=True, index=False, sep="\t")
    print(f"final output={out_file}")
    return


if __name__ == "__main__":
    _, sample, seg_ucn, chain_file, in_refvers, out_refvers, out_dir = sys.argv
    liftover_segs(sample, seg_ucn, chain_file, in_refvers, out_refvers, out_dir)
