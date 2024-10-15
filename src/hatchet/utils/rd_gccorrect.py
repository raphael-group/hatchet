import pandas as pd
import statsmodels.formula.api as smf
from pybedtools import BedTool


def rd_gccorrect(bb, ref_genome):
    """
    Function to correct GC bias in read depth data for each sample.
    Parameters:
    - bb: DataFrame containing read depth data, including columns '#CHR', 'START', 'END', 'RD', 'SAMPLE'
    - ref_genome: File path to the reference genome in FASTA format
    Return type: DataFrame with corrected read depth data
    """

    bb["CHR"] = bb["#CHR"].str.replace("chr", "")
    bb = bb.sort_values(by=["CHR", "START"]).reset_index(drop=True)

    # add GC content
    # ref_genome = '/n/fs/ragr-data/datasets/ref-genomes/GRCh38_10X/fasta/genome.fa'
    bb = bb.merge(
        BedTool.from_dataframe(bb[["#CHR", "START", "END"]].drop_duplicates())
        .nucleotide_content(fi=ref_genome)
        .to_dataframe(disable_auto_names=True)
        .rename(
            columns={
                "#1_usercol": "#CHR",
                "2_usercol": "START",
                "3_usercol": "END",
                "5_pct_gc": "GC",
            }
        )[["#CHR", "START", "END", "GC"]]
    )

    # Correcting GC bias per sample
    gccorrect = bb.groupby("SAMPLE").apply(
        lambda D: smf.quantreg("RD ~ GC + I(GC ** 2.0)", data=D).fit(q=0.5)
    )
    bb["UNCORR_RD"] = bb["RD"].copy()
    bb["GCCORR"] = bb.groupby("SAMPLE")["GC"].transform(
        lambda X: gccorrect[X.name].predict(X.to_frame("GC")).values
    )
    bb["RD"] = bb["UNCORR_RD"] / bb["GCCORR"].where(
        (bb["GCCORR"] > 0) & ~pd.isnull(bb["GCCORR"]), 1
    )
    bb["RD"] = bb["RD"] / bb.groupby("SAMPLE")["RD"].transform(lambda X: X.mean())
    print(bb.head(10))

    bb = bb.sort_values(by=["CHR", "START"]).reset_index(drop=True)

    # Drop the 'CHR' column
    bb.drop(columns=["CHR"], inplace=True)

    return bb
