import os
import logging
import argparse

import numpy as np
import pandas as pd

from hatchet.utils import read_seg_ucn_file
from hatchet.evaluate.vaf_utils import estimate_vaf, is_explained_mut, relative_error


def read_snv_vcf(vcf_path, sample="tumor"):
    """Read somatic SNVs from a VCF file via cyvcf2."""
    try:
        from cyvcf2 import VCF
    except ImportError:
        raise ImportError("cyvcf2 is required for VCF input. Use --snv_tsv for TSV input.")

    v = VCF(vcf_path)
    sidx = v.samples.index(sample) if sample in v.samples else 0

    rows = []
    for rec in v:
        gt = rec.genotypes[sidx]
        dp = rec.format("DP")
        ad = rec.format("AD")
        vaf = rec.format("VAF")

        alt_i = 0
        if gt is not None and len(gt) >= 2:
            for a in gt[:2]:
                if a is not None and a > 0:
                    alt_i = a - 1
                    break

        alt = rec.ALT[alt_i] if rec.ALT and alt_i < len(rec.ALT) else None

        ref_reads = alt_reads = None
        if ad is not None:
            ad_s = ad[sidx] if ad.ndim == 2 else ad
            if len(ad_s) >= 1:
                ref_reads = int(ad_s[0])
            j = 1 + alt_i
            if len(ad_s) > j:
                alt_reads = int(ad_s[j])

        dp_val = None
        if dp is not None:
            dp_s = dp[sidx] if dp.ndim == 2 else dp
            try:
                dp_val = int(dp_s[0]) if hasattr(dp_s, "__len__") else int(dp_s)
            except Exception:
                pass
        if dp_val is None and ref_reads is not None and alt_reads is not None:
            dp_val = ref_reads + alt_reads

        vaf_val = None
        if vaf is not None:
            vaf_s = vaf[sidx] if vaf.ndim == 2 else vaf
            try:
                vaf_val = float(vaf_s[0]) if hasattr(vaf_s, "__len__") else float(vaf_s)
            except Exception:
                pass
        if vaf_val is None and dp_val and alt_reads is not None:
            vaf_val = alt_reads / dp_val

        rows.append((rec.CHROM, rec.POS, rec.REF, alt, ref_reads, alt_reads, dp_val, vaf_val))

    return pd.DataFrame(
        rows,
        columns=["#CHR", "POS", "REF", "ALT", "ref_reads", "alt_reads", "DP", "observed_VAF"],
    )


def read_snv_tsv(tsv_path):
    """Read somatic SNVs from a TSV file.

    Expected columns: #CHR, POS, ref_reads, alt_reads.
    Optional: REF, ALT, observed_VAF.
    """
    df = pd.read_csv(tsv_path, sep="\t")
    if "observed_VAF" not in df.columns:
        df["observed_VAF"] = df["alt_reads"] / (df["ref_reads"] + df["alt_reads"])
    return df


def assign_snvs_to_segments(snv_df, segs):
    """Map each SNV to the segment containing it via interval lookup.

    Uses sorted segment intervals per chromosome and np.searchsorted.
    """
    snv_df = snv_df.copy()
    snv_df["SEG_IDX"] = -1

    for chrom, seg_grp in segs.groupby("#CHR", sort=False):
        seg_grp = seg_grp.sort_values("START").reset_index(drop=True)
        starts = seg_grp["START"].to_numpy()
        ends = seg_grp["END"].to_numpy()
        seg_indices = seg_grp.index.to_numpy()

        snv_mask = snv_df["#CHR"] == chrom
        if not snv_mask.any():
            continue

        positions = snv_df.loc[snv_mask, "POS"].to_numpy()
        # searchsorted: find which segment each position falls into
        idx = np.searchsorted(starts, positions, side="right") - 1
        valid = (idx >= 0) & (idx < len(starts)) & (positions < ends[np.clip(idx, 0, len(ends) - 1)])
        snv_df.loc[snv_mask, "SEG_IDX"] = np.where(valid, seg_indices[idx], -1)

    return snv_df[snv_df["SEG_IDX"] >= 0].reset_index(drop=True)


def evaluate_snvs(segs, clones, clone_props, snv_df, gamma=0.05):
    """Evaluate each SNV against the CN solution."""
    results = []
    clone_props = np.array(clone_props)

    for _, row in snv_df.iterrows():
        seg_idx = int(row["SEG_IDX"])
        seg = segs.iloc[seg_idx]
        cns = seg["CNP"].split(";")
        obs_vaf = row["observed_VAF"]

        best, pred_vaf, ccf, allele = estimate_vaf(obs_vaf, clones, cns, clone_props)
        if best is None:
            continue

        ref_reads = int(row["ref_reads"]) if pd.notna(row.get("ref_reads")) else 0
        alt_reads = int(row["alt_reads"]) if pd.notna(row.get("alt_reads")) else 0
        explained = is_explained_mut(ref_reads, alt_reads, pred_vaf, gamma=gamma)
        rel_err = relative_error(pred_vaf, obs_vaf) if obs_vaf > 0 else np.nan
        is_subclonal = 0 < sum(np.array(list(best))[1:] > 0) < (len(clones) - 1)

        results.append({
            "#CHR": row["#CHR"],
            "POS": row["POS"],
            "observed_VAF": round(obs_vaf, 4),
            "predicted_VAF": round(pred_vaf, 4),
            "relative_error": round(rel_err, 4) if np.isfinite(rel_err) else np.nan,
            "mutated_copies": ",".join(map(str, best)),
            "allele": allele,
            "CCF": round(ccf, 4),
            "CNP": seg["CNP"],
            "is_subclonal": is_subclonal,
            "is_explained": explained,
        })

    return pd.DataFrame(results)


def run(args=None):
    logging.info("run hatchet evaluate")
    if isinstance(args, argparse.Namespace):
        args = vars(args)

    # Read seg.ucn
    seg_file = args.get("seg")
    if seg_file is None:
        result_dir = args.get("result_dir")
        if result_dir is None:
            raise ValueError("Either --seg or --result_dir must be provided")
        seg_file = os.path.join(result_dir, "best.seg.ucn")
    segs, clones, clone_props = read_seg_ucn_file(seg_file)

    samples = segs["SAMPLE"].unique().tolist()
    gamma = args.get("gamma", 0.05)
    min_vaf = args.get("min_vaf", 0.0)
    out_dir = args.get("out_dir", ".")
    os.makedirs(out_dir, exist_ok=True)

    # Read SNVs
    snv_file = args.get("snv")
    snv_tsv = args.get("snv_tsv")
    if snv_file is not None:
        snv_df = read_snv_vcf(snv_file, sample=args.get("vcf_sample", "tumor"))
    elif snv_tsv is not None:
        snv_df = read_snv_tsv(snv_tsv)
    else:
        raise ValueError("Either --snv (VCF) or --snv_tsv (TSV) must be provided")

    if min_vaf > 0:
        snv_df = snv_df[snv_df["observed_VAF"] >= min_vaf].reset_index(drop=True)
    logging.info(f"loaded {len(snv_df)} SNVs (min_vaf={min_vaf})")

    # Evaluate per sample
    all_results = []
    summary_rows = []
    for sample in samples:
        segs_s = segs[segs["SAMPLE"] == sample].reset_index(drop=True)
        purity = round(np.sum(clone_props[1:]), 3)
        if purity <= 1e-8:
            logging.warning(f"{sample}: purity too low ({purity}), skipping")
            continue

        snv_mapped = assign_snvs_to_segments(snv_df, segs_s)
        logging.info(f"{sample}: {len(snv_mapped)}/{len(snv_df)} SNVs mapped to segments")

        result_df = evaluate_snvs(segs_s, clones, clone_props, snv_mapped, gamma=gamma)
        result_df.insert(0, "SAMPLE", sample)
        all_results.append(result_df)

        n_total = len(result_df)
        n_explained = int(result_df["is_explained"].sum())
        ratio = n_explained / n_total if n_total > 0 else 0.0
        mean_err = result_df["relative_error"].mean() if n_total > 0 else np.nan
        logging.info(
            f"{sample}: {n_explained}/{n_total} explained ({ratio:.1%}), "
            f"mean relative error={mean_err:.3f}"
        )
        summary_rows.append({
            "SAMPLE": sample,
            "n_snvs": n_total,
            "n_explained": n_explained,
            "explained_ratio": round(ratio, 4),
            "mean_relative_error": round(mean_err, 4) if np.isfinite(mean_err) else np.nan,
            "purity": purity,
        })

    if all_results:
        all_df = pd.concat(all_results, ignore_index=True)
        out_snv = os.path.join(out_dir, "somatic_snvs.tsv")
        all_df.to_csv(out_snv, sep="\t", index=False)
        logging.info(f"wrote {out_snv} ({len(all_df)} SNVs)")

    if summary_rows:
        summary_df = pd.DataFrame(summary_rows)
        out_summary = os.path.join(out_dir, "eval_summary.tsv")
        summary_df.to_csv(out_summary, sep="\t", index=False)
        logging.info(f"wrote {out_summary}")

    return


if __name__ == "__main__":
    from hatchet.hatchet_parser import add_arguments_evaluate
    from hatchet.utils import setup_logging

    parser = argparse.ArgumentParser(
        prog="HATCHet evaluate",
        description="Evaluate CN solutions against somatic SNVs",
    )
    add_arguments_evaluate(parser)
    args = parser.parse_args()
    setup_logging(args)
    run(args)
