import os
import logging
import argparse

import numpy as np
import pandas as pd

from hatchet.utils import read_seg_ucn_file, read_genome_sizes, read_region_bed
from hatchet.evaluate.evaluate_utils import (
    read_snv_vcf,
    read_snv_tsv,
    assign_snvs_to_segments,
    evaluate_snvs,
    evaluate_pool_solutions,
    plot_vaf_1d,
)


def run(args=None):
    logging.info("run hatchet evaluate")
    if isinstance(args, argparse.Namespace):
        args = vars(args)

    gamma = args.get("gamma", 0.05)
    min_vaf = args.get("min_vaf", 0.0)
    min_depth = args.get("min_depth", 5)
    out_dir = args.get("out_dir", ".")
    os.makedirs(out_dir, exist_ok=True)
    eval_all = args.get("eval_all", False)

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
    if min_depth > 0 and "DP" in snv_df.columns:
        snv_df = snv_df[snv_df["DP"] >= min_depth].reset_index(drop=True)
    logging.info(
        f"loaded {len(snv_df)} SNVs (min_vaf={min_vaf}, min_depth={min_depth})"
    )

    result_dir = args.get("result_dir")

    # Mode 1: evaluate all pool solutions
    if eval_all and result_dir is not None:
        logging.info("evaluating all pool solutions")
        pool_eval = evaluate_pool_solutions(result_dir, snv_df, gamma=gamma)
        out_path = os.path.join(out_dir, "pool_eval.tsv")
        pool_eval.to_csv(out_path, sep="\t", index=False)
        logging.info(f"wrote {out_path} ({len(pool_eval)} solutions)")
        return

    # Mode 2: evaluate best solution only
    seg_file = args.get("seg")
    if seg_file is None:
        if result_dir is None:
            raise ValueError("Either --seg or --result_dir must be provided")
        seg_file = os.path.join(result_dir, "best.seg.ucn")
    segs, clones, clone_props = read_seg_ucn_file(seg_file)

    samples = segs["SAMPLE"].unique().tolist()

    all_results = []
    summary_rows = []
    for sample in samples:
        segs_s = segs[segs["SAMPLE"] == sample].reset_index(drop=True)
        purity = round(np.sum(clone_props[1:]), 3)
        if purity <= 1e-8:
            logging.warning(f"{sample}: purity too low ({purity}), skipping")
            continue

        snv_mapped = assign_snvs_to_segments(snv_df, segs_s)
        logging.info(
            f"{sample}: {len(snv_mapped)}/{len(snv_df)} SNVs mapped to segments"
        )

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
        summary_rows.append(
            {
                "SAMPLE": sample,
                "n_snvs": n_total,
                "n_explained": n_explained,
                "explained_ratio": round(ratio, 4),
                "mean_relative_error": round(mean_err, 4)
                if np.isfinite(mean_err)
                else np.nan,
                "purity": purity,
            }
        )

    if all_results:
        all_df = pd.concat(all_results, ignore_index=True)
        out_snv = os.path.join(out_dir, "somatic_snvs.tsv")
        all_df.to_csv(out_snv, sep="\t", index=False)
        logging.info(f"wrote {out_snv} ({len(all_df)} SNVs)")

        # Plot VAF along genome
        genome_size = args.get("genome_size")
        region_bed = args.get("region_bed")
        if genome_size is not None:
            from hatchet.plot.plot_utils import get_expected_baf_fcn

            chrom_sizes = read_genome_sizes(genome_size)
            regions = read_region_bed(region_bed)
            for sample in samples:
                sample_df = all_df[all_df["SAMPLE"] == sample]
                if len(sample_df) == 0:
                    continue
                segs_s = segs[segs["SAMPLE"] == sample].reset_index(drop=True)
                cn_cols = sorted(
                    [col for col in segs_s.columns if col.startswith("cn_")],
                    key=lambda col: (0 if col == "cn_normal" else 1, col),
                )
                segs_plot = segs_s.copy()
                segs_plot["predicted_VAF"] = 0.0
                for idx, seg in segs_plot.iterrows():
                    a_b = [seg[col].split("|") for col in cn_cols]
                    states = [(int(a), int(b)) for a, b in a_b]
                    _, _, _, exp_baf = get_expected_baf_fcn(states, clone_props)
                    segs_plot.at[idx, "predicted_VAF"] = exp_baf
                out_plot = os.path.join(out_dir, f"{sample}.vaf_1d.pdf")
                plot_vaf_1d(sample_df, segs_plot, chrom_sizes, regions, out_plot)

    if summary_rows:
        summary_df = pd.DataFrame(summary_rows)
        out_summary = os.path.join(out_dir, "eval_summary.tsv")
        summary_df.to_csv(out_summary, sep="\t", index=False)
        logging.info(f"wrote {out_summary}")


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
