#!/usr/bin/env python
"""TOST equivalence test for whether each cluster's BAF = 0.5."""

import argparse
import sys
from pathlib import Path

import pandas as pd
from scipy import stats


def main():
    parser = argparse.ArgumentParser(
        description="Test cluster BAF = 0.5 via TOST equivalence test"
    )
    parser.add_argument(
        "--seg",
        default=str(Path(__file__).resolve().parent / "bbc_cna_plus_plus" / "bulk.seg"),
        help="Path to .seg file with BAF and BAF-se columns",
    )
    parser.add_argument(
        "--alpha", type=float, default=0.05, help="Significance level (default: 0.05)"
    )
    parser.add_argument(
        "--delta", type=float, default=0.03,
        help="Equivalence margin for TOST: |BAF - 0.5| < delta => balanced (default: 0.03)",
    )
    args = parser.parse_args()

    seg = pd.read_csv(args.seg, sep="\t")

    required = {"BAF", "BAF-se"}
    if not required.issubset(seg.columns):
        sys.exit(f"ERROR: seg file missing columns: {required - set(seg.columns)}")

    # One row per cluster: take first occurrence (values are constant within cluster)
    clusters = seg.drop_duplicates(subset=["#ID"]).copy()
    clusters = clusters[["#ID", "SAMPLE", "BAF", "BAF-se"]].reset_index(drop=True)

    diff = clusters["BAF"] - 0.5
    se = clusters["BAF-se"]

    # --- TOST equivalence: H0: |BAF - 0.5| >= delta vs H1: |BAF - 0.5| < delta ---
    # Two one-sided tests at alpha (not alpha/2, per standard TOST)
    z_upper = (diff - args.delta) / se   # test BAF < 0.5 + delta
    z_lower = (diff + args.delta) / se   # test BAF > 0.5 - delta
    p_upper = stats.norm.cdf(z_upper)    # want small (BAF well below upper bound)
    p_lower = stats.norm.sf(z_lower)     # want small (BAF well above lower bound)
    p_tost = pd.DataFrame({"pu": p_upper, "pl": p_lower}).max(axis=1)

    # --- Classification ---
    # "balanced":   p_TOST < alpha (reject H0 — evidence BAF is within delta of 0.5)
    # "imbalanced": p_TOST >= alpha (cannot conclude equivalence)
    clusters["p_tost"] = p_tost
    clusters["call"] = p_tost.apply(lambda p: "balanced" if p < args.alpha else "imbalanced")

    # --- Print ---
    print(f"alpha = {args.alpha} | TOST delta = {args.delta}")
    print(f"{'Cluster':<10} {'Sample':<14} {'BAF':>7} {'SE':>7} "
          f"{'p(TOST)':>9} {'Call':<12}")
    print("-" * 62)
    for _, r in clusters.iterrows():
        print(f"{r['#ID']:<10} {r['SAMPLE']:<14} {r['BAF']:7.4f} {r['BAF-se']:7.4f} "
              f"{r['p_tost']:9.2e} {r['call']:<12}")

    n = len(clusters)
    for label in ("balanced", "imbalanced"):
        cnt = (clusters["call"] == label).sum()
        if cnt:
            print(f"  {label}: {cnt}/{n}")


if __name__ == "__main__":
    main()
