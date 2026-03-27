import os
import re
import glob
import logging

import numpy as np
import pandas as pd

from hatchet.evaluate.vaf_utils import estimate_vaf, is_explained_mut, relative_error


def read_snv_vcf(vcf_path, sample="tumor"):
    """Read somatic SNVs from a VCF file via cyvcf2."""
    try:
        from cyvcf2 import VCF
    except ImportError:
        raise ImportError(
            "cyvcf2 is required for VCF input. Use --snv_tsv for TSV input."
        )

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

        rows.append(
            (rec.CHROM, rec.POS, rec.REF, alt, ref_reads, alt_reads, dp_val, vaf_val)
        )

    return pd.DataFrame(
        rows,
        columns=[
            "#CHR",
            "POS",
            "REF",
            "ALT",
            "ref_reads",
            "alt_reads",
            "DP",
            "observed_VAF",
        ],
    )


def read_snv_tsv(tsv_path):
    """Read somatic SNVs from a TSV file.

    Supported column names (auto-mapped):
      CHR / #CHR, POS, REF_COUNT / ref_reads, VAR_COUNT / alt_reads,
      VAF / observed_VAF.
    """
    df = pd.read_csv(tsv_path, sep="\t")
    col_map = {
        "CHR": "#CHR",
        "REF_COUNT": "ref_reads",
        "VAR_COUNT": "alt_reads",
        "VAF": "observed_VAF",
        "DEPTH": "DP",
    }
    df = df.rename(columns={k: v for k, v in col_map.items() if k in df.columns})
    if "#CHR" in df.columns:
        df["#CHR"] = df["#CHR"].apply(
            lambda x: x if str(x).startswith("chr") else f"chr{x}"
        )
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
        idx = np.searchsorted(starts, positions, side="right") - 1
        valid = (
            (idx >= 0)
            & (idx < len(starts))
            & (positions < ends[np.clip(idx, 0, len(ends) - 1)])
        )
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
        n_mutated_clones = int((np.array(best[1:]) > 0).sum())
        is_subclonal = 0 < n_mutated_clones < (len(clones) - 1)

        dp = ref_reads + alt_reads
        results.append(
            {
                "#CHR": row["#CHR"],
                "POS": row["POS"],
                "ref_reads": ref_reads,
                "alt_reads": alt_reads,
                "DP": dp,
                "observed_VAF": round(obs_vaf, 4),
                "predicted_VAF": round(pred_vaf, 4),
                "relative_error": round(rel_err, 4) if np.isfinite(rel_err) else np.nan,
                "mutated_copies": ",".join(map(str, best)),
                "allele": allele,
                "CCF": round(ccf, 4),
                "CNP": seg["CNP"],
                "is_subclonal": is_subclonal,
                "is_explained": explained,
            }
        )

    return pd.DataFrame(results)


def _prepare_seg_for_cnp(segs):
    """Build CNP and PROPS columns from cn_*/u_* columns in seg.ucn format."""
    cn_cols = sorted(
        [c for c in segs.columns if c.startswith("cn_")],
        key=lambda c: (0 if c == "cn_normal" else 1, c),
    )
    u_cols = sorted(
        [c for c in segs.columns if c.startswith("u_")],
        key=lambda c: (0 if c == "u_normal" else 1, c),
    )
    seg_plot = segs.copy()
    seg_plot["CNP"] = seg_plot.apply(
        lambda r: ";".join(str(r[c]) for c in cn_cols), axis=1
    )
    seg_plot["PROPS"] = seg_plot.apply(
        lambda r: ";".join(str(r[c]) for c in u_cols), axis=1
    )
    return seg_plot


def _compute_abs_coords(regions, chrs):
    """Compute absolute coordinates using the same system as plot_cnv_profile (start from 0)."""
    regions_chs = regions.groupby(by="#CHR", sort=False)
    ch_offset = 0
    ch_coords = []
    seg_coords = []
    # For each chromosome, map (wl_start, wl_end) → (abs_start, abs_end)
    intervals = {}  # {chr: [(wl_start, wl_end, abs_start, abs_end), ...]}
    for ch in chrs:
        ch_coords.append(ch_offset)
        regions_ch = regions_chs.get_group(ch)
        ch_intervals = []
        for si in range(len(regions_ch)):
            wl = regions_ch.iloc[si]
            wl_start, wl_end = wl["START"], wl["END"]
            abs_start = ch_offset
            abs_end = ch_offset + (wl_end - wl_start)
            ch_intervals.append((wl_start, wl_end, abs_start))
            ch_offset = abs_end
            if si < len(regions_ch) - 1:
                seg_coords.append(ch_offset)
        intervals[ch] = ch_intervals
    ch_coords.append(ch_offset)
    return intervals, ch_coords, seg_coords, ch_offset


def plot_vaf_1d(
    snv_result, segs, chrom_sizes, regions, out_file, dpi=500, transparent=False
):
    """Plot observed VAF along the genome with segment-level expected VAF bars and CNP profile."""
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    from hatchet.plot.plot_cn_utils import plot_cnv_profile, plot_cnv_legend

    seg_chrs = segs["#CHR"].unique().tolist()

    snv = snv_result.copy()
    snv = snv[snv["#CHR"].isin(seg_chrs)].reset_index(drop=True)
    chr_rank = {c: i for i, c in enumerate(seg_chrs)}
    snv["_chr_rank"] = snv["#CHR"].map(chr_rank)
    snv = (
        snv.sort_values(["_chr_rank", "POS"])
        .drop(columns=["_chr_rank"])
        .reset_index(drop=True)
    )

    # Compute absolute coordinates matching plot_cnv_profile
    intervals, ch_coords, seg_coords, genome_end = _compute_abs_coords(
        regions, seg_chrs
    )

    def pos_to_abs(chrom, pos):
        for wl_start, wl_end, abs_start in intervals.get(chrom, []):
            if wl_start <= pos < wl_end:
                return abs_start + (pos - wl_start)
        return None

    snv["abs_pos"] = snv.apply(lambda r: pos_to_abs(r["#CHR"], r["POS"]), axis=1)
    snv = snv.dropna(subset=["abs_pos"]).reset_index(drop=True)

    # Color: clonal explained (blue), subclonal explained (green), unexplained (red)
    def _snv_color(row):
        if not row["is_explained"]:
            return "red"
        return "seagreen" if row["is_subclonal"] else "steelblue"

    colors = snv.apply(_snv_color, axis=1).to_numpy()

    fig, axes = plt.subplots(
        nrows=3,
        ncols=1,
        figsize=(20, 6),
        gridspec_kw={"height_ratios": [3, 2, 1]},
    )
    ax_vaf, ax_cnp, ax_leg = axes

    # Scatter observed VAF
    ax_vaf.scatter(
        snv["abs_pos"],
        snv["observed_VAF"],
        s=12,
        c=colors,
        edgecolors="none",
        linewidths=0,
        zorder=2,
    )

    # Legend
    from matplotlib.lines import Line2D

    legend_handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            markerfacecolor="steelblue",
            markersize=6,
            label="clonal explained",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            markerfacecolor="seagreen",
            markersize=6,
            label="subclonal explained",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            markerfacecolor="red",
            markersize=6,
            label="unexplained",
        ),
    ]
    ax_vaf.legend(
        handles=legend_handles,
        loc="upper left",
        bbox_to_anchor=(1.01, 1),
        frameon=False,
        fontsize=8,
    )

    # Draw segment-level expected VAF bars
    exp_col = "predicted_VAF" if "predicted_VAF" in segs.columns else None
    if exp_col is not None:
        seg_lines = []
        seg_lines_comp = []
        for _, seg in segs.iterrows():
            ch = seg["#CHR"]
            for wl_start, wl_end, abs_start in intervals.get(ch, []):
                s0 = max(seg["START"], wl_start)
                s1 = min(seg["END"], wl_end)
                if s0 >= s1:
                    continue
                x0 = abs_start + (s0 - wl_start)
                x1 = abs_start + (s1 - wl_start)
                baf = seg[exp_col]
                seg_lines.append([(x0, baf), (x1, baf)])
                seg_lines_comp.append([(x0, 1.0 - baf), (x1, 1.0 - baf)])
        if seg_lines:
            ax_vaf.add_collection(
                LineCollection(
                    seg_lines,
                    linewidth=1.5,
                    colors=[(0, 0, 0, 1)] * len(seg_lines),
                    zorder=1,
                )
            )
            ax_vaf.add_collection(
                LineCollection(
                    seg_lines_comp,
                    linewidth=1.5,
                    colors=[(0, 0, 0, 0.4)] * len(seg_lines_comp),
                    zorder=1,
                )
            )

    # Chromosome boundaries and labels
    ax_vaf.vlines(
        ch_coords,
        ymin=0,
        ymax=1,
        transform=ax_vaf.get_xaxis_transform(),
        linewidth=1,
        colors="k",
    )
    for sc in seg_coords:
        ax_vaf.axvline(sc, color="k", linewidth=1, linestyle="dashed", ymin=0, ymax=1)

    ax_vaf.set_xlim(0, genome_end)
    ax_vaf.set_ylim(-0.01, 1.05)
    ax_vaf.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax_vaf.set_ylabel("VAF")
    ax_vaf.set_xticks([])
    ax_vaf.grid(False)

    # CNP profile below (same coordinate system, starts from 0)
    seg_plot = _prepare_seg_for_cnp(segs)
    plot_cnv_profile(
        ax_cnp,
        seg_plot,
        regions,
        width=20,
        height=1,
        plot_chrname=True,
        show_clone_name=True,
        show_prop=True,
    )
    plot_cnv_legend(ax_leg)

    n_explained = int(snv["is_explained"].sum())
    n_total = len(snv)
    fig.suptitle(
        f"SNV VAF evaluation: {n_explained}/{n_total} explained ({n_explained / n_total:.1%})"
    )
    plt.tight_layout()
    fig.savefig(out_file, dpi=dpi, bbox_inches="tight", transparent=transparent)
    plt.close(fig)
    logging.info(f"wrote {out_file}")


def pool_sol_to_seg(pool_tsv, bbc_df):
    """Convert a cluster-level pool solution TSV to a seg-level DataFrame.

    Joins the pool's CN/u columns with the BBC's per-bin genomic coordinates
    via cluster ID, producing a DataFrame compatible with evaluate_snvs.
    """
    pool = pd.read_csv(pool_tsv, sep="\t")
    cn_cols = sorted(
        [c for c in pool.columns if c.startswith("cn_")],
        key=lambda c: (0 if c == "cn_normal" else 1, c),
    )
    u_cols = sorted(
        [c for c in pool.columns if c.startswith("u_")],
        key=lambda c: (0 if c == "u_normal" else 1, c),
    )
    clones = [c.replace("cn_", "") for c in cn_cols]

    # Build CNP string per cluster per sample
    pool["CNP"] = pool.apply(lambda row: ";".join(str(row[c]) for c in cn_cols), axis=1)

    # Merge with BBC to get per-bin coordinates
    merged = bbc_df[["#CHR", "START", "END", "SAMPLE", "CLUSTER"]].merge(
        pool[["CLUSTER", "SAMPLE", "CNP"] + u_cols],
        on=["CLUSTER", "SAMPLE"],
        how="inner",
    )
    return merged, clones, pool, cn_cols, u_cols


def evaluate_pool_solutions(result_dir, snv_df, gamma=0.05):
    """Evaluate all pool solutions in result_dir against SNVs.

    Returns a DataFrame with one row per pool solution: ploidy, n_clones, tag,
    IMF/REG objectives (from summary.tsv), and VAF evaluation metrics.
    """
    summary_path = os.path.join(result_dir, "summary.tsv")
    bbc_path = os.path.join(result_dir, "bulk.good.bbc")
    sols_dir = os.path.join(result_dir, "sols")

    summary = pd.read_csv(summary_path, sep="\t")
    bbc_df = pd.read_csv(bbc_path, sep="\t")
    samples = bbc_df["SAMPLE"].unique().tolist()

    eval_rows = []
    for _, srow in summary.iterrows():
        ploidy = srow["ploidy"]
        n_clones = int(srow["n_clones"])
        tag = srow["tag"]

        # Parse tag to find the sol file: pool_p{pparam}_s{idx}
        m = re.match(r"pool_p([\d.]+)_s(\d+)", tag)
        if not m:
            continue
        pparam, pidx = m.group(1), m.group(2)
        sol_subdir = os.path.join(sols_dir, f"{ploidy}_n{n_clones}")
        sol_file = os.path.join(sol_subdir, f"cd_sol{pparam}_pool{pidx}.tsv")
        if not os.path.exists(sol_file):
            logging.warning(f"sol file not found: {sol_file}")
            continue

        merged, clones, pool, cn_cols, u_cols = pool_sol_to_seg(sol_file, bbc_df)

        # Evaluate per sample
        n_total = n_explained = 0
        total_err = 0.0
        for sample in samples:
            segs_s = merged[merged["SAMPLE"] == sample].reset_index(drop=True)
            if len(segs_s) == 0:
                continue

            pool_s = pool[pool["SAMPLE"] == sample].iloc[0]
            clone_props = np.array([float(pool_s[uc]) for uc in u_cols])
            purity = float(np.sum(clone_props[1:]))
            if purity <= 1e-8:
                continue

            snv_mapped = assign_snvs_to_segments(snv_df, segs_s)
            result_df = evaluate_snvs(
                segs_s, clones, clone_props, snv_mapped, gamma=gamma
            )
            n_total += len(result_df)
            n_explained += (
                int(result_df["is_explained"].sum()) if len(result_df) > 0 else 0
            )
            if len(result_df) > 0:
                total_err += result_df["relative_error"].sum()

        ratio = n_explained / n_total if n_total > 0 else 0.0
        mean_err = total_err / n_total if n_total > 0 else np.nan

        row = {
            "ploidy": ploidy,
            "n_clones": n_clones,
            "tag": tag,
            "n_snvs": n_total,
            "n_explained": n_explained,
            "explained_ratio": round(ratio, 4),
            "mean_relative_error": round(mean_err, 4)
            if np.isfinite(mean_err)
            else np.nan,
        }
        # Copy objective columns from summary
        for col in summary.columns:
            if col not in row:
                row[col] = srow[col]
        eval_rows.append(row)

    return pd.DataFrame(eval_rows)
