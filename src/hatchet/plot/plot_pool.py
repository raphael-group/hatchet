"""Plot pool CNP panel from pool_instances."""

import os
import re
import logging

import matplotlib.pyplot as plt

from hatchet.utils import (
    prepare_seg_ucn,
    read_region_bed,
)
from hatchet.plot.plot_cn_utils import (
    plot_ascn_legend,
    plot_ascn_profile,
    plot_cnv_legend,
    plot_cnv_profile,
)


def _format_pool_label(tag):
    """Convert 'pool_p0.05_s1' to 'p=0.05,s=1'."""
    m = re.match(r"pool_p([^_]+)_s(\d+)", tag)
    if m:
        return f"p={m.group(1)},s={m.group(2)}"
    return tag


def _fmt_prop(v):
    """Round proportion to 2 decimals as percent; literal '0' if zero."""
    pct = round(v * 100, 2)
    return "0" if pct == 0 else f"{pct}%"


def plot_pool_cnp(
    pool_instances,
    region_bed,
    out_dir,
    sel_df=None,
    segs=None,
    title=None,
    width=20,
    height=1,
    dpi=150,
    solve_mode=None,
    sample_names=None,
    plot_ascn=True,
):
    """Plot pool CNP panel into out_dir.

    For cnt_cd: renders each solution's tree as a separate PDF.
    For cd/ilp: plots a single multi-row CNV profile panel of Pareto solutions.

    Args:
        pool_instances: {sol_id: {"imf_obj", "reg_obj", "cA", "cB", "u", ...}}.
        region_bed: Path to the whitelist region BED file.
        out_dir: Output directory for pool plots.
        sel_df: Selection DataFrame from model_select_elbow_from_regularization.
        segs: {sol_id: seg_df} pre-computed segmentation DataFrames.
        title: Optional figure title.
        width: Figure width in inches.
        height: Height in inches per profile row.
        dpi: Output resolution.
        solve_mode: "cd", "ilp", "cnt_cd", etc.
        sample_names: Sample names for cnt_cd rendering.
    """
    os.makedirs(out_dir, exist_ok=True)

    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42
    plt.rcParams["svg.fonttype"] = "none"

    regions = read_region_bed(region_bed)

    if solve_mode == "cnt_cd":
        from hatchet.plot.plot_cn_tree import render_cnt_tree

        for sol_id, sol in pool_instances.items():
            tree = sol.get("tree")
            if tree is None:
                continue
            seg_info, _, _ = prepare_seg_ucn(segs[sol_id])
            out_file = os.path.join(out_dir, f"{sol_id}.pdf")
            render_cnt_tree(
                tree,
                seg_info,
                regions,
                out_file,
                u=sol.get("u"),
                sample_names=sample_names,
            )
        return

    selected_ids = set()
    if sel_df is not None:
        selected_ids = set(
            sel_df.loc[sel_df["selected"] == "*", "instance_id"].tolist()
        )

    pareto_ids = []
    if sel_df is not None:
        pareto_ids = sel_df.loc[sel_df["is_pareto"], "instance_id"].tolist()
    else:
        pareto_ids = sorted(pool_instances.keys())

    entries = []
    for sol_id in pareto_ids:
        sol = pool_instances[sol_id]
        is_selected = sol_id in selected_ids
        entries.append((sol_id, segs[sol_id], sol["imf_obj"], is_selected))
    entries.sort(key=lambda x: x[2])

    if not entries:
        logging.warning(f"plot_pool_cnp: no solutions to plot, skipping {out_dir}")
        return

    nrows = len(entries)
    first_seg_df = entries[0][1]
    n_clones = len([c for c in first_seg_df.columns if c.startswith("cn_")])
    n_samples = first_seg_df["SAMPLE"].nunique()
    row_h = height * max(1, n_clones - 1) + 0.2 * max(0, n_samples - 1)
    fig, axes = plt.subplots(
        nrows=nrows + 1,
        ncols=1,
        figsize=(width, row_h * nrows),
        gridspec_kw={"height_ratios": [row_h] * nrows + [2 * height]},
    )
    fig.subplots_adjust(hspace=0.6 + 0.1 * max(0, n_samples - 1))
    main_axes = axes[:-1] if nrows > 1 else [axes[0]]
    ax_leg = axes[-1]

    for i, (label, seg_df, obj, is_selected) in enumerate(entries):
        seg_info_all, clones, _ = prepare_seg_ucn(seg_df)
        samples = seg_info_all["SAMPLE"].unique().tolist()
        seg_info = seg_info_all.loc[
            seg_info_all["SAMPLE"] == samples[0], :
        ].reset_index(drop=True)

        _profile_fn = plot_ascn_profile if plot_ascn else plot_cnv_profile
        _profile_fn(
            main_axes[i],
            seg_info,
            regions,
            plot_chrname=True,
            width=width,
            height=height,
            show_clone_name=False,
            show_prop=False,
        )

        short_label = _format_pool_label(str(label))
        if is_selected:
            short_label += " *"
        prop_lines = []
        for sid in samples:
            sp = seg_info_all.loc[seg_info_all["SAMPLE"] == sid, :].reset_index(
                drop=True
            )
            cps = sp[[f"u_{c}" for c in clones]].iloc[0].tolist()
            prop_lines.append(f"{sid}:" + "|".join(_fmt_prop(c) for c in cps))
        ylabel = f"{short_label}\nimf {round(obj, 2)}\n" + "\n".join(prop_lines)
        color = "red" if is_selected else "black"
        main_axes[i].set_ylabel(
            ylabel, rotation=0, ha="right", va="center", color=color
        )

    _legend_fn = plot_ascn_legend if plot_ascn else plot_cnv_legend
    _legend_fn(ax_leg)

    if title:
        main_axes[0].set_title(title)
    out_file = os.path.join(out_dir, "pool.pdf")
    plt.savefig(out_file, dpi=dpi, bbox_inches="tight")
    plt.close()
    logging.info(f"pool CNP panel saved to {out_file}")
