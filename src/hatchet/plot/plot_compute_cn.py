"""Plotting for the compute-cn command.

Solution-pool CNP panels, the scaling-factor 2D diagnostic, and clone-tree
rendering — everything the compute-cn pipeline draws.
"""

from __future__ import annotations

import os
import re
import logging

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.collections import LineCollection

from cnplot import annotate_landmarks, plot_cnv_profile, plot_scatter_2d, set_palette
from hatchet.utils import sort_df_chr
from hatchet.plot import plot_cn as _plot_cn
from hatchet.plot.plot_utils import build_genome_axis, use_editable_fonts


def _clones_from_cn(df):
    """Ordered clone names ("normal", "clone1", ...) from a df's cn_ columns."""
    n = len([c for c in df.columns if c.startswith("cn_")])
    return ["normal"] + [f"clone{i}" for i in range(1, n)]


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
    genome_size,
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
    out_name="pool.pdf",
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
        genome_size: Chromosome-sizes path for the genome axis.
    """
    os.makedirs(out_dir, exist_ok=True)

    use_editable_fonts()
    # Restrict the axis to chromosomes present in the solutions' segments.
    keep = set()
    for sdf in (segs or {}).values():
        keep.update(sdf["#CHR"].unique())
    genome_axis = build_genome_axis(region_bed, genome_size, keep_chroms=keep or None)

    if solve_mode == "cnt_cd":
        for sol_id, sol in pool_instances.items():
            tree = sol.get("tree")
            if tree is None:
                continue
            seg_info = sort_df_chr(segs[sol_id].copy(), pos="START")
            out_file = os.path.join(out_dir, f"{sol_id}.pdf")
            render_cnt_tree(
                tree,
                seg_info,
                genome_axis,
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
        seg_info_all = sort_df_chr(seg_df.copy(), pos="START")
        clones = _clones_from_cn(seg_info_all)
        samples = seg_info_all["SAMPLE"].unique().tolist()
        seg_info = seg_info_all.loc[
            seg_info_all["SAMPLE"] == samples[0], :
        ].reset_index(drop=True)

        plot_cnv_profile(
            main_axes[i],
            seg_info,
            genome_axis,
            ax_leg=(ax_leg if i == nrows - 1 else None),
            plot_chrname=True,
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

    if title:
        main_axes[0].set_title(title)
    out_file = os.path.join(out_dir, out_name)
    plt.savefig(out_file, dpi=dpi, bbox_inches="tight")
    plt.close()
    logging.info(f"pool CNP panel saved to {out_file}")


def plot_scaling_2d(
    samples: list,
    bbcs: pd.DataFrame,
    segs: pd.DataFrame,
    scaling: dict,
    out_file: str,
    markersize: float = 3.0,
    markersize_centroid: float = 14,
    dpi: int = 300,
    transparent: bool = False,
    maxlim_rdr: int = 10,
):
    """2D RDR-vs-BAF scatter anchoring the scaling inference from get_scaling_factor.

    One PDF page per sample and WGD mode (noWGD = diploid, WGD = tetraploid,
    drawn only when a WGD scaling was inferred). Bins are colored per cluster via
    cnplot's plot_scatter_2d; anchor clusters are circled and annotated with their
    inferred (a, b) clonal states via annotate_landmarks.
    """
    use_editable_fonts()

    clusters = sorted(bbcs["CLUSTER"].unique().tolist())
    palette = set_palette(num_colors=len(clusters))
    pal = {str(c): palette[i] for i, c in enumerate(clusters)}

    panels = [("noWGD", scaling["diploid"])]
    if scaling.get("tetraploid") is not None:
        panels.append(("WGD", scaling["tetraploid"]))

    cent = segs.set_index(["#ID", "SAMPLE"])[["RD", "BAF"]]

    pdf = PdfPages(out_file)
    for sample in samples:
        sub = bbcs[bbcs["SAMPLE"] == sample]
        obs = pd.DataFrame({"BAF": sub["BAF"].to_numpy(), "RD": sub["RD"].to_numpy()})
        obs["CLUSTER"] = sub["CLUSTER"].astype(str).to_numpy()

        lim_baf = (0, 1) if obs["BAF"].max() > 0.5 else (0, 0.55)
        lim_rdr = (0, min(max(2, int(np.ceil(obs["RD"].max()))), maxlim_rdr))

        for label, info in panels:
            landmarks = []
            for c, (a, b) in info["clonal"].items():
                if (c, sample) not in cent.index:
                    continue
                landmarks.append(
                    {
                        "x": cent.loc[(c, sample), "BAF"],
                        "y": cent.loc[(c, sample), "RD"],
                        "label": f"({a},{b})",
                        "clonal": True,
                    }
                )

            p = (info.get("purities") or {}).get(sample)
            ptxt = f"  purity={p:.3f}" if p is not None else ""
            grid = plot_scatter_2d(
                obs,
                xcol="BAF",
                ycol="RD",
                hue="CLUSTER",
                palette=pal,
                xlim=lim_baf,
                ylim=lim_rdr,
                xlabel="BAF",
                ylabel="RDR",
                title=f"sample={sample}  {label}{ptxt}",
                refline_x=0.5,
                markersize=markersize,
            )
            annotate_landmarks(grid.ax_joint, landmarks, markersize=markersize_centroid)
            pdf.savefig(
                grid.figure, dpi=dpi, bbox_inches="tight", transparent=transparent
            )
            plt.close(grid.figure)
    pdf.close()
    logging.info(f"scaling 2D scatter saved to {out_file}")


def _display_name(v, tree):
    if v == tree.normal_leaf or v == tree.root:
        return "normal"
    if v in set(tree.tumor_leaves):
        return f"clone{tree.tumor_leaves.index(v) + 1}"
    return f"v{v}"


def _inorder(node, ch):
    if node not in ch:
        return [node]
    lc, r = ch[node]
    return _inorder(lc, ch) + [node] + _inorder(r, ch)


def _compute_x_weighted(node, ch, ec, x_parent=0):
    pos = {node: x_parent}
    if node in ch:
        for c in ch[node]:
            pos.update(_compute_x_weighted(c, ch, ec, x_parent + ec.get((node, c), 1)))
    return pos


def _draw_sample_tree(
    ax,
    tree,
    tumor_children,
    tumor_subtree_root,
    leaf_set,
    x_weighted,
    node_y,
    edge_cost,
    max_leaf_x,
    inactive,
    sample_name,
    normal_prop,
    leaf_props,
    fontsize,
):
    """Draw a per-sample tree with inactive clones dotted."""
    tumor_y_scale, tumor_y_offset = 0.75, 0.15

    def _sy(y):
        return tumor_y_offset + y * tumor_y_scale

    def _is_inactive(node):
        if node in leaf_set:
            return node in inactive
        if node not in tumor_children:
            return False
        lc, r = tumor_children[node]
        return _is_inactive(lc) and _is_inactive(r)

    tumor_root_sy = _sy(node_y[tumor_subtree_root])
    normal_y = tumor_y_offset + tumor_y_scale + 0.06
    root_y = (tumor_root_sy + normal_y) / 2
    root_x = -1

    ax.set_xlim(root_x - 0.5, max_leaf_x + 3.5)
    ax.set_ylim(-0.05, 1.05)
    ax.axis("off")

    # Root -> normal (top)
    ax.plot([root_x, root_x], [root_y, normal_y], "k-", lw=2)
    ax.plot([root_x, root_x + 0.5], [normal_y, normal_y], "k-", lw=2)
    ax.text(
        root_x + 0.6,
        normal_y,
        f"normal ({normal_prop:.2%})",
        ha="left",
        va="center",
        fontsize=fontsize,
    )
    # Root -> tumor subtree
    ax.plot([root_x, root_x], [root_y, tumor_root_sy], "k-", lw=2)
    ax.plot(
        [root_x, x_weighted[tumor_subtree_root]],
        [tumor_root_sy, tumor_root_sy],
        "k-",
        lw=2,
    )
    ax.text(
        root_x,
        root_y,
        "normal",
        ha="center",
        va="center",
        fontsize=fontsize,
        bbox=dict(
            boxstyle="round,pad=0.15", facecolor="white", edgecolor="black", lw=0.5
        ),
        zorder=10,
    )

    for parent, (lc, r) in tumor_children.items():
        px, py = x_weighted[parent], _sy(node_y[parent])
        for c in (lc, r):
            cx, cy = x_weighted[c], _sy(node_y[c])
            style = ":" if _is_inactive(c) else "-"
            color = "#BBBBBB" if _is_inactive(c) else "black"
            ax.plot([px, px], [py, cy], linestyle=style, color=color, lw=2)
            ax.plot([px, cx], [cy, cy], linestyle=style, color=color, lw=2)
            cost = edge_cost.get((parent, c), "")
            if cost:
                ax.text(
                    (px + cx) / 2,
                    cy + 0.008,
                    str(cost),
                    ha="center",
                    va="bottom",
                    fontsize=fontsize - 1,
                    color=color,
                )

    for v in x_weighted:
        vx, vy = x_weighted[v], _sy(node_y[v])
        if v in leaf_set:
            li = tree.tumor_leaves.index(v)
            prop = leaf_props[li]
            color = "#BBBBBB" if v in inactive else "black"
            ax.text(
                vx + 0.15,
                vy,
                f"{_display_name(v, tree)} ({prop:.2%})",
                ha="left",
                va="center",
                fontsize=fontsize,
                color=color,
            )
        else:
            color = "#BBBBBB" if _is_inactive(v) else "black"
            ax.text(
                vx,
                vy,
                _display_name(v, tree),
                ha="center",
                va="center",
                fontsize=fontsize,
                color=color,
                bbox=dict(
                    boxstyle="round,pad=0.15",
                    facecolor="white",
                    edgecolor=color,
                    lw=0.5,
                ),
                zorder=10,
            )

    ax.set_title(sample_name, fontsize=fontsize + 2)


def render_cnt_tree(
    labeled_tree,
    bin_info,
    genome_axis,
    out_path,
    u=None,
    sample_names=None,
    min_prop=0.03,
    fontsize=10,
):
    """Render a LabeledCloneTree as a multi-page PDF.

    Page 1: dendrogram + CNP profile + legend.
    Page 2 (if u provided): per-sample trees with inactive clones dotted.

    Args:
        labeled_tree: LabeledCloneTree with inferred CN and events.
        bin_info: DataFrame with #CHR, START, END, and cn_clone* columns.
        genome_axis: cnplot GenomeAxis to draw the profile on.
        out_path: output PDF path.
        u: (n_leaves, P) usage matrix. If provided, page 2 is generated.
        sample_names: list of P sample names.
        min_prop: clones with usage <= min_prop shown as dotted.
        fontsize: font size for labels.
    """
    tree = labeled_tree
    tumor_children = {k: v for k, v in tree.children.items() if k != tree.root}
    tumor_subtree_root = [c for c in tree.children[tree.root] if c != tree.normal_leaf][
        0
    ]
    leaf_set = set(tree.tumor_leaves)

    cnp_nodes = _inorder(tumor_subtree_root, tumor_children)
    num_rows = len(cnp_nodes)
    h = 1.0 / num_rows
    node_y = {cnp_nodes[i]: (num_rows - 1 - i + 0.5) * h for i in range(num_rows)}

    te_map = {(p, c): i for i, (p, c) in enumerate(tree.tumor_edges)}
    edge_cost = {}
    for p, c in tree.edges:
        edge_cost[(p, c)] = tree._edge_cost(te_map[(p, c)]) if (p, c) in te_map else 0

    x_weighted = _compute_x_weighted(tumor_subtree_root, tumor_children, edge_cost, 0)
    max_leaf_x = max((x_weighted[v] for v in tree.tumor_leaves), default=1) or 1

    # Reverse-inorder clone order maps each node to the same row plot_cnv_profile
    # drew it in under the former CNP layout (entries stacked bottom-to-top).
    cnp_rev = list(reversed(cnp_nodes))
    node_clones = [f"n{v}" for v in cnp_rev]
    leaf_cn_to_seg = {}
    for s in range(tree.a_all.shape[0]):
        key = tuple(
            (int(tree.a_all[s, v]), int(tree.b_all[s, v])) for v in tree.tumor_leaves
        )
        leaf_cn_to_seg[key] = s
    clone_cols = [f"cn_clone{i}" for i in range(1, len(tree.tumor_leaves) + 1)]
    bin_info = bin_info.copy()
    seg_idx = []
    for _, row in bin_info.iterrows():
        key = tuple(
            (int(v.split("|")[0]), int(v.split("|")[1]))
            for v in (row[c] for c in clone_cols)
        )
        seg_idx.append(leaf_cn_to_seg[key])
    seg_idx = np.asarray(seg_idx)
    for v in cnp_nodes:
        a = tree.a_all[seg_idx, v].astype(int)
        b = tree.b_all[seg_idx, v].astype(int)
        bin_info[f"cn_n{v}"] = [f"{ai}|{bi}" for ai, bi in zip(a, b)]

    with PdfPages(out_path) as pdf:
        # Page 1
        fig = plt.figure(figsize=(26, max(8, num_rows * 0.9)))
        gs = fig.add_gridspec(
            2,
            2,
            width_ratios=[1.5, 5],
            height_ratios=[5, 1],
            wspace=0.08,
            hspace=0.3,
        )
        ax_tree = fig.add_subplot(gs[0, 0])
        ax_cn = fig.add_subplot(gs[0, 1])
        ax_leg = fig.add_subplot(gs[1, :])

        plot_cnv_profile(
            ax_cn,
            bin_info,
            genome_axis,
            ax_leg=ax_leg,
            clones=node_clones,
            show_prop=False,
        )
        ax_cn.set_yticks([(num_rows - 1 - i + 0.5) * h for i in range(num_rows)])
        ax_cn.set_yticklabels(
            [_display_name(v, tree) for v in cnp_nodes], fontsize=fontsize
        )
        ax_cn.tick_params(axis="x", labelsize=fontsize - 2)

        for child in ax_cn.get_children():
            if isinstance(child, LineCollection):
                segs = child.get_segments()
                for seg in segs:
                    if len(seg) == 2 and seg[0][0] == seg[1][0] and seg[1][1] > 1.05:
                        seg[1][1] = 1.05
                child.set_segments(segs)

        ax_tree.set_xlim(-0.5, max_leaf_x + 2.0)
        ax_tree.set_ylim(0, 1)
        ax_tree.axis("off")
        for parent, (lc, r) in tumor_children.items():
            px, py = x_weighted[parent], node_y[parent]
            for c in (lc, r):
                cx, cy = x_weighted[c], node_y[c]
                ax_tree.plot([px, px], [py, cy], "k-", lw=2)
                ax_tree.plot([px, cx], [cy, cy], "k-", lw=2)
                cost = edge_cost.get((parent, c), "")
                if cost:
                    ax_tree.text(
                        (px + cx) / 2,
                        cy + 0.008,
                        str(cost),
                        ha="center",
                        va="bottom",
                        fontsize=fontsize,
                    )
        for v in x_weighted:
            if v not in leaf_set:
                ax_tree.text(
                    x_weighted[v],
                    node_y[v],
                    _display_name(v, tree),
                    ha="center",
                    va="center",
                    fontsize=fontsize,
                    bbox=dict(
                        boxstyle="round,pad=0.15",
                        facecolor="white",
                        edgecolor="black",
                        lw=0.5,
                    ),
                    zorder=10,
                )

        pdf.savefig(fig, bbox_inches="tight", dpi=150)
        plt.close(fig)

        # Page 2
        if u is not None:
            u_arr = np.array(u)
            P = u_arr.shape[1]
            if sample_names is None:
                sample_names = [f"Sample_{i}" for i in range(P)]

            fig2 = plt.figure(figsize=(12, 5 * P))
            gs2 = fig2.add_gridspec(P, 1, hspace=0.4)

            for pi in range(P):
                ax = fig2.add_subplot(gs2[pi, 0])
                inactive = {
                    tree.tumor_leaves[li]
                    for li in range(len(tree.tumor_leaves))
                    if u_arr[li + 1, pi] <= min_prop
                }
                _draw_sample_tree(
                    ax,
                    tree,
                    tumor_children,
                    tumor_subtree_root,
                    leaf_set,
                    x_weighted,
                    node_y,
                    edge_cost,
                    max_leaf_x,
                    inactive,
                    sample_names[pi],
                    u_arr[0, pi],
                    [u_arr[li + 1, pi] for li in range(len(tree.tumor_leaves))],
                    fontsize,
                )

            pdf.savefig(fig2, bbox_inches="tight", dpi=150)
            plt.close(fig2)


def run_plot_cn(args, bbc, seg, gamma_file, plot_dir, ploidy, name=None):
    """Auto-run plot-cn on a compute-cn solution; styling falls through to hatchet.yaml."""
    if not os.path.exists(bbc) or not os.path.exists(seg):
        return
    _plot_cn.run(
        {
            "bbc": bbc,
            "seg": seg,
            "genome_size": args["genome_size"],
            "region_bed": args["region_bed"],
            "gamma_file": gamma_file,
            "solfile": None,
            "patient_id": name,
            "plot_dir": plot_dir,
            "ploidy": ploidy,
        }
    )


def plot_pareto_curve(summary_df, plot_dir, reg_term, elbow_fig=None):
    """Plot REG vs IMF Pareto curves + elbow/BIC page as a multi-page PDF."""
    outfile = os.path.join(plot_dir, "model_selection.pdf")
    reg_col = reg_term if reg_term in summary_df.columns else "REG"
    ploidies = sorted(summary_df["ploidy"].unique())
    cmap = plt.get_cmap("tab10")

    n_pages = 0
    with PdfPages(outfile) as pdf:
        # One page per ploidy; overlay all n-clone solutions, each n a distinct color.
        for ploidy in ploidies:
            pdf_grp = summary_df[summary_df["ploidy"] == ploidy]
            ns = sorted(pdf_grp["n_clones"].unique())
            fig, ax = plt.subplots(figsize=(7, 5))

            # Non-pareto across all n: shared light-gray backdrop
            non_pareto = pdf_grp[~pdf_grp["is_pareto"]]
            if len(non_pareto) > 0:
                ax.scatter(
                    non_pareto[reg_col],
                    non_pareto["IMF"],
                    c="0.8",
                    s=15,
                    zorder=2,
                    alpha=0.4,
                    linewidths=0,
                )

            for ni, n_clones in enumerate(ns):
                grp = pdf_grp[pdf_grp["n_clones"] == n_clones]
                color = cmap(ni % 10)
                pareto = grp[grp["is_pareto"]].sort_values(reg_col)
                if len(pareto) > 0:
                    ax.plot(
                        pareto[reg_col],
                        pareto["IMF"],
                        "-o",
                        color=color,
                        markersize=5,
                        linewidth=1.3,
                        zorder=4,
                        label=f"n={n_clones}",
                    )
                sel = grp[grp["selected"] == "*"]
                if len(sel) > 0:
                    ax.scatter(
                        sel[reg_col],
                        sel["IMF"],
                        facecolors=color,
                        marker="*",
                        s=250,
                        zorder=5,
                        edgecolors="black",
                        linewidths=1,
                    )

            ax.set_xlabel(reg_col, fontsize=11)
            ax.set_ylabel("IMF", fontsize=11)
            ax.set_title(
                f"{ploidy} ({len(pdf_grp)} solutions)",
                fontsize=13,
                fontweight="bold",
            )
            ax.legend(fontsize=9, title="clones")
            ax.grid(True, alpha=0.3)
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)
            n_pages += 1

        # Append elbow/BIC figure as last page
        if elbow_fig is not None:
            pdf.savefig(elbow_fig)
            plt.close(elbow_fig)
            n_pages += 1

    logging.info(f"wrote {outfile} ({n_pages} pages)")
