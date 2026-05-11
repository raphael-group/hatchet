"""CNT clone-tree + CNP profile rendering."""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.collections import LineCollection

from hatchet.plot.plot_cn_utils import (
    plot_ascn_legend,
    plot_ascn_profile,
    plot_cnv_legend,
    plot_cnv_profile,
)


def _display_name(v, tree):
    if v == tree.normal_leaf or v == tree.root:
        return "normal"
    if v in set(tree.tumor_leaves):
        return f"clone{tree.tumor_leaves.index(v) + 1}"
    return f"v{v}"


def _inorder(node, ch):
    if node not in ch:
        return [node]
    l, r = ch[node]
    return _inorder(l, ch) + [node] + _inorder(r, ch)


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
        l, r = tumor_children[node]
        return _is_inactive(l) and _is_inactive(r)

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

    for parent, (l, r) in tumor_children.items():
        px, py = x_weighted[parent], _sy(node_y[parent])
        for c in (l, r):
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
    regions,
    out_path,
    u=None,
    sample_names=None,
    min_prop=0.03,
    fontsize=10,
    plot_ascn=True,
):
    """Render a LabeledCloneTree as a multi-page PDF.

    Page 1: dendrogram + CNP profile + legend.
    Page 2 (if u provided): per-sample trees with inactive clones dotted.

    Args:
        labeled_tree: LabeledCloneTree with inferred CN and events.
        bin_info: DataFrame with #CHR, START, END, CNP, PROPS.
        regions: DataFrame with #CHR, START, END.
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

    # Rebuild CNP to include all tree nodes (leaves + internal) in reverse
    # inorder, since plot_cnv_profile renders entries bottom-to-top.
    cnp_rev = list(reversed(cnp_nodes))
    leaf_cn_to_seg = {}
    for s in range(tree.a_all.shape[0]):
        key = tuple(
            (int(tree.a_all[s, v]), int(tree.b_all[s, v])) for v in tree.tumor_leaves
        )
        leaf_cn_to_seg[key] = s
    clone_cols = [f"cn_clone{i}" for i in range(1, len(tree.tumor_leaves) + 1)]
    bin_info = bin_info.copy()
    new_cnps = []
    for _, row in bin_info.iterrows():
        key = tuple(
            (int(v.split("|")[0]), int(v.split("|")[1]))
            for v in (row[c] for c in clone_cols)
        )
        s = leaf_cn_to_seg[key]
        parts = [
            f"{int(tree.a_all[s, tree.normal_leaf])}|{int(tree.b_all[s, tree.normal_leaf])}"
        ]
        for v in cnp_rev:
            parts.append(f"{int(tree.a_all[s, v])}|{int(tree.b_all[s, v])}")
        new_cnps.append(";".join(parts))
    bin_info["CNP"] = new_cnps
    prop_parts = [str(np.mean(tree.node_props[tree.normal_leaf]))]
    for v in cnp_rev:
        prop_parts.append(str(np.mean(tree.node_props[v])))
    bin_info["PROPS"] = ";".join(prop_parts)

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

        _profile_fn = plot_ascn_profile if plot_ascn else plot_cnv_profile
        _profile_fn(
            ax_cn,
            bin_info,
            regions,
            width=20,
            height=1,
            show_prop=False,
            show_clone_name=False,
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
        for parent, (l, r) in tumor_children.items():
            px, py = x_weighted[parent], node_y[parent]
            for c in (l, r):
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

        _legend_fn = plot_ascn_legend if plot_ascn else plot_cnv_legend
        _legend_fn(ax_leg)
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
