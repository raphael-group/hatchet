"""Clone tree data structure and unlabeled full binary tree enumeration.

Node ordering convention:
    v1            = normal leaf (fixed CN=(1,1)), direct child of root
    v2..vn        = tumor leaves (latent mixture components)
    v_{n+1}..v_{2n-2} = internal tumor nodes
    v_{2n-1}      = root (fixed CN=(1,1))
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass
class CloneTree:
    """Rooted full binary clone tree (unlabeled shape)."""

    n: int  # total leaves (normal + tumor)
    n_nodes: int  # 2n-1
    root: int  # v_{2n-1}
    normal_leaf: int  # v1
    parent: dict[int, int] = field(repr=False)  # child -> parent
    children: dict[int, tuple[int, int]] = field(repr=False)  # node -> (left, right)
    leaves: list[int] = field(repr=False)  # [1, 2, ..., n]
    tumor_leaves: list[int] = field(repr=False)  # [2, ..., n]
    internal_nodes: list[int] = field(repr=False)  # [n+1, ..., 2n-1]
    edges: list[tuple[int, int]] = field(repr=False)  # [(parent, child), ...]

    @property
    def n_edges(self) -> int:
        return len(self.edges)

    @property
    def n_tumor_leaves(self) -> int:
        return self.n - 1

    @property
    def tumor_edges(self) -> list[tuple[int, int]]:
        return [(p, c) for p, c in self.edges if c != self.normal_leaf]

    def label(
        self, a_all: np.ndarray, b_all: np.ndarray, events: dict[str, np.ndarray]
    ) -> LabeledCloneTree:
        """Attach inferred CN and event values to produce a LabeledCloneTree."""
        return LabeledCloneTree(
            n=self.n,
            n_nodes=self.n_nodes,
            root=self.root,
            normal_leaf=self.normal_leaf,
            parent=self.parent,
            children=self.children,
            leaves=self.leaves,
            tumor_leaves=self.tumor_leaves,
            internal_nodes=self.internal_nodes,
            edges=self.edges,
            a_all=a_all,
            b_all=b_all,
            events=events,
        )


@dataclass
class LabeledCloneTree(CloneTree):
    """CloneTree with inferred CN at all nodes and event variables on edges.

    Attributes:
        a_all: (S, V+1) haplotype-A CN at every node (1-indexed columns).
        b_all: (S, V+1) haplotype-B CN at every node.
        events: dict of (S, n_tumor_edges) arrays with keys:
            alpha_a, alpha_b, delta_a, delta_b (copy counts),
            abar_a, abar_b, dbar_a, dbar_b (interval starts).
    """

    a_all: np.ndarray = field(default=None, repr=False)
    b_all: np.ndarray = field(default=None, repr=False)
    events: dict[str, np.ndarray] = field(default=None, repr=False)

    def _edge_cost(self, ei):
        """Total interval starts on tumor edge ei."""
        return int(
            self.events["abar_a"][:, ei].sum()
            + self.events["abar_b"][:, ei].sum()
            + self.events["dbar_a"][:, ei].sum()
            + self.events["dbar_b"][:, ei].sum()
        )

    def _cn_nhx(self, v):
        """NHX annotation string for node CN profiles."""
        a = "/".join(str(int(x)) for x in self.a_all[:, v])
        b = "/".join(str(int(x)) for x in self.b_all[:, v])
        return f"A={a}:B={b}"

    def to_newick(self) -> str:
        """Export as NHX Newick with node CN annotations and edge tree costs.

        Format: standard NHX — annotations in [&&NHX:key=val:key=val] after
        branch length. Compatible with ete3, FigTree, toytree.
        """
        te_map = {(p, c): i for i, (p, c) in enumerate(self.tumor_edges)}

        def _fmt(v):
            name = f"v{v}"
            nhx = f"[&&NHX:{self._cn_nhx(v)}]"
            if v in self.children:
                left, right = self.children[v]
                inner = f"({_fmt(left)},{_fmt(right)}){name}"
            else:
                inner = name
            if v == self.root:
                return f"{inner}{nhx}"
            edge_key = (self.parent[v], v)
            cost = self._edge_cost(te_map[edge_key]) if edge_key in te_map else 0
            return f"{inner}:{cost}{nhx}"

        return _fmt(self.root) + ";"

    def to_dict(self) -> dict:
        """Full solution as a JSON-serializable dict."""
        te = self.tumor_edges
        te_map = {(p, c): i for i, (p, c) in enumerate(te)}

        nodes = {}
        for v in range(1, self.n_nodes + 1):
            nodes[f"v{v}"] = {
                "a": self.a_all[:, v].astype(int).tolist(),
                "b": self.b_all[:, v].astype(int).tolist(),
            }

        edges = {}
        for p, c in te:
            ei = te_map[(p, c)]
            edges[f"v{p}->v{c}"] = {
                "cost": self._edge_cost(ei),
                "alpha_a": self.events["alpha_a"][:, ei].astype(int).tolist(),
                "alpha_b": self.events["alpha_b"][:, ei].astype(int).tolist(),
                "delta_a": self.events["delta_a"][:, ei].astype(int).tolist(),
                "delta_b": self.events["delta_b"][:, ei].astype(int).tolist(),
                "abar_a": self.events["abar_a"][:, ei].astype(int).tolist(),
                "abar_b": self.events["abar_b"][:, ei].astype(int).tolist(),
                "dbar_a": self.events["dbar_a"][:, ei].astype(int).tolist(),
                "dbar_b": self.events["dbar_b"][:, ei].astype(int).tolist(),
            }

        return {
            "n": self.n,
            "root": self.root,
            "normal_leaf": self.normal_leaf,
            "topology": [[p, c] for p, c in self.edges],
            "nodes": nodes,
            "edges": edges,
        }


# ---------------------------------------------------------------------------
# Unlabeled shape enumeration
# ---------------------------------------------------------------------------


def _enum_shapes(k):
    """Yield all unlabeled rooted full binary tree shapes on k leaves.

    Each shape is a nested tuple: a leaf is ``None``, an internal node is
    ``(left_shape, right_shape)`` with ``left_shape <= right_shape`` to
    canonicalise mirror images.
    """
    if k == 1:
        yield None
        return
    for left_k in range(1, k // 2 + 1):
        right_k = k - left_k
        left_shapes = list(_enum_shapes(left_k))
        if left_k == right_k:
            for i, ls in enumerate(left_shapes):
                right_shapes = list(_enum_shapes(right_k))
                for j in range(i, len(right_shapes)):
                    yield (ls, right_shapes[j])
        else:
            right_shapes = list(_enum_shapes(right_k))
            for ls in left_shapes:
                for rs in right_shapes:
                    yield (ls, rs)


def _shape_to_tree(shape, next_leaf, next_internal):
    """Convert a nested-tuple shape into a children dict.

    Returns (root_id, children_dict, next_leaf, next_internal).
    """
    if shape is None:
        node = next_leaf
        return node, {}, next_leaf + 1, next_internal

    left_shape, right_shape = shape
    lr, lc, nl, ni = _shape_to_tree(left_shape, next_leaf, next_internal)
    rr, rc, nl2, ni2 = _shape_to_tree(right_shape, nl, ni)
    node = ni2
    children = {}
    children.update(lc)
    children.update(rc)
    children[node] = (lr, rr)
    return node, children, nl2, ni2 + 1


def enumerate_binary_trees(n_leaves: int) -> list[CloneTree]:
    """Enumerate all unlabeled rooted full binary trees with n_leaves total leaves.

    v1 is the normal leaf (always direct child of root).  Only tumor-subtree
    shapes (on n_leaves-1 tumor leaves) are enumerated.
    """
    n = n_leaves
    if n < 2:
        raise ValueError(f"need >= 2 leaves, got {n}")

    k = n - 1  # tumor leaves
    results = []

    if k == 1:
        # trivial: root -> (v1, v2)
        root = 2 * n - 1  # v3
        children = {root: (1, 2)}
        parent = {1: root, 2: root}
        results.append(
            CloneTree(
                n=n,
                n_nodes=2 * n - 1,
                root=root,
                normal_leaf=1,
                parent=parent,
                children=children,
                leaves=[1, 2],
                tumor_leaves=[2],
                internal_nodes=[root],
                edges=[(root, 1), (root, 2)],
            )
        )
        return results

    for shape in _enum_shapes(k):
        # Build tumor subtree: leaves v2..vn, internals v_{n+1}..v_{2n-2}
        tumor_root, tumor_children, _, _ = _shape_to_tree(
            shape, next_leaf=2, next_internal=n + 1
        )

        # Renumber so tumor subtree root = v_{2n-2}
        target_tumor_root = 2 * n - 2
        if tumor_root != target_tumor_root:
            # Swap tumor_root and target in the children dict
            swapped = {}
            for p, (l, r) in tumor_children.items():
                pp = (
                    target_tumor_root
                    if p == tumor_root
                    else (tumor_root if p == target_tumor_root else p)
                )
                ll = (
                    target_tumor_root
                    if l == tumor_root
                    else (tumor_root if l == target_tumor_root else l)
                )
                rr = (
                    target_tumor_root
                    if r == tumor_root
                    else (tumor_root if r == target_tumor_root else r)
                )
                swapped[pp] = (ll, rr)
            tumor_children = swapped
            tumor_root = target_tumor_root

        # Full tree: root v_{2n-1} -> (v1, tumor_root)
        root = 2 * n - 1
        children = dict(tumor_children)
        children[root] = (1, tumor_root)

        # Build parent map and edge list
        parent = {}
        edges = []
        for p, (l, r) in children.items():
            parent[l] = p
            parent[r] = p
            edges.append((p, l))
            edges.append((p, r))

        leaves = list(range(1, n + 1))
        tumor_leaves = list(range(2, n + 1))
        internal_nodes = sorted(children.keys())

        results.append(
            CloneTree(
                n=n,
                n_nodes=2 * n - 1,
                root=root,
                normal_leaf=1,
                parent=parent,
                children=children,
                leaves=leaves,
                tumor_leaves=tumor_leaves,
                internal_nodes=internal_nodes,
                edges=edges,
            )
        )

    return results


def parse_newick(newick_str: str) -> CloneTree:
    """Parse a Newick string into a CloneTree.

    Expected format: a rooted full binary tree where the root has exactly
    two children — one leaf labeled ``normal`` (or ``N``) and a tumor
    subtree.  All other leaves become tumor clones v2..vn.

    Example for n=4::

        ((tumor1,(tumor2,tumor3)),normal);

    Raises ValueError if the tree is not full binary or if the root does
    not have exactly one normal leaf child.
    """
    s = newick_str.strip().rstrip(";").strip()

    # recursive descent parser
    _pos = [0]

    def _parse():
        if s[_pos[0]] == "(":
            _pos[0] += 1  # skip '('
            left = _parse()
            if s[_pos[0]] != ",":
                raise ValueError(f"expected ',' at pos {_pos[0]}")
            _pos[0] += 1  # skip ','
            right = _parse()
            if s[_pos[0]] != ")":
                raise ValueError(f"expected ')' at pos {_pos[0]}")
            _pos[0] += 1  # skip ')'
            # skip optional branch length / internal label
            while _pos[0] < len(s) and s[_pos[0]] not in (",", ")", ";"):
                _pos[0] += 1
            return (left, right)
        else:
            start = _pos[0]
            while _pos[0] < len(s) and s[_pos[0]] not in (",", ")", ";", ":"):
                _pos[0] += 1
            name = s[start : _pos[0]].strip()
            # skip branch length
            if _pos[0] < len(s) and s[_pos[0]] == ":":
                _pos[0] += 1
                while _pos[0] < len(s) and s[_pos[0]] not in (",", ")", ";"):
                    _pos[0] += 1
            return name

    tree = _parse()

    # tree must be a tuple (root is internal)
    if not isinstance(tree, tuple):
        raise ValueError("Newick tree must have an internal root, got a single leaf")

    # Identify normal leaf — must be a direct child of root
    left, right = tree
    if isinstance(left, str) and left.lower() in ("normal", "n"):
        tumor_subtree = right
    elif isinstance(right, str) and right.lower() in ("normal", "n"):
        tumor_subtree = left
    else:
        raise ValueError(
            "Root must have exactly one direct leaf child named 'normal' or 'N'. "
            f"Got children: {left}, {right}"
        )

    # Collect tumor leaves in order, verify full binary
    tumor_leaves_names = []

    def _collect(node):
        if isinstance(node, str):
            tumor_leaves_names.append(node)
        else:
            if len(node) != 2:
                raise ValueError(
                    f"Tree must be full binary, got node with {len(node)} children"
                )
            _collect(node[0])
            _collect(node[1])

    _collect(tumor_subtree)

    k = len(tumor_leaves_names)
    n = k + 1  # total leaves including normal

    # Build the tree with canonical node ids
    _next_internal = [n + 1]

    def _build(node):
        """Returns node_id and updates children dict."""
        if isinstance(node, str):
            idx = tumor_leaves_names.index(node) + 2  # v2..vn
            return idx
        left_id = _build(node[0])
        right_id = _build(node[1])
        nid = _next_internal[0]
        _next_internal[0] += 1
        _children[nid] = (left_id, right_id)
        return nid

    _children = {}
    tumor_root = _build(tumor_subtree)

    # Renumber so tumor subtree root = v_{2n-2}
    target = 2 * n - 2
    if tumor_root != target:
        swapped = {}
        for p, (l, r) in _children.items():
            pp = target if p == tumor_root else (tumor_root if p == target else p)
            ll = target if l == tumor_root else (tumor_root if l == target else l)
            rr = target if r == tumor_root else (tumor_root if r == target else r)
            swapped[pp] = (ll, rr)
        _children = swapped
        tumor_root = target

    # Full tree: root v_{2n-1} -> (v1, tumor_root)
    root = 2 * n - 1
    children = dict(_children)
    children[root] = (1, tumor_root)

    parent = {}
    edges = []
    for p, (l, r) in children.items():
        parent[l] = p
        parent[r] = p
        edges.append((p, l))
        edges.append((p, r))

    return CloneTree(
        n=n,
        n_nodes=2 * n - 1,
        root=root,
        normal_leaf=1,
        parent=parent,
        children=children,
        leaves=list(range(1, n + 1)),
        tumor_leaves=list(range(2, n + 1)),
        internal_nodes=sorted(children.keys()),
        edges=edges,
    )
