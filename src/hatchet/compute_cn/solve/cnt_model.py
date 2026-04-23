"""CNT-CD model builders: C-step MILP (per chromosome) and U-step LP.

C-step: fix U, infer haplotype CN + CNT interval events. Two-stage lexicographic:
  1. C-step IMF stage: minimize fit loss (L1 or CI hinge)
  2. C-step CNT stage: minimize tree parsimony cost subject to near-optimal fit

U-step: fix leaf CN, optimize clone usages by LP.
"""

from __future__ import annotations

import logging

import numpy as np
from pyomo import environ as pe

from hatchet.compute_cn.solve.cnt_tree import CloneTree
from hatchet.compute_cn.solve.datatypes import SolverParams, SolverInputs


# ── C-step (per chromosome) ─────────────────────────────────────────────


def build_cnt_c_model(
    tree: CloneTree,
    params: SolverParams,
    inputs: SolverInputs,
    fixed_u: np.ndarray,
    seg_indices: np.ndarray,
) -> tuple[pe.ConcreteModel, dict]:
    """Build per-chromosome C-step MILP."""
    S = len(seg_indices)
    P = inputs.k
    cn_max = params.cn_max
    big_M = 2 * cn_max
    V = tree.n_nodes

    fa_obs = inputs.f_a.values[seg_indices]
    fb_obs = inputs.f_b.values[seg_indices]
    w_vals = inputs.w
    cluster_ids = inputs.cluster_ids

    model = pe.ConcreteModel()
    model.constraints = pe.ConstraintList()

    all_nodes = list(range(1, V + 1))
    fixed_nodes = {tree.normal_leaf, tree.root}
    free_nodes = [v for v in all_nodes if v not in fixed_nodes]

    # CN variables
    model.a = pe.Var(range(S), all_nodes, bounds=(0, cn_max), domain=pe.Integers)
    model.b = pe.Var(range(S), all_nodes, bounds=(0, cn_max), domain=pe.Integers)
    for s in range(S):
        model.a[s, tree.root].fix(1)
        model.b[s, tree.root].fix(1)
        model.a[s, tree.normal_leaf].fix(1)
        model.b[s, tree.normal_leaf].fix(1)

    # Nonzero indicators
    model.za = pe.Var(range(S), free_nodes, bounds=(0, 1), domain=pe.Binary)
    model.zb = pe.Var(range(S), free_nodes, bounds=(0, 1), domain=pe.Binary)
    for s in range(S):
        for v in free_nodes:
            model.constraints.add(model.a[s, v] <= cn_max * model.za[s, v])
            model.constraints.add(model.a[s, v] >= model.za[s, v])
            model.constraints.add(model.b[s, v] <= cn_max * model.zb[s, v])
            model.constraints.add(model.b[s, v] >= model.zb[s, v])

    # Leaf domain constraints
    seg_to_cluster = [cluster_ids[seg_indices[s]] for s in range(S)]
    w_total = sum(w_vals)
    balanced_set = set(inputs.balanced_clusters) if inputs.balanced_clusters else set()
    for s in range(S):
        cid = seg_to_cluster[s]
        for v in tree.tumor_leaves:
            model.constraints.add(model.a[s, v] + model.b[s, v] <= cn_max)
            if w_vals[cid] / w_total >= params.zero_cn_thres:
                model.constraints.add(model.a[s, v] + model.b[s, v] >= 1)
            if cid in balanced_set:
                model.constraints.add(model.a[s, v] == model.b[s, v])

    # Event variables on tumor edges
    tumor_edges = [(p, c) for p, c in tree.edges if c != tree.normal_leaf]
    te_idx = list(range(len(tumor_edges)))
    model.alpha_a = pe.Var(range(S), te_idx, bounds=(0, cn_max), domain=pe.Integers)
    model.alpha_b = pe.Var(range(S), te_idx, bounds=(0, cn_max), domain=pe.Integers)
    model.delta_a = pe.Var(range(S), te_idx, bounds=(0, cn_max), domain=pe.Integers)
    model.delta_b = pe.Var(range(S), te_idx, bounds=(0, cn_max), domain=pe.Integers)
    model.abar_a = pe.Var(range(S), te_idx, bounds=(0, cn_max), domain=pe.Integers)
    model.abar_b = pe.Var(range(S), te_idx, bounds=(0, cn_max), domain=pe.Integers)
    model.dbar_a = pe.Var(range(S), te_idx, bounds=(0, cn_max), domain=pe.Integers)
    model.dbar_b = pe.Var(range(S), te_idx, bounds=(0, cn_max), domain=pe.Integers)

    def _z(hap, s, v):
        return (
            1
            if v in fixed_nodes
            else (model.za[s, v] if hap == "a" else model.zb[s, v])
        )

    def _c(hap, s, v):
        return model.a[s, v] if hap == "a" else model.b[s, v]

    # CNT edge feasibility
    for ei, (u_node, v_node) in enumerate(tumor_edges):
        for hap, alpha_var, delta_var in [
            ("a", model.alpha_a, model.delta_a),
            ("b", model.alpha_b, model.delta_b),
        ]:
            for s in range(S):
                zu, zv = _z(hap, s, u_node), _z(hap, s, v_node)
                cu, cv = _c(hap, s, u_node), _c(hap, s, v_node)
                alpha, delta = alpha_var[s, ei], delta_var[s, ei]
                model.constraints.add(zv <= zu)
                model.constraints.add(cv + delta - cu - alpha <= big_M * (2 - zu - zv))
                model.constraints.add(cu + alpha - cv - delta <= big_M * (2 - zu - zv))
                model.constraints.add(delta <= cu - 1 + (cn_max + 1) * (2 - zu - zv))
                model.constraints.add(delta >= cu - big_M * (1 - zu + zv))

    # Interval-start counting
    for ei in te_idx:
        for abar_v, alpha_v, dbar_v, delta_v in [
            (model.abar_a, model.alpha_a, model.dbar_a, model.delta_a),
            (model.abar_b, model.alpha_b, model.dbar_b, model.delta_b),
        ]:
            for s in range(S):
                if s == 0:
                    model.constraints.add(abar_v[s, ei] >= alpha_v[s, ei])
                    model.constraints.add(dbar_v[s, ei] >= delta_v[s, ei])
                else:
                    model.constraints.add(
                        abar_v[s, ei] >= alpha_v[s, ei] - alpha_v[s - 1, ei]
                    )
                    model.constraints.add(
                        dbar_v[s, ei] >= delta_v[s, ei] - delta_v[s - 1, ei]
                    )

    # Mixture prediction (leaves only)
    model.fA = pe.Var(range(S), range(P), bounds=(0, cn_max), domain=pe.Reals)
    model.fB = pe.Var(range(S), range(P), bounds=(0, cn_max), domain=pe.Reals)
    leaf_to_col = {v: v - 1 for v in tree.leaves}
    for s in range(S):
        for p in range(P):
            sA = sum(
                model.a[s, v] * float(fixed_u[leaf_to_col[v], p]) for v in tree.leaves
            )
            sB = sum(
                model.b[s, v] * float(fixed_u[leaf_to_col[v], p]) for v in tree.leaves
            )
            model.constraints.add(model.fA[s, p] == sA)
            model.constraints.add(model.fB[s, p] == sB)

    # Fit loss
    if params.obj_type == "imf":
        model.yA = pe.Var(range(S), range(P), bounds=(0, None), domain=pe.Reals)
        model.yB = pe.Var(range(S), range(P), bounds=(0, None), domain=pe.Reals)
        for s in range(S):
            for p in range(P):
                model.constraints.add(
                    model.yA[s, p] >= float(fa_obs[s, p]) - model.fA[s, p]
                )
                model.constraints.add(
                    model.yA[s, p] >= model.fA[s, p] - float(fa_obs[s, p])
                )
                model.constraints.add(
                    model.yB[s, p] >= float(fb_obs[s, p]) - model.fB[s, p]
                )
                model.constraints.add(
                    model.yB[s, p] >= model.fB[s, p] - float(fb_obs[s, p])
                )
    else:
        model.hA = pe.Var(range(S), range(P), bounds=(0, None), domain=pe.Reals)
        model.hB = pe.Var(range(S), range(P), bounds=(0, None), domain=pe.Reals)
        fa_lo, fa_hi = (
            inputs.fa_lo.values[seg_indices],
            inputs.fa_hi.values[seg_indices],
        )
        fb_lo, fb_hi = (
            inputs.fb_lo.values[seg_indices],
            inputs.fb_hi.values[seg_indices],
        )
        for s in range(S):
            for p in range(P):
                model.constraints.add(
                    model.hA[s, p] >= model.fA[s, p] - float(fa_hi[s, p])
                )
                model.constraints.add(
                    model.hA[s, p] >= float(fa_lo[s, p]) - model.fA[s, p]
                )
                model.constraints.add(
                    model.hB[s, p] >= model.fB[s, p] - float(fb_hi[s, p])
                )
                model.constraints.add(
                    model.hB[s, p] >= float(fb_lo[s, p]) - model.fB[s, p]
                )

    # Objective expressions
    def _w(s):
        return float(w_vals[cluster_ids[seg_indices[s]]])

    if params.obj_type == "imf":
        fit_expr = sum(
            _w(s) * (model.yA[s, p] + model.yB[s, p])
            for s in range(S)
            for p in range(P)
        )
    else:
        fit_expr = sum(
            _w(s) * (model.hA[s, p] + model.hB[s, p])
            for s in range(S)
            for p in range(P)
        )

    tree_cost_expr = sum(
        model.abar_a[s, ei]
        + model.abar_b[s, ei]
        + model.dbar_a[s, ei]
        + model.dbar_b[s, ei]
        for s in range(S)
        for ei in te_idx
    )

    return model, {
        "fit_expr": fit_expr,
        "tree_cost_expr": tree_cost_expr,
        "seg_indices": seg_indices,
        "tumor_edges": tumor_edges,
        "te_idx": te_idx,
    }


def solve_cnt_c_lexi(model, aux, solver, eps_fit, eta=1e-6, timelimit=None):
    """Solve C-step MILP: IMF stage then CNT stage.

    Returns dict with F_star, F_actual, T_star, or None if infeasible.
    """
    fit_expr, tree_cost_expr = aux["fit_expr"], aux["tree_cost_expr"]
    opts = {"TimeLimit": timelimit} if timelimit else {}

    # C-step IMF stage
    model.obj = pe.Objective(expr=fit_expr, sense=pe.minimize)
    result = solver.solve(model, tee=False, options=opts)
    status = result.solver.termination_condition
    if status not in (
        pe.TerminationCondition.optimal,
        pe.TerminationCondition.maxTimeLimit,
    ):
        logging.warning(f"C-step IMF stage infeasible: {status}")
        return None
    if status == pe.TerminationCondition.maxTimeLimit and model.obj() is None:
        return None

    F_star = pe.value(model.obj)

    # C-step CNT stage
    eta_val = eta * max(1.0, abs(F_star))
    model.fit_bound = pe.Constraint(expr=fit_expr <= F_star + eps_fit + eta_val)
    model.del_component(model.obj)
    model.obj = pe.Objective(expr=tree_cost_expr, sense=pe.minimize)

    result = solver.solve(model, tee=False, warmstart=True, options=opts)
    status = result.solver.termination_condition
    if status not in (
        pe.TerminationCondition.optimal,
        pe.TerminationCondition.maxTimeLimit,
    ):
        logging.warning(
            f"C-step CNT stage infeasible ({status}), retrying with relaxed eta"
        )
        model.fit_bound.deactivate()
        model.fit_bound_relaxed = pe.Constraint(
            expr=fit_expr <= F_star + 10 * eta_val + eps_fit
        )
        result = solver.solve(model, tee=False, warmstart=True, options=opts)
        if result.solver.termination_condition not in (
            pe.TerminationCondition.optimal,
            pe.TerminationCondition.maxTimeLimit,
        ):
            logging.warning("C-step CNT stage still infeasible after relaxation")
            return None

    return {
        "F_star": F_star,
        "F_actual": pe.value(fit_expr),
        "T_star": pe.value(tree_cost_expr),
    }


def extract_cnt_c(model, tree, aux):
    """Extract CN and event variables from solved C-step model."""
    S = len(aux["seg_indices"])
    V = tree.n_nodes
    all_nodes = list(range(1, V + 1))
    te_idx = aux["te_idx"]
    n_te = len(te_idx)

    a = np.zeros((S, V + 1))
    b = np.zeros((S, V + 1))
    for s in range(S):
        for v in all_nodes:
            a[s, v] = round(pe.value(model.a[s, v]))
            b[s, v] = round(pe.value(model.b[s, v]))

    result = {"a": a, "b": b}
    for name, var in [
        ("alpha_a", model.alpha_a),
        ("alpha_b", model.alpha_b),
        ("delta_a", model.delta_a),
        ("delta_b", model.delta_b),
        ("abar_a", model.abar_a),
        ("abar_b", model.abar_b),
        ("dbar_a", model.dbar_a),
        ("dbar_b", model.dbar_b),
    ]:
        arr = np.zeros((S, n_te))
        for s in range(S):
            for ei in te_idx:
                arr[s, ei] = round(pe.value(var[s, ei]))
        result[name] = arr
    return result


# ── U-step (global) ──────────────────────────────────────────────────────


def build_cnt_u_model(params, inputs, fixed_a, fixed_b):
    """Build U-step LP for fixed leaf copy numbers.

    Args:
        fixed_a: (S, n) leaf haplotype-A CN.
        fixed_b: (S, n) leaf haplotype-B CN.
    """
    S, P, n = inputs.m, inputs.k, params.n
    fa_obs, fb_obs = inputs.f_a.values, inputs.f_b.values

    model = pe.ConcreteModel()
    model.constraints = pe.ConstraintList()

    model.u = pe.Var(range(n), range(P), bounds=(0, 1), domain=pe.Reals)
    for p in range(P):
        model.constraints.add(sum(model.u[i, p] for i in range(n)) == 1)

    if inputs.purities is not None:
        for p, sid in enumerate(inputs.sample_ids):
            if sid in inputs.purities:
                model.u[0, p].fix(1 - inputs.purities[sid])

    if params.minprop > 0:
        model.x = pe.Var(range(1, n), range(P), domain=pe.Binary)
        for i in range(1, n):
            for p in range(P):
                model.constraints.add(model.x[i, p] >= model.u[i, p])
                model.constraints.add(model.u[i, p] >= params.minprop * model.x[i, p])

    model.fA = pe.Var(range(S), range(P), bounds=(0, params.cn_max), domain=pe.Reals)
    model.fB = pe.Var(range(S), range(P), bounds=(0, params.cn_max), domain=pe.Reals)
    for s in range(S):
        for p in range(P):
            model.constraints.add(
                model.fA[s, p]
                == sum(float(fixed_a[s, i]) * model.u[i, p] for i in range(n))
            )
            model.constraints.add(
                model.fB[s, p]
                == sum(float(fixed_b[s, i]) * model.u[i, p] for i in range(n))
            )

    if params.obj_type == "imf":
        model.yA = pe.Var(range(S), range(P), bounds=(0, None), domain=pe.Reals)
        model.yB = pe.Var(range(S), range(P), bounds=(0, None), domain=pe.Reals)
        for s in range(S):
            for p in range(P):
                model.constraints.add(
                    model.yA[s, p] >= float(fa_obs[s, p]) - model.fA[s, p]
                )
                model.constraints.add(
                    model.yA[s, p] >= model.fA[s, p] - float(fa_obs[s, p])
                )
                model.constraints.add(
                    model.yB[s, p] >= float(fb_obs[s, p]) - model.fB[s, p]
                )
                model.constraints.add(
                    model.yB[s, p] >= model.fB[s, p] - float(fb_obs[s, p])
                )
    else:
        model.hA = pe.Var(range(S), range(P), bounds=(0, None), domain=pe.Reals)
        model.hB = pe.Var(range(S), range(P), bounds=(0, None), domain=pe.Reals)
        fa_lo, fa_hi = inputs.fa_lo.values, inputs.fa_hi.values
        fb_lo, fb_hi = inputs.fb_lo.values, inputs.fb_hi.values
        for s in range(S):
            for p in range(P):
                model.constraints.add(
                    model.hA[s, p] >= model.fA[s, p] - float(fa_hi[s, p])
                )
                model.constraints.add(
                    model.hA[s, p] >= float(fa_lo[s, p]) - model.fA[s, p]
                )
                model.constraints.add(
                    model.hB[s, p] >= model.fB[s, p] - float(fb_hi[s, p])
                )
                model.constraints.add(
                    model.hB[s, p] >= float(fb_lo[s, p]) - model.fB[s, p]
                )

    cids = inputs.cluster_ids
    if params.obj_type == "imf":
        fit_expr = sum(
            float(inputs.w[cids[s]]) * (model.yA[s, p] + model.yB[s, p])
            for s in range(S)
            for p in range(P)
        )
    else:
        fit_expr = sum(
            float(inputs.w[cids[s]]) * (model.hA[s, p] + model.hB[s, p])
            for s in range(S)
            for p in range(P)
        )

    model.obj = pe.Objective(expr=fit_expr, sense=pe.minimize)
    return model, {"fit_expr": fit_expr}


def extract_u(model, params, inputs):
    """Extract U matrix (n x P) from solved U-step model."""
    u = np.zeros((params.n, inputs.k))
    for i in range(params.n):
        for p in range(inputs.k):
            u[i, p] = pe.value(model.u[i, p])
    return u
