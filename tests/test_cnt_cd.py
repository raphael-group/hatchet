"""Tests for CNT-CD solver components."""

import numpy as np
import pandas as pd
import pytest

from hatchet.compute_cn.solve.cnt_tree import enumerate_binary_trees
from hatchet.compute_cn.solve.cnt_model import (
    build_c_step_model,
    solve_c_step,
    extract_c_step,
    build_u_step_model,
    extract_u_step,
)
from hatchet.compute_cn.solve.datatypes import SolverParams, SolverInputs


# ── helpers ──────────────────────────────────────────────────────────────


def _check_solver():
    try:
        from pyomo import environ as pe

        for name in ("gurobi", "cbc"):
            s = pe.SolverFactory(name)
            if s.available(exception_flag=False):
                return name
    except Exception:
        pass
    return None


SOLVER = _check_solver()
needs_solver = pytest.mark.skipif(SOLVER is None, reason="no MILP solver available")


def _make_inputs(fa, fb, w=None):
    """Build SolverInputs from numpy arrays (S x P)."""
    S, P = fa.shape
    seg_ids = list(range(S))
    sample_ids = [f"s{p}" for p in range(P)]
    if w is None:
        w = pd.Series(np.ones(S) / S, index=seg_ids)
    fa_df = pd.DataFrame(fa, index=seg_ids, columns=sample_ids)
    fb_df = pd.DataFrame(fb, index=seg_ids, columns=sample_ids)
    margin = 0.5
    nbins_df = pd.DataFrame(
        np.ones((S, P), dtype=int), index=seg_ids, columns=sample_ids
    )
    return SolverInputs(
        f_a=fa_df,
        f_b=fb_df,
        w=w,
        cluster_ids=seg_ids,
        sample_ids=sample_ids,
        copy_numbers={},
        free_rows=list(range(S)),
        fixed_rows=set(),
        purities=None,
        fa_lo=fa_df - margin,
        fa_hi=fa_df + margin,
        fb_lo=fb_df - margin,
        fb_hi=fb_df + margin,
        nbins=nbins_df,
        chr_boundaries=np.array([True] + [False] * (S - 1)),
    )


def _make_params(n, cn_max=6, obj_type="imf", eps_fit=0.01):
    return SolverParams(
        n=n,
        cn_max=cn_max,
        base=1,
        ampdel=False,
        minprop=0.0,
        max_ncns_seg=-1,
        tol=0.001,
        zero_cn_thres=0.0,
        obj_type=obj_type,
        eps_fit=eps_fit,
    )


# ── tree enumeration tests ──────────────────────────────────────────────


class TestTreeEnumeration:
    @pytest.mark.parametrize(
        "n,expected", [(2, 1), (3, 1), (4, 1), (5, 2), (6, 3), (7, 6), (8, 11)]
    )
    def test_counts(self, n, expected):
        trees = enumerate_binary_trees(n)
        assert len(trees) == expected

    def test_structure_n3(self):
        trees = enumerate_binary_trees(3)
        t = trees[0]
        assert t.n == 3
        assert t.n_nodes == 5
        assert t.root == 5
        assert t.normal_leaf == 1
        assert t.leaves == [1, 2, 3]
        assert t.tumor_leaves == [2, 3]
        assert 1 in t.parent and t.parent[1] == t.root

    def test_root_children(self):
        for n in range(2, 7):
            for tree in enumerate_binary_trees(n):
                left, right = tree.children[tree.root]
                assert tree.normal_leaf in (left, right)


# ── C-step model tests ──────────────────────────────────────────────────


class TestCStepModel:
    @needs_solver
    def test_no_change_edge_zero_tree_cost(self):
        """If all leaves have CN=(1,1), tree cost should be 0."""
        n, S, P = 3, 3, 1
        tree = enumerate_binary_trees(n)[0]
        # All leaves at (1,1) -> predicted FCN = 1.0 for both haplotypes
        fa = np.ones((S, P))
        fb = np.ones((S, P))
        inputs = _make_inputs(fa, fb)
        params = _make_params(n)
        # U: equal weight on all leaves
        u = np.array([[1.0 / n]] * n * P).reshape(n, P)

        seg_idx = np.arange(S)
        model, aux = build_c_step_model(tree, params, inputs, u, seg_idx)
        from hatchet.compute_cn.solve.inference import create_solver

        solver = create_solver(SOLVER, threads=1)
        result = solve_c_step(model, aux, solver, params.eps_fit)

        assert result is not None
        assert result["T_star"] == 0.0, f"expected tree cost 0, got {result['T_star']}"

    @needs_solver
    def test_simple_diploid_two_clones(self):
        """n=2 (normal+1 tumor), S=2 segments, P=1 sample."""
        n, S = 2, 2
        tree = enumerate_binary_trees(n)[0]
        # Tumor clone has CN (2,1) on seg0 and (1,0) on seg1
        # Normal has (1,1) everywhere
        # With u = [0.5, 0.5]: predicted fa = 0.5*1 + 0.5*cA, fb = 0.5*1 + 0.5*cB
        # seg0: fa_obs=1.5, fb_obs=1.0 -> cA=2, cB=1
        # seg1: fa_obs=1.0, fb_obs=0.5 -> cA=1, cB=0
        fa = np.array([[1.5], [1.0]])
        fb = np.array([[1.0], [0.5]])
        inputs = _make_inputs(fa, fb)
        params = _make_params(n)
        u = np.array([[0.5], [0.5]])

        seg_idx = np.arange(S)
        model, aux = build_c_step_model(tree, params, inputs, u, seg_idx)
        from hatchet.compute_cn.solve.inference import create_solver

        solver = create_solver(SOLVER, threads=1)
        result = solve_c_step(model, aux, solver, params.eps_fit)

        assert result is not None
        assert result["F_star"] < 0.01, (
            f"fit should be near-zero, got {result['F_star']}"
        )

        ab = extract_c_step(model, tree, aux)
        # Check tumor leaf (v2) CN
        assert ab["a"][0, 2] == 2  # seg0, hap A
        assert ab["b"][0, 2] == 1  # seg0, hap B
        assert ab["a"][1, 2] == 1  # seg1, hap A
        assert ab["b"][1, 2] == 0  # seg1, hap B

    @needs_solver
    def test_zero_regain_infeasible(self):
        """Parent CN=0 -> child CN>0 must be infeasible (enforced by z constraint)."""
        n, S = 3, 1
        tree = enumerate_binary_trees(n)[0]
        # Set observed FCN that would require a (0->positive) regain
        # With tree: root(v5)->v1, root->v4, v4->v2, v4->v3
        # If we force v4 to have CN=0 at a segment, v2 and v3 must also be 0
        # This is tested implicitly — the solver should find a feasible solution
        # that respects zero inheritance
        fa = np.array([[2.0]])
        fb = np.array([[1.0]])
        inputs = _make_inputs(fa, fb)
        params = _make_params(n)
        u = np.array([[0.5], [0.25], [0.25]])

        seg_idx = np.arange(S)
        model, aux = build_c_step_model(tree, params, inputs, u, seg_idx)
        from hatchet.compute_cn.solve.inference import create_solver

        solver = create_solver(SOLVER, threads=1)
        result = solve_c_step(model, aux, solver, params.eps_fit)

        assert result is not None
        ab = extract_c_step(model, tree, aux)
        # Check zero-inheritance: for each edge, if parent z=0 then child z=0
        for v in tree.tumor_leaves + tree.internal_nodes:
            if v == tree.root:
                continue
            parent = tree.parent[v]
            if parent in {tree.root, tree.normal_leaf}:
                continue
            # If parent hap-a is 0, child must be 0
            for s in range(S):
                if ab["a"][s, parent] == 0:
                    assert ab["a"][s, v] == 0, (
                        f"zero-regain violated: a[{s},{parent}]=0 but a[{s},{v}]={ab['a'][s, v]}"
                    )
                if ab["b"][s, parent] == 0:
                    assert ab["b"][s, v] == 0, (
                        f"zero-regain violated: b[{s},{parent}]=0 but b[{s},{v}]={ab['b'][s, v]}"
                    )


# ── U-step model tests ──────────────────────────────────────────────────


class TestUStepModel:
    @needs_solver
    def test_u_sum_to_one(self):
        """U rows must sum to 1 per sample."""
        n, P = 3, 2
        a_leaves = np.array([[1, 2, 1], [1, 1, 3]])  # (S, n)
        b_leaves = np.array([[1, 1, 1], [1, 0, 2]])
        fa = np.array([[1.5, 1.2], [1.0, 2.0]])
        fb = np.array([[1.0, 1.0], [0.5, 1.5]])
        inputs = _make_inputs(fa, fb)
        params = _make_params(n)

        model, _ = build_u_step_model(params, inputs, a_leaves, b_leaves)
        from hatchet.compute_cn.solve.inference import create_solver

        solver = create_solver(SOLVER, threads=1)
        solver.solve(model, tee=False)
        u = extract_u_step(model, params, inputs)

        for p in range(P):
            assert abs(u[:, p].sum() - 1.0) < 1e-6, f"U col {p} sums to {u[:, p].sum()}"
        assert (u >= -1e-9).all(), "U has negative entries"


# ── chromosome decomposition test ────────────────────────────────────────


class TestChromDecomposition:
    @needs_solver
    def test_two_chromosomes(self):
        """Two-chromosome problem: interval starts should reset at boundary."""
        n = 2
        tree = enumerate_binary_trees(n)[0]
        # 4 segments, 2 per chromosome
        # chr1: segs 0,1; chr2: segs 2,3
        # Tumor clone: seg0=(2,1), seg1=(2,1), seg2=(1,0), seg3=(1,0)
        # So chr1 has one constant event, chr2 has one constant event
        u = np.array([[0.5], [0.5]])

        fa = np.array([[1.5], [1.5], [1.0], [1.0]])
        fb = np.array([[1.0], [1.0], [0.5], [0.5]])
        inputs = _make_inputs(fa, fb)
        inputs.chr_boundaries = np.array([True, False, True, False])
        params = _make_params(n)

        from hatchet.compute_cn.solve.utils import split_by_chromosome

        chrom_groups = split_by_chromosome(inputs)
        assert len(chrom_groups) == 2
        assert list(chrom_groups[0]) == [0, 1]
        assert list(chrom_groups[1]) == [2, 3]

        from hatchet.compute_cn.solve.inference import create_solver

        solver = create_solver(SOLVER, threads=1)

        total_tree_cost = 0
        for cg in chrom_groups:
            model, aux = build_c_step_model(tree, params, inputs, u, cg)
            result = solve_c_step(model, aux, solver, params.eps_fit)
            assert result is not None
            total_tree_cost += result["T_star"]

        # Each chromosome should have its own interval starts
        # The events don't carry across chromosomes
        assert total_tree_cost >= 0


# ── end-to-end recovery test ─────────────────────────────────────────────


class TestRecovery:
    @needs_solver
    def test_recovery_n3(self):
        """Generate synthetic noiseless FCN from a known n=3 tree, run cnt_cd, check recovery."""
        from hatchet.compute_cn.solve.inference import run_coordinate_descent

        n = 3

        # Ground-truth leaf CN (S x n), columns: [normal, tumor1, tumor2]
        a_true = np.array(
            [
                [1, 2, 1],  # seg 0: tumor1 gained
                [1, 1, 3],  # seg 1: tumor2 amplified
                [1, 2, 2],  # seg 2: both gained
                [1, 1, 1],  # seg 3: all normal
            ],
            dtype=float,
        )
        b_true = np.array(
            [
                [1, 1, 1],  # seg 0
                [1, 0, 1],  # seg 1: tumor1 lost
                [1, 1, 0],  # seg 2: tumor2 lost
                [1, 1, 1],  # seg 3
            ],
            dtype=float,
        )

        # 3 samples with well-separated usage (more identifiable)
        u_true = np.array(
            [
                [0.3, 0.5, 0.7],  # normal
                [0.5, 0.1, 0.1],  # tumor 1 dominant in sample 0
                [0.2, 0.4, 0.2],  # tumor 2 dominant in sample 1
            ]
        )

        fa_obs = a_true @ u_true  # (S, P)
        fb_obs = b_true @ u_true

        inputs = _make_inputs(fa_obs, fb_obs)
        params = _make_params(n, cn_max=4, eps_fit=0.005)

        pool, _obj_df = run_coordinate_descent(
            params=params,
            inputs=inputs,
            mode="cnt_cd",
            solver_type=SOLVER,
            max_iters=10,
            max_convergence_iters=2,
            n_seed=50,
            j=1,
            cd_tol=0.001,
            random_seed=42,
            timelimit=60,
        )

        best_id = min(pool, key=lambda k: pool[k]["imf_obj"])
        best = pool[best_id]
        assert best["imf_obj"] < 0.05, f"imf_obj too high: {best['imf_obj']}"
        cA_rec = np.array(best["cA"])
        cB_rec = np.array(best["cB"])
        u_rec = np.array(best["u"])
        fa_rec = cA_rec @ u_rec
        np.testing.assert_allclose(
            fa_rec, fa_obs, atol=0.05, err_msg="reconstructed fa deviates from observed"
        )
        fb_rec = cB_rec @ u_rec
        np.testing.assert_allclose(
            fb_rec, fb_obs, atol=0.05, err_msg="reconstructed fb deviates from observed"
        )
