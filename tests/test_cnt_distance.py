"""Unit tests for CNT distance computation."""

import numpy as np
import pytest

from hatchet.compute_cn.solve.cnt_distance import (
    _cnt_distance_1d,
    compute_cnt_distances,
)


class TestCntDistance1d:
    """Test the single-allele, single-chromosome helper."""

    def test_zero_diff(self):
        assert _cnt_distance_1d(np.array([0, 0, 0])) == 0

    def test_uniform_amplification(self):
        # One contiguous amplification covering all 3 segments
        assert _cnt_distance_1d(np.array([1, 1, 1])) == 1

    def test_two_separate_amplifications(self):
        # Two separate amplification events (gap in the middle)
        assert _cnt_distance_1d(np.array([1, 0, 1])) == 2

    def test_nested_amplification(self):
        # Two nested amplification events: outer covers all 3, inner adds 1 more to first
        assert _cnt_distance_1d(np.array([2, 1, 0])) == 2

    def test_mixed_amp_and_del(self):
        # d = [-1, 2, -1]: d+=[0,2,0] → amp=2, d-=[1,0,1] → del=2, total=4
        assert _cnt_distance_1d(np.array([-1, 2, -1])) == 4

    def test_uniform_deletion(self):
        assert _cnt_distance_1d(np.array([-2, -2, -2])) == 2

    def test_single_segment(self):
        assert _cnt_distance_1d(np.array([3])) == 3
        assert _cnt_distance_1d(np.array([-2])) == 2
        assert _cnt_distance_1d(np.array([0])) == 0


class TestComputeCntDistances:
    """Test the full pairwise distance computation."""

    def test_identical_clones(self):
        cA = np.array([[1, 1], [1, 1]], dtype=int)
        cB = np.array([[1, 1], [1, 1]], dtype=int)
        boundaries = np.array([True, False])
        dist = compute_cnt_distances(cA, cB, boundaries)
        assert dist[0, 1] == 0
        assert dist[1, 0] == 0

    def test_symmetry(self):
        rng = np.random.default_rng(42)
        m, n = 20, 4
        cA = rng.integers(0, 4, size=(m, n))
        cB = rng.integers(0, 4, size=(m, n))
        boundaries = np.zeros(m, dtype=bool)
        boundaries[0] = True
        boundaries[10] = True
        dist = compute_cnt_distances(cA, cB, boundaries)
        np.testing.assert_array_equal(dist, dist.T)

    def test_diagonal_zero(self):
        rng = np.random.default_rng(7)
        m, n = 10, 3
        cA = rng.integers(0, 3, size=(m, n))
        cB = rng.integers(0, 3, size=(m, n))
        boundaries = np.array([True] + [False] * 9)
        dist = compute_cnt_distances(cA, cB, boundaries)
        np.testing.assert_array_equal(np.diag(dist), 0)

    def test_multi_chromosome(self):
        # Two chromosomes, two clones
        # Chr1 (2 segments): clone0 A=[1,1], clone1 A=[2,2] → d=[1,1] → 1 amp
        # Chr2 (2 segments): clone0 A=[1,1], clone1 A=[1,1] → d=[0,0] → 0
        # Same for B allele: all 1s both clones → 0
        cA = np.array([[1, 2], [1, 2], [1, 1], [1, 1]], dtype=int)
        cB = np.ones((4, 2), dtype=int)
        boundaries = np.array([True, False, True, False])
        dist = compute_cnt_distances(cA, cB, boundaries)
        assert dist[0, 1] == 1
        assert dist[1, 0] == 1

    def test_allele_specific(self):
        # Single chromosome, two clones
        # A: clone0=[1], clone1=[2]
        # B: clone0=[1], clone1=[0]
        # Reverse direction B: source=0, target=1 → infeasible → inf
        cA = np.array([[1, 2]], dtype=int)
        cB = np.array([[1, 0]], dtype=int)
        boundaries = np.array([True])
        dist = compute_cnt_distances(cA, cB, boundaries)
        assert dist[0, 1] == np.inf

    def test_allele_specific_feasible(self):
        # Both directions feasible (no zero-to-positive)
        # A: clone0=[1], clone1=[2] → 1 event
        # B: clone0=[1], clone1=[1] → 0 events
        cA = np.array([[1, 2]], dtype=int)
        cB = np.array([[1, 1]], dtype=int)
        boundaries = np.array([True])
        dist = compute_cnt_distances(cA, cB, boundaries)
        assert dist[0, 1] == 1
