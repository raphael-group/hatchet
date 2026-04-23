"""Unit tests for CNT distance computation."""

import numpy as np

from hatchet.compute_cn.solve.cnt_distance import (
    _cnt_distance_1d,
    compute_cnt_distances,
)


class TestCntDistance1d:
    """Test the single-allele, single-chromosome helper."""

    def test_zero_diff(self):
        assert _cnt_distance_1d(np.array([0, 0, 0])) == 0

    def test_uniform_amplification(self):
        assert _cnt_distance_1d(np.array([2, 2, 2])) == 2

    def test_two_separate_amplifications(self):
        assert _cnt_distance_1d(np.array([1, 0, 1])) == 2

    def test_nested_amplification(self):
        assert _cnt_distance_1d(np.array([1, 2, 1])) == 2

    def test_mixed_amp_and_del(self):
        d = np.array([-1, 2, -1])
        assert _cnt_distance_1d(d) == 4

    def test_uniform_deletion(self):
        assert _cnt_distance_1d(np.array([-3, -3])) == 3

    def test_single_segment(self):
        assert _cnt_distance_1d(np.array([5])) == 5
        assert _cnt_distance_1d(np.array([-3])) == 3


class TestComputeCntDistances:
    """Test the pairwise distance matrix builder."""

    def test_identical_clones(self):
        cA = np.array([[1, 1], [2, 2]], dtype=int)
        cB = np.array([[1, 1], [1, 1]], dtype=int)
        boundaries = np.array([True, False])
        dist = compute_cnt_distances(cA, cB, boundaries)
        assert dist[0, 1] == 0
        assert dist[1, 0] == 0

    def test_asymmetry(self):
        # CND is asymmetric when LOH is present
        # Clone0: A=[1], B=[1]; Clone1: A=[2], B=[0]
        # 0→1: A 1→2 (1 amp), B 1→0 (1 del) → 2
        # 1→0: A 2→1 (1 del), B 0→1 (infeasible) → inf
        cA = np.array([[1, 2]], dtype=int)
        cB = np.array([[1, 0]], dtype=int)
        boundaries = np.array([True])
        dist = compute_cnt_distances(cA, cB, boundaries)
        assert dist[0, 1] == 2
        assert dist[1, 0] == np.inf

    def test_symmetric_when_no_loh(self):
        # Without zeros, CND is symmetric (same events forward and backward)
        cA = np.array([[1, 2], [2, 1]], dtype=int)
        cB = np.array([[1, 1], [1, 1]], dtype=int)
        boundaries = np.array([True, False])
        dist = compute_cnt_distances(cA, cB, boundaries)
        assert dist[0, 1] == dist[1, 0]

    def test_diagonal_zero(self):
        rng = np.random.default_rng(7)
        m, n = 10, 3
        cA = rng.integers(0, 3, size=(m, n))
        cB = rng.integers(0, 3, size=(m, n))
        boundaries = np.array([True] + [False] * 9)
        dist = compute_cnt_distances(cA, cB, boundaries)
        np.testing.assert_array_equal(np.diag(dist), 0)

    def test_multi_chromosome(self):
        # Two chromosomes, two clones, no LOH
        # Chr1 (2 segments): clone0 A=[1,1], clone1 A=[2,2] → 1 amp
        # Chr2 (2 segments): clone0 A=[1,1], clone1 A=[1,1] → 0
        # B allele: all 1s → 0
        cA = np.array([[1, 2], [1, 2], [1, 1], [1, 1]], dtype=int)
        cB = np.ones((4, 2), dtype=int)
        boundaries = np.array([True, False, True, False])
        dist = compute_cnt_distances(cA, cB, boundaries)
        assert dist[0, 1] == 1
        assert dist[1, 0] == 1

    def test_allele_specific_feasible(self):
        # Both directions feasible (no zero-to-positive)
        # A: clone0=[1], clone1=[2] → 1 event
        # B: clone0=[1], clone1=[1] → 0 events
        cA = np.array([[1, 2]], dtype=int)
        cB = np.array([[1, 1]], dtype=int)
        boundaries = np.array([True])
        dist = compute_cnt_distances(cA, cB, boundaries)
        assert dist[0, 1] == 1
