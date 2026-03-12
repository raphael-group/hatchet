"""Integration tests for HATCHet3 compute-cn step."""

import os
import numpy as np
import pandas as pd
import pytest

from conftest import CBC_AVAILABLE

pytestmark = pytest.mark.skipif(not CBC_AVAILABLE, reason="CBC solver not available")


class TestComputeCnOutputFiles:
    """Verify compute-cn produces expected output files."""

    def test_best_bbc_ucn_exists(self, compute_cn_result):
        result_dir, _ = compute_cn_result
        assert os.path.isfile(os.path.join(result_dir, "best.bbc.ucn"))

    def test_best_seg_ucn_exists(self, compute_cn_result):
        result_dir, _ = compute_cn_result
        assert os.path.isfile(os.path.join(result_dir, "best.seg.ucn"))

    def test_gammas_exists(self, compute_cn_result):
        result_dir, _ = compute_cn_result
        assert os.path.isfile(os.path.join(result_dir, "gammas.tsv"))


class TestComputeCnOutputFormat:
    """Verify UCN file formats are correct."""

    def test_bbc_ucn_cn_columns(self, compute_cn_result):
        result_dir, _ = compute_cn_result
        bbc = pd.read_table(os.path.join(result_dir, "best.bbc.ucn"), sep="\t")
        assert "cn_normal" in bbc.columns, "Missing cn_normal column"
        assert "cn_clone1" in bbc.columns, "Missing cn_clone1 column"
        assert "u_normal" in bbc.columns, "Missing u_normal column"
        assert "u_clone1" in bbc.columns, "Missing u_clone1 column"

    def test_cn_format_pipe_separated(self, compute_cn_result):
        """CN values should be 'a|b' format."""
        result_dir, _ = compute_cn_result
        bbc = pd.read_table(os.path.join(result_dir, "best.bbc.ucn"), sep="\t")
        for cn_val in bbc["cn_clone1"].unique():
            parts = str(cn_val).split("|")
            assert len(parts) == 2, f"CN format should be 'a|b', got '{cn_val}'"
            assert parts[0].isdigit() and parts[1].isdigit(), f"CN parts should be integers: '{cn_val}'"

    def test_proportions_sum_to_one(self, compute_cn_result):
        """Clone proportions should sum to ~1.0 for each row."""
        result_dir, _ = compute_cn_result
        bbc = pd.read_table(os.path.join(result_dir, "best.bbc.ucn"), sep="\t")
        u_cols = [c for c in bbc.columns if c.startswith("u_")]
        prop_sums = bbc[u_cols].sum(axis=1)
        assert np.allclose(prop_sums, 1.0, atol=0.05), (
            f"Proportions should sum to ~1.0, got range [{prop_sums.min():.3f}, {prop_sums.max():.3f}]"
        )

    def test_cn_normal_is_diploid(self, compute_cn_result):
        """Normal clone should always be 1|1."""
        result_dir, _ = compute_cn_result
        bbc = pd.read_table(os.path.join(result_dir, "best.bbc.ucn"), sep="\t")
        assert (bbc["cn_normal"] == "1|1").all(), "Normal clone CN should be 1|1 everywhere"


class TestComputeCnRecovery:
    """Verify CN recovery against ground truth."""

    def test_balanced_region_cn(self, compute_cn_result):
        """Bins from balanced diploid regions should get cn_clone1='1|1'."""
        result_dir, ground_truth = compute_cn_result
        bbc = pd.read_table(os.path.join(result_dir, "best.bbc.ucn"), sep="\t")

        # chr1 p-arm bins (CN=(1,1) in ground truth) — check majority
        chr1p_bins = bbc[(bbc["#CHR"] == "chr1") & (bbc["END"] <= 121700000)]
        if len(chr1p_bins) > 0:
            frac_11 = (chr1p_bins["cn_clone1"] == "1|1").mean()
            assert frac_11 > 0.5, (
                f"Expected majority of chr1 p-arm bins to be 1|1, got {frac_11:.2f}"
            )

    def test_gamma_approximately_correct(self, compute_cn_result):
        """Gamma should be approximately 2.0 (since balanced cluster RDR ~ 1.0)."""
        result_dir, ground_truth = compute_cn_result
        gammas = pd.read_table(
            os.path.join(result_dir, "gammas.tsv"),
            sep="\t",
            header=None,
            names=["sample", "gamma_diploid", "gamma_tetraploid"],
        )
        gamma_dip = gammas["gamma_diploid"].iloc[0]
        assert 1.5 < gamma_dip < 3.0, (
            f"Expected gamma near 2.0, got {gamma_dip:.3f}"
        )
