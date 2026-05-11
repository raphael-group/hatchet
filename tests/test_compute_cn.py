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
            assert parts[0].isdigit() and parts[1].isdigit(), (
                f"CN parts should be integers: '{cn_val}'"
            )

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
        assert (bbc["cn_normal"] == "1|1").all(), (
            "Normal clone CN should be 1|1 everywhere"
        )


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
        assert 1.5 < gamma_dip < 3.0, f"Expected gamma near 2.0, got {gamma_dip:.3f}"


class TestCDMode:
    """Test coordinate descent solver with different regularization terms."""

    @pytest.fixture(scope="class")
    def cd_base_args(self, cluster_bins_result, synthetic_data, tmp_path_factory):
        if not CBC_AVAILABLE:
            pytest.skip("CBC solver not available")
        bbc_dir, _ = cluster_bins_result
        _, genome_sizes, regions_bed, _ = synthetic_data
        # Only override defaults that differ from src/hatchet/hatchet.yaml.
        return {
            "bbc": os.path.join(bbc_dir, "bulk.bbc"),
            "seg": os.path.join(bbc_dir, "bulk.seg"),
            "genome_size": genome_sizes,
            "region_bed": regions_bed,
            "mode": "cd",
            "solver": "cbc",
            "timelimit": 30,
            "maxClone": 2,
            "diploid": True,
            "force": True,
            "reg_steps": 2,
            "reg_bound": 0.1,
            "diploidcmax": 6,
            "cd_niters": 5,
            "cd_nseeds": 10,
            "cd_njobs": 1,
            "zero_cn_thres": 0.005,
        }

    @pytest.mark.parametrize(
        "reg_term", ["RAW", "MAXCN", "DBOX_L1", "DBOX_L0", "DROOT_SUM", "DADJ_SUM"]
    )
    def test_cd_reg_term(self, cd_base_args, reg_term, tmp_path):
        """CD mode should produce a valid UCN file for each reg term."""
        from hatchet.compute_cn.compute_cn import run as run_compute_cn

        result_dir = str(tmp_path / f"cd_{reg_term}")
        args = {**cd_base_args, "result_dir": result_dir, "reg_term": reg_term}
        run_compute_cn(args)

        ucn = os.path.join(result_dir, "best.bbc.ucn")
        assert os.path.isfile(ucn), f"CD with {reg_term} did not produce best.bbc.ucn"
        bbc = pd.read_table(ucn, sep="\t")
        assert (bbc["cn_normal"] == "1|1").all(), (
            f"Normal clone not 1|1 with {reg_term}"
        )


class TestCntCDMode:
    """Test CNT-CD solver mode through the full pipeline."""

    def test_cnt_cd(self, cluster_bins_result, synthetic_data, tmp_path):
        from hatchet.compute_cn.compute_cn import run as run_compute_cn

        bbc_dir, _ = cluster_bins_result
        _, genome_sizes, regions_bed, _ = synthetic_data
        result_dir = str(tmp_path / "cnt_cd")
        # Only override defaults that differ from src/hatchet/hatchet.yaml.
        # tree_file is a path arg (kept in argparse, not YAML) — must be supplied here.
        args = {
            "bbc": os.path.join(bbc_dir, "bulk.bbc"),
            "seg": os.path.join(bbc_dir, "bulk.seg"),
            "result_dir": result_dir,
            "genome_size": genome_sizes,
            "region_bed": regions_bed,
            "mode": "cnt_cd",
            "solver": "cbc",
            "timelimit": 30,
            "maxClone": 2,
            "diploid": True,
            "force": True,
            "reg_term": "RAW",
            "reg_steps": 1,
            "reg_bound": 0.0,
            "diploidcmax": 6,
            "cd_niters": 5,
            "cd_nseeds": 10,
            "cd_njobs": 1,
            "zero_cn_thres": 0.005,
            "tree_file": None,
        }
        run_compute_cn(args)

        ucn = os.path.join(result_dir, "best.bbc.ucn")
        assert os.path.isfile(ucn), "cnt_cd did not produce best.bbc.ucn"
        bbc = pd.read_table(ucn, sep="\t")
        assert (bbc["cn_normal"] == "1|1").all(), "Normal clone not 1|1"
        u_cols = [c for c in bbc.columns if c.startswith("u_")]
        prop_sums = bbc[u_cols].sum(axis=1)
        assert np.allclose(prop_sums, 1.0, atol=0.05), (
            f"Proportions don't sum to 1: [{prop_sums.min():.3f}, {prop_sums.max():.3f}]"
        )
