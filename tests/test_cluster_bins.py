"""Integration tests for HATCHet3 cluster-bins step."""

import os
import numpy as np
import pandas as pd
import pytest


class TestClusterBinsOutputFiles:
    """Verify cluster-bins produces expected output files."""

    def test_bbc_file_exists(self, cluster_bins_result):
        bbc_dir, _ = cluster_bins_result
        assert os.path.isfile(os.path.join(bbc_dir, "bulk.bbc"))

    def test_seg_file_exists(self, cluster_bins_result):
        bbc_dir, _ = cluster_bins_result
        assert os.path.isfile(os.path.join(bbc_dir, "bulk.seg"))

    def test_model_scores_exists(self, cluster_bins_result):
        bbc_dir, _ = cluster_bins_result
        assert os.path.isfile(os.path.join(bbc_dir, "model_scores.tsv"))


class TestClusterBinsOutputFormat:
    """Verify BBC and SEG file formats are correct."""

    def test_bbc_columns(self, cluster_bins_result):
        bbc_dir, _ = cluster_bins_result
        bbc = pd.read_table(os.path.join(bbc_dir, "bulk.bbc"), sep="\t")
        required = {
            "#CHR",
            "START",
            "END",
            "SAMPLE",
            "#SNPS",
            "CLUSTER",
            "RD",
            "COV",
            "BAF",
        }
        assert required.issubset(set(bbc.columns)), (
            f"Missing columns: {required - set(bbc.columns)}"
        )

    def test_bbc_rdr_positive(self, cluster_bins_result):
        bbc_dir, _ = cluster_bins_result
        bbc = pd.read_table(os.path.join(bbc_dir, "bulk.bbc"), sep="\t")
        assert (bbc["RD"] > 0).all(), "All RDR values should be positive"

    def test_bbc_baf_range(self, cluster_bins_result):
        bbc_dir, _ = cluster_bins_result
        bbc = pd.read_table(os.path.join(bbc_dir, "bulk.bbc"), sep="\t")
        assert (bbc["BAF"] >= 0).all() and (bbc["BAF"] <= 1).all(), (
            "BAF should be in [0, 1]"
        )

    def test_seg_columns(self, cluster_bins_result):
        bbc_dir, _ = cluster_bins_result
        seg = pd.read_table(os.path.join(bbc_dir, "bulk.seg"), sep="\t")
        required = {"#ID", "SAMPLE", "#BINS", "#SNPS", "LENGTH", "BAF", "RD"}
        assert required.issubset(set(seg.columns)), (
            f"Missing columns: {required - set(seg.columns)}"
        )

    def test_model_scores_format(self, cluster_bins_result):
        bbc_dir, _ = cluster_bins_result
        scores = pd.read_table(os.path.join(bbc_dir, "model_scores.tsv"), sep="\t")
        assert "K" in scores.columns
        assert "bic" in scores.columns or "icl" in scores.columns


class TestClusterBinsClusterSeparation:
    """Verify clusters recover distinct CN states."""

    def test_at_least_3_clusters(self, cluster_bins_result):
        bbc_dir, _ = cluster_bins_result
        seg = pd.read_table(os.path.join(bbc_dir, "bulk.seg"), sep="\t")
        n_clusters = seg["#ID"].nunique()
        assert n_clusters >= 3, f"Expected >= 3 clusters, got {n_clusters}"

    def test_cluster_rdr_spread(self, cluster_bins_result):
        """Cluster median RDRs should span a range > 0.3 (distinct CN states)."""
        bbc_dir, _ = cluster_bins_result
        bbc = pd.read_table(os.path.join(bbc_dir, "bulk.bbc"), sep="\t")
        medians = bbc.groupby("CLUSTER")["RD"].median()
        rdr_range = medians.max() - medians.min()
        assert rdr_range > 0.3, (
            f"Cluster RDR range {rdr_range:.3f} too narrow (expected > 0.3)"
        )
