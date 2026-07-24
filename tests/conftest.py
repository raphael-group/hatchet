"""Pytest fixtures for HATCHet3 integration tests."""

import json
import os

import matplotlib

matplotlib.use("Agg")

import pytest

DATA_DIR = os.path.join(os.path.dirname(__file__), "data")


@pytest.fixture(scope="session")
def synthetic_data():
    """Load pre-generated synthetic data from tests/data/."""
    bb_dir = os.path.join(DATA_DIR, "bb_dir")
    genome_sizes = os.path.join(DATA_DIR, "genome.sizes")
    regions_bed = os.path.join(DATA_DIR, "regions.bed")
    with open(os.path.join(DATA_DIR, "ground_truth.json")) as f:
        ground_truth = json.load(f)
    return bb_dir, genome_sizes, regions_bed, ground_truth


@pytest.fixture(scope="session")
def cluster_bins_result(synthetic_data, tmp_path_factory):
    """Run cluster-bins once per session, return (bbc_dir, ground_truth)."""
    bb_dir, genome_sizes, regions_bed, ground_truth = synthetic_data
    bbc_dir = str(tmp_path_factory.mktemp("cluster_bins_out"))

    from hatchet.cluster_bins.cluster_bins import run as run_cluster_bins

    # Only override defaults that differ from src/hatchet/hatchet.yaml.
    args = {
        "bb_dir": bb_dir,
        "bbc_dir": bbc_dir,
        "genome_size": genome_sizes,
        "region_bed": regions_bed,
        "maxK": 5,
        "restarts": 3,
        "top_restarts": 2,
        "n_local_trials": 2,
        "niters": 5,
        "tau_iters": 1,
        "decode_method": "viterbi",
        "force": True,
        "verbosity": 1,
    }
    run_cluster_bins(args)
    return bbc_dir, ground_truth


def _check_cbc_available():
    try:
        from pyomo import environ as pe

        return pe.SolverFactory("cbc").available(exception_flag=False)
    except Exception:
        return False


CBC_AVAILABLE = _check_cbc_available()


@pytest.fixture(scope="session")
def compute_cn_result(cluster_bins_result, synthetic_data, tmp_path_factory):
    """Run compute-cn once per session, return (result_dir, ground_truth)."""
    if not CBC_AVAILABLE:
        pytest.skip("CBC solver not available")

    bbc_dir, ground_truth = cluster_bins_result
    bb_dir, genome_sizes, regions_bed, _ = synthetic_data
    result_dir = str(tmp_path_factory.mktemp("compute_cn_out"))

    bbc_file = os.path.join(bbc_dir, "bulk.bbc")
    seg_file = os.path.join(bbc_dir, "bulk.seg")

    from hatchet.compute_cn.compute_cn import run as run_compute_cn

    # Only override defaults that differ from src/hatchet/hatchet.yaml.
    args = {
        "bbc": bbc_file,
        "seg": seg_file,
        "result_dir": result_dir,
        "genome_size": genome_sizes,
        "region_bed": regions_bed,
        "solver": "cbc",
        "timelimit": 60,
        "maxClone": 2,
        "diploid": True,
        "force": True,
        "reg_steps": 3,
        "diploidcmax": 6,
        "cd_njobs": 1,
        "zero_cn_thres": 0.005,
        "verbosity": 1,
    }
    run_compute_cn(args)
    return result_dir, ground_truth
