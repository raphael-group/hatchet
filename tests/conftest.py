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

    args = {
        "bb_dir": bb_dir,
        "bbc_dir": bbc_dir,
        "genome_size": genome_sizes,
        "minK": 3,
        "maxK": 5,
        "restarts": 3,
        "top_restarts": 2,
        "n_local_trials": 2,
        "niters": 5,
        "t": 1e-6,
        "min_tau": 50,
        "max_tau": 200,
        "baf_eps": 1e-3,
        "min_covar": 1e-3,
        "ig_alpha": 10.0,
        "tau_iters": 1,
        "seed": 42,
        "decode_method": "viterbi",
        "score_method": "icl",
        "log_rdr": False,
        "init_method": "cna_plus_plus",
        "force": True,
        "verbosity": 1,
        "bal_lrt_alpha": 0.05,
        "bal_margin": 0.03,
        "filter_std": 2.0,
        "min_nbins": 10,
        "ub_nbins": 50,
        "skip_mhbafs": False,
        "training_method": "baum_welch",
        "free_baf_c0": False,
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

    args = {
        "bbc": bbc_file,
        "seg": seg_file,
        "result_dir": result_dir,
        "genome_size": genome_sizes,
        "region_bed": regions_bed,
        "mode": "ilp",
        "solver": "cbc",
        "timelimit": 60,
        "minClone": 2,
        "maxClone": 2,
        "diploid": True,
        "tetraploid": False,
        "fcn_ci_alpha": 0.05,
        "min_ci_margin": 0.1,
        "obj_type": "imf",
        "model_select": "bic",
        "force": True,
        "reg_term": "MAXCN",
        "reg_steps": 3,
        "reg_bound": 0.15,
        "no_ampdel": False,
        "num_cnstates": -1,
        "diploidcmax": 6,
        "tetraploidcmax": 12,
        "min_prop": 0.01,
        "purities": None,
        "cd_niters": 10,
        "cd_convergence_iters": 2,
        "cd_nseeds": 400,
        "cd_njobs": 1,
        "cd_seed": 42,
        "u_init": "dirichlet",
        "u_dir_alpha": 0.3,
        "solver_threads": None,
        "zero_cn_thres": 0.005,
        "cd_tol": 0.001,
        "fix_cn_dip": {},
        "fix_cn_tet": {},
        "verbosity": 1,
    }
    run_compute_cn(args)
    return result_dir, ground_truth
