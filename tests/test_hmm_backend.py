"""Verify both HMM backends: the C++ extension (_hmm_cpp) and the Numba fallback.

The backend is chosen by ``hatchet.cluster_bins.hmm._USE_CPP``. These tests assert
the C++ path is actually selected when the extension is compiled (guarding against a
silent fallback to Numba), and that forcing the Numba path still yields valid output.
"""

import os

import pandas as pd
import pytest

import hatchet.cluster_bins.hmm as hmm


def _run_cluster_bins(bb_dir, genome_sizes, regions_bed, out_dir):
    """Run cluster-bins on the synthetic data with the shared test settings."""
    from hatchet.cluster_bins.cluster_bins import run as run_cluster_bins

    args = {
        "bb_dir": bb_dir,
        "bbc_dir": out_dir,
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
    return out_dir


@pytest.mark.skipif(not hmm._CPP_IMPORTABLE, reason="_hmm_cpp extension not compiled")
def test_cpp_backend_active():
    """When the extension is compiled, the C++ backend must be selected."""
    assert hmm._USE_CPP is True, (
        "_hmm_cpp is importable but the Numba fallback was selected"
    )


def test_numba_backend_runs(synthetic_data, tmp_path, monkeypatch):
    """Force the Numba/scipy fallback and verify cluster-bins still runs correctly."""
    monkeypatch.setattr(hmm, "_USE_CPP", False)
    bb_dir, genome_sizes, regions_bed, _ = synthetic_data
    out = str(tmp_path / "numba_out")

    _run_cluster_bins(bb_dir, genome_sizes, regions_bed, out)

    bbc = pd.read_table(os.path.join(out, "bulk.bbc"), sep="\t")
    assert (bbc["RD"] > 0).all()
    assert (bbc["BAF"] >= 0).all() and (bbc["BAF"] <= 1).all()
    seg = pd.read_table(os.path.join(out, "bulk.seg"), sep="\t")
    assert seg["#ID"].nunique() >= 3


@pytest.mark.skipif(not hmm._CPP_IMPORTABLE, reason="_hmm_cpp extension not compiled")
def test_backends_agree_on_cluster_count(synthetic_data, tmp_path, monkeypatch):
    """C++ and Numba backends should recover a comparable number of clusters."""
    bb_dir, genome_sizes, regions_bed, _ = synthetic_data

    monkeypatch.setattr(hmm, "_USE_CPP", True)
    cpp_out = _run_cluster_bins(
        bb_dir, genome_sizes, regions_bed, str(tmp_path / "cpp")
    )

    monkeypatch.setattr(hmm, "_USE_CPP", False)
    nb_out = _run_cluster_bins(bb_dir, genome_sizes, regions_bed, str(tmp_path / "nb"))

    n_cpp = pd.read_table(os.path.join(cpp_out, "bulk.seg"), sep="\t")["#ID"].nunique()
    n_nb = pd.read_table(os.path.join(nb_out, "bulk.seg"), sep="\t")["#ID"].nunique()
    assert abs(n_cpp - n_nb) <= 1, f"cluster counts diverge: C++={n_cpp}, Numba={n_nb}"
