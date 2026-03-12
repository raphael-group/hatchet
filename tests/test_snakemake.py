"""Integration test for the full HATCHet3 Snakemake pipeline."""

import os
import subprocess

import yaml
import pandas as pd
import pytest

from conftest import CBC_AVAILABLE
from simulate_bb_dir import simulate_bb_dir

pytestmark = pytest.mark.skipif(not CBC_AVAILABLE, reason="CBC solver not available")

SNAKEFILE = os.path.join(os.path.dirname(__file__), "..", "Snakefile")


@pytest.fixture(scope="module")
def snakemake_result(tmp_path_factory):
    """Run the full Snakemake pipeline on synthetic data."""
    base = tmp_path_factory.mktemp("snakemake")
    bb_dir = str(base / "bb_dir")
    genome_sizes = str(base / "genome.sizes")
    regions_bed = str(base / "regions.bed")
    workdir = str(base / "output")
    os.makedirs(workdir, exist_ok=True)

    ground_truth = simulate_bb_dir(bb_dir, genome_sizes, regions_bed)

    config = {
        "bb_dir": bb_dir,
        "genome_size": genome_sizes,
        "region_bed": regions_bed,
        "bbc_dir": "bbc",
        "result_dir": "results",
        "plot_dir": "summary",
        "log_dir": "logs",
        "verbosity": 1,
        "threads": 1,
        "cluster_bins": {
            "minK": 3, "maxK": 5, "t": 1e-6,
            "restarts": 3, "top_restarts": 2, "niters": 5,
            "bb_quantile": 0.2, "min_tau": 50, "max_tau": 500,
            "baf_eps": 1e-3, "min_covar": 1e-3, "ig_alpha": 10.0,
            "tau_iters": 1, "log_rdr": False,
            "decode_method": "viterbi", "score_method": "icl",
            "init_method": "cna_plus_plus",
        },
        "compute_cn": {
            "k": None, "mode": "ilp", "solver": "cbc", "timelimit": 60,
            "filter_cluster": False, "filter_std": 2.0,
            "balanced_baf_tol": 0.04, "toleranceRDR": 0.08, "toleranceBAF": 0.04,
            "minClone": 2, "maxClone": 2,
            "diploid": True, "tetraploid": False,
            "reg_term": "MAXCN", "reg_steps": 3, "reg_stepsize": 0.01,
            "no_ampdel": False, "num_cnstates": -1,
            "diploidcmax": 6, "tetraploidcmax": 12, "min_prop": 0.01,
            "purities": None,
            "cd_niters": 10, "cd_convergence_iters": 2,
            "cd_nseeds": 400, "cd_njobs": 1, "cd_seed": 42,
        },
        "plot_cn": {
            "dpi": 100, "img_type": "png", "transparent": False,
            "keep_gap": False, "tail_alpha": 0.8, "center_alpha": 1.0,
            "onetail_area": 0.025, "maxlim_fcn": 30,
        },
    }

    config_file = str(base / "config.yaml")
    with open(config_file, "w") as f:
        yaml.dump(config, f)

    cmd = [
        "snakemake", "-p", "--cores", "1",
        "-s", os.path.abspath(SNAKEFILE),
        "--configfile", config_file,
        "--directory", workdir,
    ]

    result = subprocess.run(
        cmd, capture_output=True, text=True, timeout=300,
    )

    return workdir, ground_truth, result


class TestSnakemakePipeline:
    """Verify the full Snakemake pipeline runs and produces correct outputs."""

    def test_snakemake_exit_code(self, snakemake_result):
        workdir, _, result = snakemake_result
        assert result.returncode == 0, (
            f"Snakemake failed with exit code {result.returncode}\n"
            f"STDOUT:\n{result.stdout[-2000:]}\n"
            f"STDERR:\n{result.stderr[-2000:]}"
        )

    def test_bbc_output_exists(self, snakemake_result):
        workdir, _, result = snakemake_result
        if result.returncode != 0:
            pytest.skip("Snakemake failed")
        assert os.path.isfile(os.path.join(workdir, "bbc", "bulk.bbc"))
        assert os.path.isfile(os.path.join(workdir, "bbc", "bulk.seg"))

    def test_best_ucn_exists(self, snakemake_result):
        workdir, _, result = snakemake_result
        if result.returncode != 0:
            pytest.skip("Snakemake failed")
        assert os.path.isfile(os.path.join(workdir, "results", "best.bbc.ucn"))
        assert os.path.isfile(os.path.join(workdir, "results", "best.seg.ucn"))

    def test_gammas_exists(self, snakemake_result):
        workdir, _, result = snakemake_result
        if result.returncode != 0:
            pytest.skip("Snakemake failed")
        assert os.path.isfile(os.path.join(workdir, "results", "gammas.tsv"))

    def test_ucn_format(self, snakemake_result):
        workdir, _, result = snakemake_result
        if result.returncode != 0:
            pytest.skip("Snakemake failed")
        bbc = pd.read_table(
            os.path.join(workdir, "results", "best.bbc.ucn"), sep="\t"
        )
        assert "cn_normal" in bbc.columns
        assert "cn_clone1" in bbc.columns
        assert (bbc["cn_normal"] == "1|1").all()

    def test_logs_created(self, snakemake_result):
        workdir, _, result = snakemake_result
        if result.returncode != 0:
            pytest.skip("Snakemake failed")
        assert os.path.isfile(os.path.join(workdir, "logs", "cluster_bins.log"))
        assert os.path.isfile(os.path.join(workdir, "logs", "compute_cn.log"))
