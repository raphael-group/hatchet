# HATCHet: Holistic Allele-specific Tumor Copy-number Heterogeneity

<!-- [![CI](https://github.com/RunpengLuo/hatchet-long-read/actions/workflows/ci.yml/badge.svg?branch=hatchet3-dev)](https://github.com/RunpengLuo/hatchet-long-read/actions/workflows/ci.yml) -->
<!-- [![codecov](https://codecov.io/gh/RunpengLuo/hatchet-long-read/branch/hatchet3-dev/graph/badge.svg)](https://codecov.io/gh/RunpengLuo/hatchet-long-read) -->
[![Version](https://img.shields.io/badge/version-3.0.0b1-blue.svg)](VERSION)
[![Python](https://img.shields.io/badge/python-3.11%2B-blue.svg)](https://www.python.org/downloads/)
[![Code style: Ruff](https://img.shields.io/badge/code%20style-ruff-000000.svg)](https://github.com/astral-sh/ruff)
![Status](https://img.shields.io/badge/status-under%20development-yellow)

HATCHet is an algorithm to infer allele- and clone-specific copy-number aberrations (CNAs), clone proportions, and whole-genome duplications (WGD) for several tumor clones jointly from multiple bulk-tumor samples of the same patient. HATCHet supports both **short-read** WGS/WES and **long-read** (e.g., PacBio HiFi, Oxford Nanopore) sequencing data.

## Table of Contents
1. [Installation](#installation) <br>
    - [Option 1. Quick Install](#option-1-quick-install-recommended) <br>
    - [Option 2. Manual Install](#option-2-manual-install) <br>
    - [Setup ILP Solver](#setup-ilp-solver) <br>
2. [Usage](#usage) <br>
3. [Documentation](#documentation) <br>

## Installation

HATCHet requires a 64-bit Linux or macOS system and Python 3.11 or higher.

### Option 1. Quick Install (Recommended)

Install HATCHet from Bioconda channel:  

```sh
# or conda
mamba install -c conda-forge -c bioconda hatchet3
```

> [!IMPORTANT]
> TODO: HATCHet (v3) is not yet published on bioconda. Use Option 2 instead.

### Option 2. Manual Install

For developers, build the conda environment and install HATCHet via `Pip`. `--no-build-isolation` ensures pip do not re-fetch dependencies from PyPI. Add `-e` to install in editable mode (`pip install --no-build-isolation -e .`).

```sh
mamba env create -f ./environment.yaml -p /path/to/envs/hatchet-env
conda activate /path/to/envs/hatchet-env
pip install --no-build-isolation .
```

> [!NOTE]
> The `pip install` step compiles the C++ pybind11 extension. OpenMP support is provided per platform by `cxx-compiler` (`llvm-openmp` on macOS, `libgomp` on Linux), so no OpenMP package is pinned explicitly. Verify the install with `hatchet --help`.

### Setup ILP Solver

HATCHet requires a [Pyomo](https://pyomo.readthedocs.io/) supported solver for integer copy-number factorization step (`compute-cn`, via option `--solver`).

> [!NOTE]
> the conda environment already includes both `cbc>=2.10` and `gurobi>=13.0.0`.

#### Gurobi

**Gurobi** is a commercial optimization solver and we recommend user to use the latest version if available. The following table lists compatible Python and Gurobi versions:

| Python Version | Compatible Gurobi Versions | Notes |
|---|---|---|
| **3.11** | **10.0, 11.0, 12.0, 12.0.1-12.0.3, 13.0.0, 13.0.1** | **Minimum supported** |
| 3.12 | 11.0, 12.0, 12.0.1-12.0.3, 13.0.0, 13.0.1 | |

User can obtain the license via two options:
- a single-host license where the license is tied to a single computer.
- a network license for use in a compute cluster (using a license server in the cluster).
Both options are freely and [easily available](http://www.gurobi.com/academia/academia-center) for users in academia. Set the `GRB_LICENSE_FILE` environment variable to point to your license file:

```sh
export GRB_LICENSE_FILE="/path/to/gurobi.lic"
```

#### CBC

**CBC** is an open-source [COIN-OR Branch-and-Cut](https://github.com/coin-or/Cbc) solver, user can directly set  `--solver cbc` when running `compute-cn`.


## Usage
HATCHet (v3) inputs genomic bin by sample read-depth ratio (RDR), phased B-allele counts, and total-allele count matrices preprocessed by [Universal-Genotyping-Pipeline](https://github.com/raphael-group/Universal-Genotyping-Pipeline), see [tutorial](https://github.com/raphael-group/Universal-Genotyping-Pipeline/docs/bulk_genotyping.md) for preprocessing details and [Input](./docs/reference.md#input) for input data formats.

### Running Snakemake Pipeline
We include a [Snakemake](https://snakemake.readthedocs.io/) pipeline (version 9 or newer) that runs the full HATCHet (v3) pipeline for a single patient. First, copy and modify the Snakemake configuration file from [config/snakemake-hatchet.yaml](config/snakemake-hatchet.yaml):
```yaml
patient_id: sample                # output filename prefix
bb_dir: "/path/to/bb"             # input directory of preprocessed matrices (see Input)
genome_size: "/path/to/hg38.chrom.sizes"   # reference chromosome sizes
region_bed: "/path/to/regions.bed"         # BED of whitelist regions
```

See [config/snakemake-hatchet.yaml](config/snakemake-hatchet.yaml) for output paths and per-step parameters
(`cluster-bins`, `compute-cn`). Then run:

```sh
snakemake -p --cores <ncores> -s ./Snakefile \
    --configfile /path/to/my_config.yaml \
    --directory <output_dir>
```

> [!TIP]
> Always preview first with a **dry run** - add `--dry-run` (`-n`). It lists the jobs Snakemake would run without running any, so you can confirm the plan and catch config or input-path mistakes early:
> ```sh
> snakemake -p --cores <ncores> -s ./Snakefile \
>     --configfile /path/to/my_config.yaml \
>     --directory <output_dir> \
>     -n
> ```

### Modules

| Order | Step | Description |
|---|---|---|
| (1) | [`cluster-bins`](docs/modules/cluster-bins.md) | Local-global genome segmentation using a Gaussian RDR + Beta-Binomial BAF multi-sample factorial HMM with phase switch correction. |
| (2) | [`compute-cn`](docs/modules/compute-cn.md) | Allele-specific integer copy numbers and clone proportions deconvolution with regularization using integer linear programming (ILP) or coordinate descent. |
| (3) | [`plot-cn`](docs/modules/plot-cn.md) | Genome-wide copy-number profiles and RDR-vs-BAF plots for a single CN solution. |
| (4) | [`plot-panel`](docs/modules/plot-panel.md) | Multi-sample copy-number panel composed from several per-sample CN solutions. |

## Documentation

| Document | Description |
|----------|-------------|
| [docs/reference.md](docs/reference.md) | Output directory structure and (hyper-)parameter reference |
| [CHANGELOG.md](CHANGELOG.md) | Change logs |
| [config/snakemake-hatchet.yaml](config/snakemake-hatchet.yaml) | Snakemake pipeline configuration |
