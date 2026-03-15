"""Simulate a synthetic bb_dir for HATCHet3 integration tests.

Generates NPZ count matrices and metadata files that match the exact format
expected by cluster-bins (see cluster_bins.py:84-141).

Default scenario: 1 tumor sample, 1 tumor clone at 80% purity, with
4 genomic regions (chr1 p-arm, chr1 q-arm, chr22 p-arm, chr22 q-arm)
carrying distinct CN states.
"""

import argparse
import os

import numpy as np
import pandas as pd


# hg38 chromosome sizes (subset)
DEFAULT_GENOME_SIZES = {
    "chr1": 248956422,
    "chr22": 50818468,
}

# hg38 arm boundaries (approximate, from UCSC cytoBand)
DEFAULT_REGIONS = [
    ("chr1", 0, 121700000, "chr1_p-arm"),
    ("chr1", 121700000, 248956422, "chr1_q-arm"),
    ("chr22", 0, 13000000, "chr22_p-arm"),
    ("chr22", 13000000, 50818468, "chr22_q-arm"),
]

# Default tumor CN profile: (cn_a, cn_b) per region
DEFAULT_CN_PROFILE = {
    "chr1_p-arm": (1, 1),  # balanced diploid
    "chr1_q-arm": (2, 1),  # gain
    "chr22_p-arm": (1, 1),  # balanced diploid
    "chr22_q-arm": (1, 0),  # LOH
}

DEFAULT_CLONE_PROPORTIONS = [0.2, 0.8]  # [normal, tumor1]

DEFAULT_NOISE_PARAMS = {
    "sigma_rdr": 0.03,
    "tau_baf": 100.0,
    "depth_lambda": 100,
    "normal_depth": 30.0,
}


def simulate_bb_dir(
    bb_dir,
    genome_sizes_file,
    regions_bed_file,
    bin_size=1_000_000,
    genome_sizes=None,
    regions=None,
    cn_profile=None,
    clone_proportions=None,
    noise_params=None,
    seed=42,
):
    """Generate a synthetic bb_dir with NPZ matrices and metadata.

    Returns a ground-truth dict with per-region CN states, purity, and gamma.
    """
    genome_sizes = genome_sizes or DEFAULT_GENOME_SIZES
    regions = regions or DEFAULT_REGIONS
    cn_profile = cn_profile or DEFAULT_CN_PROFILE
    clone_proportions = clone_proportions or DEFAULT_CLONE_PROPORTIONS
    noise_params = {**DEFAULT_NOISE_PARAMS, **(noise_params or {})}
    rng = np.random.RandomState(seed)

    os.makedirs(bb_dir, exist_ok=True)

    sigma_rdr = noise_params["sigma_rdr"]
    tau_baf = noise_params["tau_baf"]
    depth_lambda = noise_params["depth_lambda"]
    normal_depth = noise_params["normal_depth"]

    prop_normal = clone_proportions[0]
    prop_tumor = clone_proportions[1]

    # Build bins
    bins = []
    for chrom, start, end, name in regions:
        n_bins = (end - start) // bin_size
        for i in range(n_bins):
            b_start = start + i * bin_size
            b_end = b_start + bin_size
            bins.append((chrom, b_start, b_end, name))

    n_bins_total = len(bins)

    # Arrays for NPZ matrices
    # Columns: [normal, tumor1] → 2 columns
    rdr_tumor = np.zeros((n_bins_total, 1), dtype=np.float32)
    depth_mat = np.zeros((n_bins_total, 2), dtype=np.float32)
    a_allele = np.zeros((n_bins_total, 2), dtype=np.int32)
    b_allele = np.zeros((n_bins_total, 2), dtype=np.int32)
    t_allele = np.zeros((n_bins_total, 2), dtype=np.int32)

    chrs = []
    starts = []
    ends = []
    region_ids = []
    switchprobs = []
    snps_list = []

    for i, (chrom, b_start, b_end, region_name) in enumerate(bins):
        cn_a, cn_b = cn_profile[region_name]

        # Mixed signal (normal CN is always (1,1))
        true_rdr = prop_normal * 1.0 + prop_tumor * (cn_a + cn_b) / 2.0
        total_cn_mixed = prop_normal * 2 + prop_tumor * (cn_a + cn_b)
        true_baf = (prop_normal * 1 + prop_tumor * cn_b) / total_cn_mixed

        # Observed RDR with Gaussian noise
        obs_rdr = max(true_rdr + rng.normal(0, sigma_rdr), 0.01)
        rdr_tumor[i, 0] = obs_rdr

        # Depths
        depth_mat[i, 0] = normal_depth
        depth_mat[i, 1] = obs_rdr * normal_depth

        # BAF: Beta-Binomial noise
        total_reads = rng.poisson(depth_lambda)
        total_reads = max(total_reads, 10)  # floor

        # Beta-Binomial: draw p ~ Beta(tau*baf, tau*(1-baf)), then B ~ Binom(total, p)
        alpha_bb = tau_baf * true_baf
        beta_bb = tau_baf * (1 - true_baf)
        p_obs = rng.beta(max(alpha_bb, 0.01), max(beta_bb, 0.01))
        b_count = rng.binomial(total_reads, p_obs)
        a_count = total_reads - b_count

        # Tumor allele counts
        a_allele[i, 1] = a_count
        b_allele[i, 1] = b_count
        t_allele[i, 1] = total_reads

        # Normal allele counts (balanced)
        normal_total = rng.poisson(depth_lambda)
        normal_total = max(normal_total, 10)
        normal_b = rng.binomial(normal_total, 0.5)
        a_allele[i, 0] = normal_total - normal_b
        b_allele[i, 0] = normal_b
        t_allele[i, 0] = normal_total

        chrs.append(chrom)
        starts.append(b_start)
        ends.append(b_end)
        region_ids.append(region_name)
        switchprobs.append(1e-6)
        snps_list.append(total_reads)

    # Write bb.tsv.gz
    bb_df = pd.DataFrame(
        {
            "#CHR": chrs,
            "START": starts,
            "END": ends,
            "region_id": region_ids,
            "switchprobs": switchprobs,
            "#SNPS": snps_list,
        }
    )
    bb_df.to_csv(
        os.path.join(bb_dir, "bb.tsv.gz"),
        sep="\t",
        index=False,
        compression="gzip",
    )

    # Write sample_ids.tsv
    sample_df = pd.DataFrame(
        {
            "SAMPLE": ["normal", "tumor1"],
            "sample_type": ["normal", "tumor"],
        }
    )
    sample_df.to_csv(os.path.join(bb_dir, "sample_ids.tsv"), sep="\t", index=False)

    # Write NPZ files
    # bb.rdr.npz: tumor RDR only, shape (N, n_tumor_samples)
    np.savez(os.path.join(bb_dir, "bb.rdr.npz"), mat=rdr_tumor)

    # bb.depth.npz: shape (N, 2) [normal, tumor]
    np.savez(os.path.join(bb_dir, "bb.depth.npz"), mat=depth_mat)

    # Allele counts: shape (N, 2) [normal, tumor]
    np.savez(os.path.join(bb_dir, "bb.Aallele.npz"), mat=a_allele)
    np.savez(os.path.join(bb_dir, "bb.Ballele.npz"), mat=b_allele)
    np.savez(os.path.join(bb_dir, "bb.Tallele.npz"), mat=t_allele)

    # Write genome.sizes
    os.makedirs(os.path.dirname(genome_sizes_file) or ".", exist_ok=True)
    with open(genome_sizes_file, "w") as f:
        for chrom, size in genome_sizes.items():
            f.write(f"{chrom}\t{size}\n")

    # Write regions.bed
    os.makedirs(os.path.dirname(regions_bed_file) or ".", exist_ok=True)
    with open(regions_bed_file, "w") as f:
        for chrom, start, end, name in regions:
            f.write(f"{chrom}\t{start}\t{end}\t{name}\n")

    # Ground truth
    ground_truth = {
        "cn_profile": cn_profile,
        "clone_proportions": clone_proportions,
        "purity": prop_tumor,
        "gamma": 2.0 / 1.0,  # balanced cluster RDR=1.0 → gamma = 2/1 = 2.0
        "n_bins": n_bins_total,
        "regions": regions,
        "expected_rdrs": {},
        "expected_bafs": {},
    }
    for name, (cn_a, cn_b) in cn_profile.items():
        true_rdr = prop_normal * 1.0 + prop_tumor * (cn_a + cn_b) / 2.0
        total_cn = prop_normal * 2 + prop_tumor * (cn_a + cn_b)
        true_baf = (prop_normal * 1 + prop_tumor * cn_b) / total_cn
        ground_truth["expected_rdrs"][name] = true_rdr
        ground_truth["expected_bafs"][name] = true_baf

    return ground_truth


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Simulate synthetic bb_dir for HATCHet3"
    )
    parser.add_argument(
        "--output_dir",
        default="tests/fixtures/synthetic_bb",
        help="Output directory (bb_dir will be created inside)",
    )
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--bin_size", type=int, default=1_000_000)
    args = parser.parse_args()

    out = args.output_dir
    os.makedirs(out, exist_ok=True)
    gt = simulate_bb_dir(
        bb_dir=os.path.join(out, "bb_dir"),
        genome_sizes_file=os.path.join(out, "genome.sizes"),
        regions_bed_file=os.path.join(out, "regions.bed"),
        bin_size=args.bin_size,
        seed=args.seed,
    )
    print(f"Simulated {gt['n_bins']} bins in {out}/bb_dir")
    print(f"Ground truth CN profile: {gt['cn_profile']}")
    print(f"Expected RDRs: {gt['expected_rdrs']}")
    print(f"Expected BAFs: {gt['expected_bafs']}")
