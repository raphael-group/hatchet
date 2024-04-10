import numpy as np
import pandas as pd
from scipy.special import softmax
from scipy.stats import binom, norm


def compute_baf_task(bin_snps, blocksize, max_snps_per_block, test_alpha, multisample_em):
    """
    Estimates the BAF for the bin containing exactly <bin_snps> SNPs.
    <bin_snps> is a dataframe with at least ALT and REF columns containing read counts.
    <blocksize>, <max_snps_per_block>, and <test_alpha> are used only for constructing phase blocks.
    """

    samples = sorted(bin_snps.SAMPLE.unique())
    bafs, alphas, betas = [0] * len(samples), [0] * len(samples), [0] * len(samples)
    result = {}

    phasing = 'PHASE' in bin_snps.columns
    if phasing:
        # TODO: select the highest coverage sample to use for constructing phase blocks?
        # or maybe use all samples for BAF test?
        all_phase_data = [
            phase_blocks_sequential(
                d,
                blocksize=blocksize,
                max_snps_per_block=max_snps_per_block,
                alpha=test_alpha,
            )
            for _, d in bin_snps.groupby('SAMPLE')
        ]
        phase_data = merge_phasing(bin_snps, all_phase_data)
        grouped_snps = collapse_blocks(bin_snps, *phase_data, bin_snps.iloc[0].CHR)
    else:
        # each snp is it's on haplo block. [1] represents the minor count.
        # After EM algorithm, the haplostring tracks the haplotype with the minor
        # allele count (beta)
        grouped_snps = bin_snps.copy()
        grouped_snps['START'] = grouped_snps['POS']
        grouped_snps['HAPLO'] = len(grouped_snps) * [[1]]
    reference_based_haplo_blocks = grouped_snps[grouped_snps.SAMPLE == samples[0]].HAPLO

    totals = grouped_snps.pivot(index='SAMPLE', columns='START', values='TOTAL').to_numpy().astype(np.uint64)
    alts = grouped_snps.pivot(index='SAMPLE', columns='START', values='ALT').to_numpy().astype(np.uint64)

    # compute phasing once using all samples
    if multisample_em:
        bafs, alphas, betas, _, haploflat, haplostring = apply_EM(
            totals, alts, reference_based_haplo_blocks, multisample_em
        )

    for i in range(len(samples)):
        sample = samples[i]
        sample_snps = bin_snps[bin_snps.SAMPLE == sample]

        # compute phases if not multisample
        if not multisample_em:
            bafs[i], alphas[i], betas[i], _, haploflat, haplostring = apply_EM(
                totals[i].reshape(1, -1), alts[i].reshape(1, -1), reference_based_haplo_blocks, multisample_em
            )
            bafs[i] = bafs[i][0]
            alphas[i] = alphas[i][0]
            betas[i] = betas[i][0]

        assert np.sum(np.choose(haploflat, [sample_snps.ALT, sample_snps.REF])) == betas[i]

        snp_pos = ','.join(map(str, sample_snps.POS))
        snp_ref_counts = ','.join(map(str, sample_snps.REF))
        snp_alt_counts = ','.join(map(str, sample_snps.ALT))
        n_snps = len(sample_snps)

        baf = bafs[i]
        cov = (alphas[i] + betas[i]) / n_snps
        tot = alphas[i] + betas[i]

        result[sample] = (
            n_snps,
            cov,
            baf,
            alphas[i],
            betas[i],
            tot,
            snp_pos,
            snp_ref_counts,
            snp_alt_counts,
            haplostring,
        )
    return result


def EM(totals_in, alts_in, start, tol=1e-6):
    """
    Adapted from chisel/Combiner.py
    """
    totals = np.array(totals_in).reshape(-1)
    alts = np.array(alts_in).reshape(-1)
    assert totals.size == alts.size and 0 < start < 1 and np.all(totals >= alts)
    if np.all(np.logical_or(totals == alts, alts == 0)):
        return 0.0, np.ones(len(totals_in)) * 0.5, 0.0
    else:
        baf = start
        prev = None
        while prev is None or abs(prev - baf) >= tol:
            prev = baf

            # E-step
            assert 0 + tol < baf < 1 - tol, (baf, totals, alts, start)
            M = (totals - 2.0 * alts) * np.log(baf) + (2.0 * alts - totals) * np.log(1.0 - baf)
            M = M.astype(float)
            M = np.exp(np.clip(a=M, a_min=-100, a_max=100))
            phases = np.reciprocal(1 + M)

            # M-step
            baf = float(np.sum(totals * (1 - phases) + alts * (2.0 * phases - 1))) / float(np.sum(totals))

        assert 0 + tol < baf < 1 - tol, (baf, totals, alts, start)
        lpmf = binom.logpmf
        log_likelihood = float(
            np.sum(
                phases * lpmf(k=list(alts), n=list(totals), p=baf)
                + (1 - phases) * lpmf(k=list(alts), n=list(totals), p=1 - baf)
            )
        )
        return [baf], phases, log_likelihood


def multisample_EM(totals, alts, start, tol=10e-6):

    refs = totals - alts
    assert refs.shape == alts.shape, 'Alternate and reference count arrays must have the same shape'
    assert 0 < start <= 0.5, 'Initial estimate must be in (0, 0.5]'

    n_samples, n_snps = alts.shape
    totals_sum = np.sum(totals, axis=1)
    assert np.all(totals_sum > 0), 'Every bin must have >0 SNP-covering reads in each sample!'

    e_arr = np.zeros((2, n_snps))

    if np.all(np.logical_or(refs == 0, alts == 0)):
        return np.array([0.0] * n_samples), np.ones(n_snps) * 0.5, 0.0
    else:
        theta = np.repeat(start, n_samples)
        prev = None
        while prev is None or np.all(np.abs(prev - theta) >= tol):
            prev = theta

            # Ensure that theta is not exactly 0 or 1 to avoid log(0)
            theta = np.clip(theta, tol, 1 - tol)

            # E-step in log space
            e_arr[0] = np.log(theta) @ refs + np.log(1 - theta) @ alts
            e_arr[1] = np.log(1 - theta) @ refs + np.log(theta) @ alts
            phases = softmax(e_arr, axis=0)[0, :]
            assert not np.any(np.isnan(phases)), (phases, e_arr)

            # M-step
            t1 = refs @ phases + alts @ (1 - phases)
            theta = t1 / totals_sum
            assert not np.any(np.isnan(theta)), (theta, t1)

    # If mean(BAF) > 0.5, flip phases accordingly
    if np.mean(theta) > 0.5:
        theta = np.clip(theta, tol, 1 - tol)

        theta = 1 - theta
        t1 = np.sum(np.log(theta) @ refs + np.log(1 - theta) @ alts, axis=0)
        t2 = np.sum(np.log(1 - theta) @ refs + np.log(theta) @ alts, axis=0)
        phases = softmax(np.vstack([t1, t2]), axis=0)[0, :]

    t1 = np.sum(np.log(theta) * (phases * refs).T, axis=0)
    t2 = np.sum(np.log(1 - theta) * (phases * alts).T, axis=0)
    t3 = np.sum(np.log(theta) * ((1 - phases) * alts).T, axis=0)
    t4 = np.sum(np.log(1 - theta) * ((1 - phases) * refs).T, axis=0)
    log_likelihood = np.sum(t1 + t2 + t3 + t4)

    return theta, phases, log_likelihood


def apply_EM(totals, alts, reference_based_haplo_blocks, multisample_em):
    if multisample_em:
        runs = (multisample_EM(totals, alts, b) for b in np.arange(0.05, 0.5, 0.05))
    else:
        runs = (EM(totals, alts, start=st) for st in np.linspace(0.01, 0.49, 50))
    baf, phases, _ = max(runs, key=lambda x: x[-1])
    refs = totals - alts
    phases = phases.round().astype(np.int8).reshape(-1)
    inverse_reference_haplo = pd.Series([[1 - ph for ph in hap] for hap in reference_based_haplo_blocks])
    haplo = np.choose(phases, [reference_based_haplo_blocks, inverse_reference_haplo])
    haploflat = [int(item) for sublist in haplo for item in sublist]
    haplostring = ','.join(list(map(str, haploflat)))
    # if not np.all(phases == (alts < refs)):
    #     sp.log(totals)
    return (
        baf,
        [int(np.sum(np.choose(phases, [alts[i], refs[i]]))) for i in range(len(alts))],
        [int(np.sum(np.choose(phases, [refs[i], alts[i]]))) for i in range(len(alts))],
        haplo,
        haploflat,
        haplostring,
    )


def merge_phasing(_, all_phase_data):
    """
    Merge phasing results across all samples:
    if a pair of SNPs is split in any sample, they won't be split.
    """

    if len(all_phase_data) == 1:
        return all_phase_data[0]

    singletons = set()
    orphans = set()
    breaks = set()
    for i in range(len(all_phase_data)):
        # Merge singletons and orphans
        singletons = singletons.union(all_phase_data[i][1])
        orphans = orphans.union(all_phase_data[i][2])

        # Identify breaks
        blocking = all_phase_data[i][0]
        j = 0
        for j in range(len(blocking) - 1):
            if blocking[j][-1] == blocking[j + 1][0] - 1:
                breaks.add(blocking[j][-1])

    # Split blocking in sample 1 to match other breaks
    blocks = []
    for b in all_phase_data[0][0]:
        br = np.argwhere([a in breaks or a in orphans for a in b[:-1]])
        blocks.extend([list(a) for a in np.split(b, br.flatten() + 1)])

    orphans = orphans.union([a[0] for a in blocks if len(a) == 1])
    blocks = [b for b in blocks if len(b) > 1]

    return blocks, singletons, orphans


def phase_blocks_sequential(df, blocksize=50e3, max_snps_per_block=10, alpha=0):
    if len(df) == 0:
        return []

    """NOTE: this currently only supports single-sample data"""
    assert len(df.SAMPLE.unique()) == 1
    assert len(df.CHR.unique()) == 1

    df = df.copy()
    df = df.reset_index(drop=True)

    # Identify contiguous segments of SNPs with phasing information
    segments = consecutive(np.where(~df.PHASE.isna())[0])

    orphans = set()
    if alpha < 1:
        ### Use binomial test to split adjacent SNPs with very different phased BAFs
        # (further dividing the segments identified above)
        # identify indices of adj SNPs that are significantly different
        df_merged = pd.concat([df, df.shift(-1).add_prefix('next_')], axis=1)
        # no_merge is the left index of each such pair
        no_merge = np.where(
            df_merged[:-1].apply(
                lambda x: binom_prop_test(
                    x.ALT,
                    x.REF,
                    x.FLIP,
                    x.next_ALT,
                    x.next_REF,
                    x.next_FLIP,
                    alpha=alpha,
                ),
                axis=1,
            )
        )[0]
        # stack the adjacent indices (top is left, bottom is right = left + 1)
        nm = np.vstack([no_merge, no_merge + 1])

        new_segments = []
        for seg in segments:
            # identify adjacent pairs that violate test (p<alpha) and are both in this segment
            violations = np.all(np.isin(nm, seg), axis=0)

            # map violations to indices in the segment
            keys = nm[1, np.nonzero(violations)][0]
            to_split = [np.where(seg == k)[0][0] for k in keys]
            splitted = np.split(seg, to_split)
            orphans = orphans.union([a[0] for a in splitted if len(a) == 1])
            new_segments.extend([a for a in splitted if len(a) > 1])

        segments = new_segments

    blocks = []

    for segment in segments:
        if len(segment) == 0:
            pass
        elif len(segment) == 1:
            blocks.append(list(segment))
        else:
            my_df = df.iloc[segment]
            ## Run blocking on my_df
            my_blocks = block_segment(my_df, blocksize, max_snps_per_block)
            mapping = {i: my_df.iloc[i].name for i in range(len(my_df))}
            blocks.extend([[mapping[a] for a in b] for b in my_blocks])

    blocks = sorted(blocks, key=lambda x: x[0])
    singletons = set(np.where(df.PHASE.isna())[0])

    return blocks, singletons, orphans


def collapse_blocks(df, blocks, singletons, orphans, ch):
    ### Construct blocked 1-bed table
    # First, merge blocks for those SNPs that are in the same block

    rows = []
    for sample, df0 in df.groupby('SAMPLE'):
        i = 0
        j = 0
        while i < len(df0):
            if i in singletons or i in orphans:
                r = df0.iloc[i]
                rows.append([ch, r.POS, r.POS, sample, r.ALT, r.REF, r.TOTAL, 1, [1]])
                i += 1
            else:
                block = blocks[j]
                assert i in block, (block, i)
                i += len(block)
                j += 1

                my_snps = df0.iloc[block]
                start = my_snps.POS.min()
                end = my_snps.POS.max()
                alt = np.sum(my_snps.FLIP * my_snps.REF + (1 - my_snps.FLIP) * my_snps.ALT).astype(np.uint64)
                ref = np.sum(my_snps.FLIP * my_snps.ALT + (1 - my_snps.FLIP) * my_snps.REF).astype(np.uint64)
                total = np.sum(my_snps.TOTAL)
                n_snps = len(my_snps)

                rows.append(
                    [
                        ch,
                        start,
                        end,
                        sample,
                        alt,
                        ref,
                        total,
                        n_snps,
                        list(my_snps.NOFLIP),
                    ]
                )

    return pd.DataFrame(
        rows,
        columns=[
            'CHR',
            'START',
            'END',
            'SAMPLE',
            'ALT',
            'REF',
            'TOTAL',
            '#SNPS',
            'HAPLO',
        ],
    ).sort_values(by=['CHR', 'START', 'SAMPLE'])


def block_segment(df, blocksize, max_snps_per_block):
    """
    Given a pandas dataframe containing read counts for a contiguous segment of SNPs,
    collects SNPs into phase blocks of size at most <blocksize> containing at most
    <max_snps_per_block> SNPs each.
    """
    splits = []
    df = df.copy().reset_index(drop=True)

    # walk along segment until either blocksize or max_snps_per_block is reached
    my_start = df.iloc[0].POS
    n = 1
    for i, r in df.iterrows():
        if i == 0:
            continue

        if n >= max_snps_per_block or r.POS - my_start > blocksize:
            splits.append(i)
            my_start = r.POS
            n = 1
        else:
            n += 1
    return np.split(np.arange(len(df)), splits)


def consecutive(data, stepsize=1):
    return np.split(data, np.where(np.diff(data) != stepsize)[0] + 1)


def binom_prop_test(alt1, ref1, flip1, alt2, ref2, flip2, alpha=0.1):
    """
    Returns True if there is sufficient evidence that SNPs 1 and 2 should not be merged, False otherwise.
    """
    if np.isnan(flip1) or np.isnan(flip2):
        return True

    x1 = ref1 if flip1 > 0 else alt1
    x2 = ref2 if flip2 > 0 else alt2
    n1 = max(1, alt1 + ref1)
    n2 = max(1, alt2 + ref2)

    p1 = x1 / n1
    p2 = x2 / n2

    p = (x1 + x2) / (n1 + n2)
    if p == 0:
        p = 0.001
    elif p == 1:
        p = 0.999
    denom = np.sqrt(p * (1 - p) * (1 / n1 + 1 / n2))
    z = (p1 - p2) / denom

    if z > 0:
        return 2 * (1 - norm.cdf(z)) < alpha
    else:
        return 2 * norm.cdf(z) < alpha
