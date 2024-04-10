import numpy as np
import hatchet.utils.Supporting as sp


def correct_haplotypes(
    orig_bafs,
    min_prop_switch=0.01,
    n_segments=10,
    min_switch_density=0.1,
    min_mean_baf=0.45,
    minmax_al_imb=0.02,
):
    # Count switches using only samples with mean allelic imbalance above <minmax_al_imb>
    imb_samples = np.where(np.mean(np.abs(orig_bafs - 0.5), axis=0) > minmax_al_imb)[0]

    if len(imb_samples) == 0:
        sp.log(
            msg=f'No sample with avg. allelic imbalance above [{minmax_al_imb}], skipping correction.\n',
            level='INFO',
        )
        return orig_bafs, None

    bafs = orig_bafs[:, imb_samples]

    # Look for haplotype switches
    above_mid = bafs > 0.5
    is_alternating = np.concatenate(
        [
            np.zeros((1, bafs.shape[1]), dtype=bool),
            np.logical_xor(above_mid[1:], above_mid[:-1]),
        ]
    )
    haplotype_switches = np.where(np.all(is_alternating, axis=1))[0]
    prop_switched = len(haplotype_switches) / len(is_alternating)
    sp.log(
        msg=f'Found haplotype switches in {prop_switched *100 : .2f}% of bins\n',
        level='INFO',
    )

    if prop_switched > min_prop_switch:
        # If sufficient switches are found, run segmentation to identify segments w/ many switches

        sp.log(
            msg=f'Checking haplotype switching using [{n_segments}] segments\n',
            level='INFO',
        )
        # Segment using the mean BAFs only - faster and more reliable than using all samples
        A, bp = segmented_piecewise(np.mean(bafs, axis=1).reshape(-1, 1), pieces=n_segments)

        ts = backtrack(bp[:, :n_segments])

        for idx in np.where(np.diff(ts) == 0)[0]:
            del ts[idx]

        segments = [orig_bafs[ts[i] : ts[i + 1]] for i in range(len(ts) - 1)]

        # Identify problematic segments as those with many switches and mean near 0.5
        # (note that mean(BAF_i) across samples i is always <= 0.5 by def. from EM function)
        switch_densities = np.array(
            [
                len(haplotype_switches[(haplotype_switches >= ts[i]) & (haplotype_switches < ts[i + 1])])
                / (ts[i + 1] - ts[i])
                for i in range(len(ts) - 1)
            ]
        )
        segment_means = np.array([(np.mean(s) if len(s) > 0 else 0.5) for s in segments])

        # ALSO only correct segments with allelic imbalance at least <min_al_imb> in at least 1 sample
        segment_imbalances = np.array(
            [(np.max(np.abs(0.5 - np.mean(np.minimum(s, 1 - s), axis=0))) if len(s) > 0 else 0) for s in segments]
        )

        """
        segment_lengths = [len(s) for s in segments]
        [sp.log(msg=f'Segment {i}: length {a},\tmean {b:.3f},\timbalance {c:.3f},\tswitch prop. {d:.3f}\n',
                level = 'INFO')
               for i, (a,b,c, d)
               in enumerate(zip(segment_lengths, segment_means, segment_imbalances, switch_densities))]
        """

        bad_segments = np.where(
            np.logical_and(
                np.logical_and(
                    switch_densities >= min_switch_density,
                    segment_means >= min_mean_baf,
                ),
                segment_imbalances >= minmax_al_imb,
            )
        )[0]

        sp.log(
            msg=f'Identified {len(bad_segments)} segments with haplotype switching\n',
            level='INFO',
        )
        final_segments = []

        for s_idx in range(len(segments)):
            if s_idx in bad_segments:
                seg = segments[s_idx].copy()
                assert len(seg) > 0, 'Haplotype switch correction flagged a length-0 segment'

                # Identify the sample with the most extreme allelic imbalance for this segment
                minseg = np.minimum(seg, 1 - seg)
                mmeans = np.mean(minseg, axis=0)
                extreme_sample = np.argmin(mmeans)

                # Flip all of those bins with mhBAF > 0.5 in the corresponding sample
                my_flips = np.where(seg[:, extreme_sample] > 0.5)

                seg[my_flips] = 1 - seg[my_flips]

                my_mean = np.mean(seg)

                # If this ends up choosing the haplotype that is more abundant on average,
                if my_mean > 0.5:
                    # flip all bins to select the minor haplotype instead
                    seg = 1 - seg

                final_segments.append(seg)
            else:
                final_segments.append(segments[s_idx])

        return np.concatenate(final_segments), ts
    else:
        sp.log(
            msg=f'Insufficient haplotype switching detected (<[{min_prop_switch}]), skipping correction.\n',
            level='INFO',
        )
        return orig_bafs, None


def backtrack(bp):
    n, p = bp.shape

    starts = [bp[-1, -1]]
    for i in range(p - 1, 0, -1):
        starts.append(bp[starts[-1] - 1, i - 1])

    thresholds = starts[::-1] + [n]
    return thresholds


def segmented_piecewise(X, pieces=2):
    n, s = X.shape
    segcost_memo = {}

    def segment_cost(i, j):
        if (i, j) in segcost_memo:
            return segcost_memo[i, j]
        else:
            my_mean = np.mean(X[i:j], axis=0) if j > i else 0
            result = np.sum(np.square(X[i:j] - my_mean))
            segcost_memo[i, j] = result
            return result

    # penalty A[i, p] is the minimum error possible for fitting X[0] with p+1 pieces
    A = np.zeros((n, pieces))
    backpoint = np.zeros((n, pieces), dtype=int)

    A[:, 0] = [segment_cost(0, i) for i in range(n)]

    for p in range(1, pieces):
        for t in range(1, n):
            # search over t' < t
            best_cost = np.inf
            best_tprime = None
            for tprime in range(t + 1):
                # compute cost for making a new segment from tprime through t
                cost = A[tprime, p - 1] + segment_cost(tprime, t)
                if cost < best_cost:
                    best_cost = cost
                    best_tprime = tprime

            A[t, p] = best_cost
            backpoint[t, p] = best_tprime   # if this throws a TypeError for int(None), then X may have NaNs
    return A, backpoint
