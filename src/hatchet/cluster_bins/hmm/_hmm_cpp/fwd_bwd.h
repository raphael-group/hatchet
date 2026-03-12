#pragma once
#include <cstdint>

/**
 * Parallel forward-backward over all segments.
 *
 * All arrays are C-contiguous (row-major) float64.
 *
 * lls0         (N, K)   log-likelihoods under haplotype h=0
 * lls1         (N, K)   log-likelihoods under haplotype h=1
 * seg_starts   (S,)     start bin index for each segment
 * seg_lengths  (S,)     number of bins in each segment
 * log_startprobs (K, 2) log start probabilities
 * log_switchprobs (N,)  log phase-switch probability per bin
 * log_stayprobs   (N,)  log phase-stay   probability per bin
 * log_transmat  (K, K)  log cluster transition matrix
 * posts        (N, K, 2) output posteriors — written in place
 *
 * Returns total log-likelihood (sum over segments).
 */
double forward_backward_cpp(
    const double*   lls0,
    const double*   lls1,
    const int64_t*  seg_starts,
    const int64_t*  seg_lengths,
    const double*   log_startprobs,
    const double*   log_switchprobs,
    const double*   log_stayprobs,
    const double*   log_transmat,
    double*         posts,
    int N, int K, int S
);
