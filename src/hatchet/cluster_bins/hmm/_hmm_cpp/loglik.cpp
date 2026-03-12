/**
 * Log-likelihood kernel for the 2-mixture BAF+RDR HMM.
 *
 * BAF: Beta-Binomial (lgamma-based log betaln).
 * RDR: Gaussian (closed-form log-normal constant).
 * Parallelised over bins with OpenMP; per-bin (k, m) loop is sequential.
 */

#include "loglik.h"

#include <cmath>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

static const double LOG2PI = std::log(2.0 * M_PI);

static inline double betaln_cpp(double a, double b) {
    return std::lgamma(a) + std::lgamma(b) - std::lgamma(a + b);
}


void compute_loglik_cpp(
    const double* X_rdrs,
    const double* X_alphas,
    const double* X_betas,
    const double* X_totals,
    const double* rdr_means,
    const double* rdr_vars,
    const double* baf_means,
    const double* baf_taus,
    double*       lls0,
    double*       lls1,
    int N, int K, int M)
{
    // Precompute per-(k,m) constants: bb_alpha, bb_beta, bb_delta, log_norm_const.
    // These are small (K*M) and shared across all bins.
    std::vector<double> bb_alpha(K * M);
    std::vector<double> bb_beta(K * M);
    std::vector<double> bb_delta(K * M);
    std::vector<double> log_norm_const(K, 0.0);

    for (int k = 0; k < K; ++k) {
        for (int m = 0; m < M; ++m) {
            double tau  = baf_taus[m];
            double mean = baf_means[(long)k * M + m];
            double a    = tau * mean;
            double b    = tau * (1.0 - mean);
            bb_alpha[(long)k * M + m] = a;
            bb_beta[(long)k * M + m]  = b;
            bb_delta[(long)k * M + m] = betaln_cpp(a, b);
            log_norm_const[k] += 0.5 * (LOG2PI + std::log(rdr_vars[(long)k * M + m]));
        }
    }

    // Precompute per-(n,m) k-invariant quantities: log_bc and lgamma(total+tau).
    // log_bc_nm[n*M+m]           = lgamma(total+1) - lgamma(alpha+1) - lgamma(beta+1)
    // lgamma_total_tau_nm[n*M+m] = lgamma(total + tau_m)  [shared denominator in betaln]
    std::vector<double> log_bc_nm(N * M);
    std::vector<double> lgamma_total_tau_nm(N * M);
    // baf_taus indexed by m only — precompute once per call (cheap: N*M lgammas)
    for (int n = 0; n < N; ++n) {
        for (int m = 0; m < M; ++m) {
            double alpha_nm = X_alphas[(long)n * M + m];
            double beta_nm  = X_betas[(long)n * M + m];
            double total_nm = X_totals[(long)n * M + m];
            log_bc_nm[(long)n * M + m] = (std::lgamma(total_nm + 1.0)
                                          - std::lgamma(alpha_nm + 1.0)
                                          - std::lgamma(beta_nm  + 1.0));
            lgamma_total_tau_nm[(long)n * M + m] = std::lgamma(total_nm + baf_taus[m]);
        }
    }

    // Main loop: parallelise over bins.
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int n = 0; n < N; ++n) {
        for (int k = 0; k < K; ++k) {
            double ll_rdr    = -log_norm_const[k];
            double ll_baf_h0 = 0.0;
            double ll_baf_h1 = 0.0;

            for (int m = 0; m < M; ++m) {
                double alpha_nm    = X_alphas[(long)n * M + m];
                double beta_nm     = X_betas[(long)n * M + m];
                double rdr_nm      = X_rdrs[(long)n * M + m];
                double log_bc      = log_bc_nm[(long)n * M + m];
                double lg_tot_tau  = lgamma_total_tau_nm[(long)n * M + m];

                // Gaussian RDR contribution for sample m
                double diff = rdr_nm - rdr_means[(long)k * M + m];
                double var  = rdr_vars[(long)k * M + m];
                ll_rdr -= 0.5 * diff * diff / var;

                // Beta-Binomial BAF contribution for sample m
                // betaln(x+a, y+b) = lgamma(x+a) + lgamma(y+b) - lgamma(x+y+tau)
                // lgamma(x+y+tau) = lg_tot_tau (k-invariant, precomputed above)
                double a     = bb_alpha[(long)k * M + m];
                double b     = bb_beta[(long)k * M + m];
                double delta = bb_delta[(long)k * M + m];

                ll_baf_h0 += log_bc + std::lgamma(alpha_nm + a) + std::lgamma(beta_nm  + b) - lg_tot_tau - delta;
                ll_baf_h1 += log_bc + std::lgamma(beta_nm  + a) + std::lgamma(alpha_nm + b) - lg_tot_tau - delta;
            }

            lls0[(long)n * K + k] = ll_rdr + ll_baf_h0;
            lls1[(long)n * K + k] = ll_rdr + ll_baf_h1;
        }
    }
}
