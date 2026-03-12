/**
 * BAF M-step: MLE for Beta-Binomial means via Boost.Math Brent minimization.
 *
 * For each (k, m) pair, maximizes the posterior-weighted Beta-Binomial
 * log-likelihood over p in (eps, 1-eps):
 *
 *   Q(p) = C + Σ_n [ w0[n]*(lgamma(α_n+τp) + lgamma(β_n+τ(1-p)))
 *                   + w1[n]*(lgamma(β_n+τp) + lgamma(α_n+τ(1-p))) ]
 *             - Wk * (lgamma(τp) + lgamma(τ(1-p)))
 *
 * where C = Σ_n (w0+w1)*(-lgamma(α_n+β_n+τ)) + Wk*lgamma(τ) is a constant
 * w.r.t. p precomputed once per (k,m) outside Brent.
 *
 * The outer (k, m) loop is parallelised with OpenMP collapse(2).
 */

#include "m_steps.h"

#include <cmath>
#include <vector>
#include <boost/math/tools/minima.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif


void update_baf_means_cpp(
    const double* alphas_mn,   // (M, N)
    const double* betas_mn,    // (M, N)
    const double* posts_kn2,   // (K, N, 2)
    const double* baf_taus,    // (M,)
    double*       p_km,        // (K, M)
    int N, int K, int M, double eps)
{
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(dynamic)
#endif
    for (int k = 0; k < K; ++k) {
        for (int m = 0; m < M; ++m) {
            const double* alpha_m = alphas_mn + (long)m * N;
            const double* beta_m  = betas_mn  + (long)m * N;
            // posts_kn2[k, n, h] = posts_kn2[k*N*2 + n*2 + h]
            const double* pkn2    = posts_kn2 + (long)k * N * 2;
            double tau = baf_taus[m];
            double lgamma_tau = std::lgamma(tau);

            // Precompute per-n constants and aggregate weights
            // lgt[n] = lgamma(alpha_n + beta_n + tau)  [constant w.r.t. p]
            // Wk = Σ_n (w0[n] + w1[n])
            // C  = Σ_n (w0[n]+w1[n]) * (-lgt[n]) + Wk * lgamma(tau)
            std::vector<double> lgt(N);
            double Wk = 0.0;
            double C  = 0.0;

            for (int n = 0; n < N; ++n) {
                double w0n = pkn2[n * 2 + 0];
                double w1n = pkn2[n * 2 + 1];
                double wn  = w0n + w1n;
                lgt[n] = std::lgamma(alpha_m[n] + beta_m[n] + tau);
                Wk += wn;
                C  -= wn * lgt[n];
            }
            C += Wk * lgamma_tau;

            // neg_Q(p): with constants hoisted, inner loop needs 4 lgamma per n
            auto neg_Q = [&](double p) -> double {
                double a = tau * p;
                double b = tau * (1.0 - p);
                // betaln(a,b) normalizer contribution (no n-dependence)
                double norm_term = Wk * (std::lgamma(a) + std::lgamma(b));
                double Q = C - norm_term;
                for (int n = 0; n < N; ++n) {
                    double alpha_n = alpha_m[n];
                    double beta_n  = beta_m[n];
                    double w0n = pkn2[n * 2 + 0];
                    double w1n = pkn2[n * 2 + 1];
                    // h=0 contribution: lgamma(α_n+a) + lgamma(β_n+b) - lgt[n]
                    // h=1 contribution: lgamma(β_n+a) + lgamma(α_n+b) - lgt[n]
                    // lgt[n] absorbed into C above
                    Q += w0n * (std::lgamma(alpha_n + a) + std::lgamma(beta_n + b))
                       + w1n * (std::lgamma(beta_n  + a) + std::lgamma(alpha_n + b));
                }
                return -Q;
            };

            auto result = boost::math::tools::brent_find_minima(neg_Q, eps, 1.0 - eps, 32);
            p_km[(long)k * M + m] = result.first;
        }
    }
}


// ---------- RDR M-step ----------

void update_rdr_params_cpp(
    const double* X_rdrs,
    const double* posts_nk,
    double*       rdr_means,
    double*       rdr_vars,
    int N, int K, int M, double min_covar,
    double ig_alpha,
    const double* ig_beta)
{
    // n-outer loop: posts_nk[n*K + k] reads are sequential across k,
    // avoiding the strided-by-K cache misses of the old k-outer layout.
    // Thread-local (K*M) accumulators are reduced at the end.
    std::vector<double> global_sum1(K * M, 0.0);
    std::vector<double> global_sum2(K * M, 0.0);
    std::vector<double> global_Nk(K, 0.0);

#ifdef _OPENMP
#pragma omp parallel
    {
        std::vector<double> loc1(K * M, 0.0);
        std::vector<double> loc2(K * M, 0.0);
        std::vector<double> locNk(K, 0.0);

#pragma omp for nowait schedule(static)
        for (int n = 0; n < N; ++n) {
            for (int k = 0; k < K; ++k) {
                double p = posts_nk[(long)n * K + k];
                locNk[k] += p;
                for (int m = 0; m < M; ++m) {
                    double r = X_rdrs[(long)n * M + m];
                    loc1[(long)k * M + m] += p * r;
                    loc2[(long)k * M + m] += p * r * r;
                }
            }
        }

#pragma omp critical
        {
            for (int km = 0; km < K * M; ++km) {
                global_sum1[km] += loc1[km];
                global_sum2[km] += loc2[km];
            }
            for (int k = 0; k < K; ++k)
                global_Nk[k] += locNk[k];
        }
    }
#else
    for (int n = 0; n < N; ++n) {
        for (int k = 0; k < K; ++k) {
            double p = posts_nk[(long)n * K + k];
            global_Nk[k] += p;
            for (int m = 0; m < M; ++m) {
                double r = X_rdrs[(long)n * M + m];
                global_sum1[(long)k * M + m] += p * r;
                global_sum2[(long)k * M + m] += p * r * r;
            }
        }
    }
#endif

    bool use_ig = (ig_alpha > 0 && ig_beta != nullptr);
    for (int k = 0; k < K; ++k) {
        double nk = std::max(global_Nk[k], 1e-10);
        for (int m = 0; m < M; ++m) {
            double mean = global_sum1[(long)k * M + m] / nk;
            double weighted_var = global_sum2[(long)k * M + m] / nk - mean * mean;
            rdr_means[(long)k * M + m] = mean;
            if (use_ig) {
                double raw_SS = weighted_var * nk;
                double beta_m = ig_beta[m];
                rdr_vars[(long)k * M + m] = std::max(
                    (raw_SS + 2.0 * beta_m) / (nk + 2.0 * (ig_alpha + 1.0)),
                    min_covar);
            } else {
                rdr_vars[(long)k * M + m] = std::max(weighted_var, min_covar);
            }
        }
    }
}


// ---------- BAF tau M-step ----------

void update_baf_tau_cpp(
    const double* alphas_nm,
    const double* betas_nm,
    const double* posts_nk,
    double*       baf_taus,
    int N, int K, int M,
    double p_fixed,
    double min_tau, double max_tau)
{
    // Find bins where argmax_k posts_nk[n, k] == 0  (cluster 0 = diploid anchor)
    std::vector<int> mask0;
    mask0.reserve(N);
    for (int n = 0; n < N; ++n) {
        int argmax    = 0;
        double maxval = posts_nk[(long)n * K + 0];
        for (int k = 1; k < K; ++k) {
            double v = posts_nk[(long)n * K + k];
            if (v > maxval) { maxval = v; argmax = k; }
        }
        if (argmax == 0) mask0.push_back(n);
    }

    int N0 = (int)mask0.size();
    if (N0 < 2) return;  // not enough bins

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (int m = 0; m < M; ++m) {
        // Gather alpha/beta for cluster-0 bins
        std::vector<double> a_vals(N0);
        std::vector<double> b_vals(N0);
        for (int i = 0; i < N0; ++i) {
            int n      = mask0[i];
            a_vals[i]  = alphas_nm[(long)n * M + m];
            b_vals[i]  = betas_nm[(long)n * M + m];
        }

        double p = p_fixed;

        // Optimise neg log-likelihood in log(tau) space — matches Python exactly.
        auto neg_ll_logtau = [&](double log_tau) -> double {
            double tau    = std::exp(log_tau);
            double a0     = tau * p;
            double b0     = tau * (1.0 - p);
            double lg_a0  = std::lgamma(a0);
            double lg_b0  = std::lgamma(b0);
            double lg_ab0 = std::lgamma(a0 + b0);
            double ll = 0.0;
            for (int i = 0; i < N0; ++i) {
                double a1 = a_vals[i] + a0;
                double b1 = b_vals[i] + b0;
                ll += (std::lgamma(a1) + std::lgamma(b1)
                       - std::lgamma(a1 + b1)
                       - lg_a0 - lg_b0 + lg_ab0);
            }
            return -ll;
        };

        auto result = boost::math::tools::brent_find_minima(
            neg_ll_logtau, std::log(min_tau), std::log(max_tau), 32);
        baf_taus[m] = std::exp(result.first);
    }
}
