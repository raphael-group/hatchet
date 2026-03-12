/**
 * pybind11 module exposing run_hmm (full EM loop) and omp_get_max_threads.
 *
 * All inputs must be C-contiguous float64 (or int64 for X_lengths).
 * forcecast ensures automatic conversion if needed.
 */

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "run_hmm.h"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace py = pybind11;

using f64arr = py::array_t<double,  py::array::c_style | py::array::forcecast>;
using i64arr = py::array_t<int64_t, py::array::c_style | py::array::forcecast>;


static py::dict run_hmm_py(
    int    K,
    f64arr X_rdrs,
    f64arr X_alphas,
    f64arr X_betas,
    f64arr X_totals,
    i64arr X_lengths,
    f64arr log_switchprobs,
    f64arr log_stayprobs,
    f64arr log_transmat,
    f64arr rdr_means0,
    f64arr rdr_vars0,
    f64arr baf_means0,
    f64arr baf_taus0,
    int    n_iter,
    double min_covar,
    double tol_ll,
    double tol,
    int    tau_iters,
    double min_tau,
    double max_tau,
    double baf_eps,
    double ig_alpha,
    f64arr ig_beta_arr)
{
    int N = (int)X_rdrs.shape(0);
    int M = (int)X_rdrs.shape(1);
    int S = (int)X_lengths.shape(0);

    RunHMMResult res = run_hmm_cpp(
        K, N, M, S,
        X_rdrs.data(),
        X_alphas.data(),
        X_betas.data(),
        X_totals.data(),
        X_lengths.data(),
        log_switchprobs.data(),
        log_stayprobs.data(),
        log_transmat.data(),
        rdr_means0.data(),
        rdr_vars0.data(),
        baf_means0.data(),
        baf_taus0.data(),
        n_iter, min_covar, tol_ll, tol,
        tau_iters, min_tau, max_tau, baf_eps,
        ig_alpha, ig_beta_arr.data());

    // Helper: copy a flat vector into a shaped numpy array.
    auto make_arr = [](const std::vector<double>& v,
                       std::vector<py::ssize_t> shape) -> py::array_t<double> {
        auto arr = py::array_t<double>(shape);
        std::copy(v.begin(), v.end(), arr.mutable_data());
        return arr;
    };

    py::dict d;
    d["RDR_means"]      = make_arr(res.rdr_means,      {K, M});
    d["RDR_vars"]       = make_arr(res.rdr_vars,        {K, M});
    d["BAF_means"]      = make_arr(res.baf_means,       {K, M});
    d["BAF_taus"]       = make_arr(res.baf_taus,        {M});
    d["log_startprobs"] = make_arr(res.log_startprobs,  {K, 2});
    d["full_posts"]     = make_arr(res.posts,            {N, K, 2});
    d["lls0"]           = make_arr(res.lls0,             {N, K});
    d["lls1"]           = make_arr(res.lls1,             {N, K});
    d["elbo_trace"]     = py::cast(res.elbo_trace);
    d["model_ll"]       = py::cast(res.loglik);
    d["n_iters_done"]   = py::cast(res.n_iters_done);
    return d;
}


PYBIND11_MODULE(_hmm_cpp, m) {
    m.doc() = "C++ HMM kernel: full EM loop (run_hmm)";
    m.def(
        "run_hmm",
        &run_hmm_py,
        py::arg("K"),
        py::arg("X_rdrs"),
        py::arg("X_alphas"),
        py::arg("X_betas"),
        py::arg("X_totals"),
        py::arg("X_lengths"),
        py::arg("log_switchprobs"),
        py::arg("log_stayprobs"),
        py::arg("log_transmat"),
        py::arg("rdr_means0"),
        py::arg("rdr_vars0"),
        py::arg("baf_means0"),
        py::arg("baf_taus0"),
        py::arg("n_iter")    = 10,
        py::arg("min_covar") = 1e-3,
        py::arg("tol_ll")    = 1e-4,
        py::arg("tol")       = 1e-6,
        py::arg("tau_iters") = 1,
        py::arg("min_tau")   = 50.0,
        py::arg("max_tau")   = 100.0,
        py::arg("baf_eps")   = 1e-6,
        py::arg("ig_alpha")  = 10.0,
        py::arg("ig_beta"),
        R"doc(
Full C++ EM loop for the 2-mixture BAF+RDR HMM.

Parameters
----------
K : int
    Number of cluster states.
X_rdrs : (N, M) float64
X_alphas : (N, M) float64
X_betas : (N, M) float64
X_totals : (N, M) float64
X_lengths : (S,) int64  segment lengths
log_switchprobs : (N,) float64
log_stayprobs   : (N,) float64
log_transmat  : (K, K) float64
rdr_means0 : (K, M) float64  initial RDR means
rdr_vars0  : (K, M) float64  initial RDR variances
baf_means0 : (K, M) float64  initial BAF means
baf_taus0  : (M,) float64    initial BB dispersion
n_iter     : int   max EM iterations
min_covar  : float min RDR variance floor
tol_ll     : float convergence threshold on log-likelihood delta
tol        : float min effective cluster size / start-prob floor
tau_iters  : int   iterations during which tau is updated
min_tau, max_tau : float  Brent search bounds for tau
baf_eps    : float  BAF mean search bounds [eps, 1-eps]

Returns
-------
dict with keys: RDR_means, RDR_vars, BAF_means, BAF_taus, log_startprobs,
                full_posts (N,K,2), lls0 (N,K), lls1 (N,K),
                elbo_trace (list), model_ll (float), n_iters_done (int)
)doc"
    );
    m.def("omp_get_max_threads", []() -> int {
#ifdef _OPENMP
        return omp_get_max_threads();
#else
        return 1;
#endif
    }, "Return the maximum number of OpenMP threads (1 if not compiled with OpenMP).");
}
