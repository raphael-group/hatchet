import logging

try:
    from hatchet.cluster_bins.hmm._hmm_cpp import (
        omp_get_max_threads as _cpp_omp_get_max_threads,
        run_hmm as _cpp_run_hmm,
    )

    _USE_CPP = True
except ImportError:
    _cpp_omp_get_max_threads = None
    _cpp_run_hmm = None
    _USE_CPP = False

if _USE_CPP:
    logging.info(
        f"HMM backend: C++ (_hmm_cpp), omp_threads={_cpp_omp_get_max_threads()}"
    )
else:
    import numba as _numba

    logging.info(
        f"HMM backend: Numba/scipy (C++ extension not found), numba_threads={_numba.get_num_threads()}"
    )
