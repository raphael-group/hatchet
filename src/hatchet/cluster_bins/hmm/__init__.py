import os
import logging

try:
    from hatchet.cluster_bins.hmm._hmm_cpp import (
        omp_get_max_threads as _cpp_omp_get_max_threads,
        run_hmm as _cpp_run_hmm,
    )

    _CPP_IMPORTABLE = True
except ImportError:
    _cpp_omp_get_max_threads = None
    _cpp_run_hmm = None
    _CPP_IMPORTABLE = False


def _cpp_disabled_by_env():
    """True if HATCHET_DISABLE_CPP forces the Numba/scipy fallback."""
    return os.environ.get("HATCHET_DISABLE_CPP", "0").lower() in ("1", "true", "yes")


# Runtime-overridable so tests can force either backend via monkeypatch.
_USE_CPP = _CPP_IMPORTABLE and not _cpp_disabled_by_env()

if _USE_CPP:
    logging.info(
        f"HMM backend: C++ (_hmm_cpp), omp_threads={_cpp_omp_get_max_threads()}"
    )
else:
    import numba as _numba

    logging.info(
        f"HMM backend: Numba/scipy (C++ extension not found), numba_threads={_numba.get_num_threads()}"
    )
