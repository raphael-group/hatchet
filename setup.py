import os
import sys
from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension, build_ext

_cpp_dir = "src/hatchet/cluster_bins/hmm/_hmm_cpp"

# Compiler flags
compile_args = ["-O3", "-march=native", "-ffast-math"]
link_args = []
include_dirs = []
library_dirs = []

if sys.platform == "darwin":
    # macOS: libomp is provided by llvm-openmp (conda-forge) or homebrew libomp.
    # Use -Xpreprocessor flag required by Apple clang; link -lomp.
    conda_prefix = os.environ.get("CONDA_PREFIX", "")
    if conda_prefix:
        include_dirs = [os.path.join(conda_prefix, "include")]
        library_dirs = [os.path.join(conda_prefix, "lib")]
    compile_args += ["-Xpreprocessor", "-fopenmp"]
    link_args += ["-lomp"]
else:
    compile_args += ["-fopenmp"]
    link_args += ["-fopenmp"]

ext_modules = [
    Pybind11Extension(
        "hatchet.cluster_bins.hmm._hmm_cpp",
        sources=[
            f"{_cpp_dir}/fwd_bwd.cpp",
            f"{_cpp_dir}/m_steps.cpp",
            f"{_cpp_dir}/loglik.cpp",
            f"{_cpp_dir}/run_hmm.cpp",
            f"{_cpp_dir}/bindings.cpp",
        ],
        extra_compile_args=compile_args,
        extra_link_args=link_args,
        include_dirs=include_dirs,
        library_dirs=library_dirs,
    )
]

setup(ext_modules=ext_modules, cmdclass={"build_ext": build_ext})
