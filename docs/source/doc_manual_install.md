### Manual Installation

If you wish to install `HATCHet` directly from this repository, the steps are a bit more involved.

Note that the complexity of manual installation is largely because the `compute-cn` step (determination of
allele-specific copy numbers) of the HATCHet pipeline uses custom-written C++11 code that uses the
[Gurobi](http://www.gurobi.com/) optimizer. If you do not have a valid Gurobi license (though it is
[easily available](http://www.gurobi.com/academia/academia-center) for users in academia), then the C++ parts of
HATCHet do not necessarily need to be compiled, and you can read the
[Compiling HATCHet without the built-in Gurobi optimizer](#withoutgurobi) section of this page.

#### Compiling HATCHet with the built-in Gurobi optimizer

The core optimization module of HATCHet is written in C++11 and thus requires a modern C++ compiler (GCC >= 4.8.1, or Clang).
As long as you have a recent version of GCC or Clang installed, `setuptools` should automatically be able to download a
recent version of `cmake` and compile the Hatchet code into a working package.

The installation process can be broken down into the following steps:

1. **Get [Gurobi](http://www.gurobi.com/)** (v9.0.2)

    The coordinate-method applied by HATCHet is based on several integer linear programming (ILP) formulations. Gurobi is a commercial ILP solver with two licensing options: (1) a single-host license where the license is tied to a single computer and (2) a network license for use in a compute cluster (using a license server in the cluster). Both options are freely and [easily available](http://www.gurobi.com/academia/academia-center) for users in academia.
[Download](https://www.gurobi.com/downloads/gurobi-optimizer-eula) Gurobi for your specific platform.


2. **Set GUROBI_HOME environment variable**
    ```shell
    $ export GUROBI_HOME=/path/to/gurobi902
    ```
    Set `GUROBI_HOME` to where you download Gurobi. Here `XXX` is the 3-digit version of gurobi.


3. **Build Gurobi**
    ```shell
    $ cd "${GUROBI_HOME}"
    $ cd linux64/src/build/
    $ make
    $ cp libgurobi_c++.a ../../lib
    ```
    Substitute `mac64` for `linux64` if using the Mac OSX platform.


4. **Create a new venv/conda environment for Hatchet**

    `Hatchet` is a Python 3 package. Unless you want to compile/install it in your default Python 3 environment, you will
want to create either a new Conda environment for Python 3 and activate it:
    ```
    conda create --name hatchet python=3.8
    conda activate hatchet
    ```
    or use `virtualenv` through `pip`:
    ```
    python3 -m pip virtualenv env
    source env/bin/activate
    ```


5. **Install basic packages**

    It is **highly recommended** that you upgrade your `pip` and `setuptools` versions to the latest, using:
    ```shell
    pip install -U pip
    pip install -U setuptools
    ```


6. **Build and install HATCHet**

    Execute the following commands from the root of HATCHet's repository.
    ```shell
    $ pip install .
    ```

    **NOTE**: If you experience a failure of compilation with an error message like:
    ```
    _undefined reference to symbol 'pthread_create@@GLIBC_2.2.5'_.
    ```

    you may need to set `CXXFLAGS` to `-pthread` before invoking the command:
    ```shell
    $ CXXFLAGS=-pthread pip install .
    ```

    When the compilation process fails or when the environment has special requirements, you may have to manually specify the required paths to Gurobi by following the [detailed instructions](doc_compilation.md).


7. **Install required utilities**

    For reading BAM files, read counting, allele counting, and SNP calling, you need to install [SAMtools, BCFtools, and tabix](http://www.htslib.org/doc/) as well as [mosdepth](https://github.com/brentp/mosdepth).
    If you want to perform reference-based phasing, you must also install [shapeit](https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html), [picard](https://broadinstitute.github.io/picard/), and [bgzip](http://www.htslib.org/doc/). The easiest way to install these is via `conda`, as all are available from the `bioconda` channel (except `shapeit` which is available from channel `dranew`).


#### Compiling HATCHet *without* the built-in Gurobi optimizer
<a name="withoutgurobi"></a>

If you wish to use an alternate ILP optimizer, then you do not need a C++ compiler.

In this case, set the `HATCHET_BUILD_NOEXT` environment variable to `1` (using `export HATCHET_BUILD_NOEXT=1`),
set the environment variable `HATCHET_COMPUTE_CN_SOLVER` to a Pyomo-supported solver (see the
[Using a different Pyomo-supported solver](README.md#usingasolver_other) section of the README for more details)
and proceed directly to step (4) above.
