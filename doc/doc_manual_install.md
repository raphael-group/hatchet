### Manual Installation

If you wish to install `HATCHet` directly from this repository, the steps are a bit more involved.
The core module of HATCHet is written in C++11 and thus requires a modern C++ compiler (GCC >= 4.8.1, or Clang).
As long as you have a recent version of GCC or Clang installed, `setuptools` should automatically be able to download a
recent version of `cmake` and compile the Hatchet code into a working package.

The installation process can be broken down into the following steps:

1. **Get [Gurobi](http://www.gurobi.com/)** (>= 6.0)

    The coordinate-method applied by HATCHet is based on several integer linear programming (ILP) formulations. Gurobi is a commercial ILP solver with two licensing options: (1) a single-host license where the license is tied to a single computer and (2) a network license for use in a compute cluster (using a license server in the cluster). Both options are freely and [easily available](http://www.gurobi.com/academia/academia-center) for users in academia.
[Download](https://www.gurobi.com/downloads/gurobi-optimizer-eula) Gurobi for your specific platform.

2. **Set GUROBI_HOME environment variable**
    ```shell
    $ export GUROBI_HOME=/path/to/gurobiXXX
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

    When the compilation process fails or when the environment has special requirements, you may have to manually specify the required paths to Gurobi by following the [detailed intructions](doc/doc_compilation.md).

7. **Install required utilities**

    For reading BAM files, read counting, allele counting, and SNP calling, you need to install [SAMtools and BCFtools](http://www.htslib.org/doc/).

