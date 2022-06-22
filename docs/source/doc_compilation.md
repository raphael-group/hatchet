# Detailed compilation

To perform the compilation, execute the following commands from the root of HATCHet's repository.

```shell
$ mkdir build
$ cd build/
$ ccmake ..
$ make
```

HATCHet's compilation process attempts to automatically find the following Gurobi's paths.

| Name | Path | Comment |
|------|------|---------|
| `GUROBI_CPP_LIB` | `/to/gurobiXXX/YY/lib/libgurobi_c++.a`  | <ul><li>`/to/gurobi` is the path to Gurobi's home, typically `/opt/gurobiXXX` for linux and `/Library/gurobiXXX` for mac</li><li>`XXX` is the Gurobi full version, e.g. 702 or 751</li><li>`YY` depends on os, typically `linux64` for linux or `mac64` for mac</li></ul> |
| `GUROBI_INCLUDE_DIR` | `/to/gurobiXXX/YY/include`  | <ul><li>`/to/gurobi` is the path to Gurobi's home</li><li>`XXX` is the Gurobi full version</li><li>`YY` depends on os</li></ul> |
| `GUROBI_LIB` | `/to/gurobiXXX/YY/lib/libgurobiZZ.QQ`  | <ul><li>`/to/gurobiXXX` is the path to Gurobi's home</li><li>`XXX` is the Gurobi full version</li><li>`YY` depends on os</li><li>`ZZ` are typically the first 2 numbers of `XXX`</li><li>`QQ` is typically `so` but becomes `dylib` for MAC version since version 8.10</li></ul> |

If the automatic compilation fails to find the Gurobi's paths, these need to be specified directly. First, user needs to verify the existence of each of these 3 files. Next, user can specify these paths directly by either using

```shell
$ ccmake ..
```

or by directly running `CMake` with proper flags as following

```shell
$ cmake .. \
        -DGUROBI_CPP_LIB=/to/gurobiXXX/YY/lib/libgurobi_c++.a \
        -DGUROBI_INCLUDE_DIR=/to/gurobiXXX/YY/include \
        -DGUROBI_LIB=/to/gurobiXXX/YY/lib/libgurobiZZ.QQ
```
