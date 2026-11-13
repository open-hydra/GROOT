# Installation

This document describes how to obtain and build **GROOT**. The instructions describe the current `install.sh` script, CMake configuration and library layout.

!!! note
    GROOT has a dual nature: it is both a library and an executable. The installation process described here produces both the static library `libGROOT.a` and the main executable `bin/GROOT`. If you are only interested in using GROOT as a library, you can link against `libGROOT.a` without caring about the executable.

## Prerequisites

Before attempting to build GROOT make sure your system provides the following external tools and compilers:

- **CMake** – 3.23 or newer.
- **Fortran compiler** – either the GNU toolchain (`gfortran`) or Intel/oneAPI (`ifort`/`ifx`) are supported..
- **OpenMP** – needed for optional parallelization support.

### Libraries

GROOT builds against two libraries. The releases it is tested with are in
[Versioning](../development/versioning.md).

**ORION**, by release. [ORION](https://github.com/MarcoGrossi92/ORION) 2.1.0 or a later 2.x
(`GROOT_ORION_VERSION` in `CMakeLists.txt`). The first configuration selects the newest such release on GitHub and
checks it out in `lib/ORION` (git-ignored); the file `deps.lock` (git-ignored) keeps it, so that the next
configurations use it without the network. `-DSEMVER_UPDATE=ALL` (`./install.sh build --update-deps`) selects the
newest release again. `-DFETCHCONTENT_SOURCE_DIR_ORION=<path>` (`--include-orion=<path>`) builds an external clone
instead, which must be a release GROOT accepts.

**FiNeR**, the INI file parser, as the Git submodule `lib/FiNeR` (v2.3.0).
`-DFETCHCONTENT_SOURCE_DIR_FINER=<path>` (`--include-finer=<path>`) builds an external clone of that release instead.
FiNeR's own libraries (PENF, FACE, FLAP, BeFoR64, StringiFor) are cloned by the configure, in `lib/FiNeR/src/third_party`,
from their default branches. If an upstream change breaks the build, configure with `-DGROOT_FINER_DEPS=BACKUP` to use
the releases GROOT is tested with.

If the configure stops because `lib/FiNeR` is at another commit (after a `git checkout` of another branch), run
`git submodule update --init` and configure again.

**radlib** *(optional)*, `lib/third_party/radlib` (`https://github.com/BYUignite/radlib`): the SLW spectral model
(RC-SLW). It is a Git submodule marked `update = none`, so it is **not** pulled by `git submodule update`: it is fetched
on demand only when you build with `--use-radlib` (see below).


## Build methods

First clone the repository with its submodules:

```bash
git clone --recurse-submodules https://github.com/open-hydra/GROOT.git
cd GROOT
```

To fully install GROOT, you may either use the bundled install script or invoke CMake manually. The script is convenient and is the preferred route for most users.

### Build with `install.sh` (recommended)

The script exposes two commands: `build` and `compile`.  It also maintains a `CMakePresets.json` file that records the configuration used for
the most recent `build` invocation.

```bash
./install.sh [GLOBAL_OPTIONS] COMMAND [COMMAND_OPTIONS]
```

**Global options**

* `-v`, `--verbose` – enable verbose logging.

**`build` command**

Performs a clean configure+build cycle.  Example usage:

```bash
# minimal GNU build with OpenMP enabled
./install.sh build --compilers=gnu --use-openmp

# full configuration with Intel compilers and all optional features
./install.sh build --compilers=intel \
                  --use-openmp\
```

Options accepted by `build`

* `--compilers=<gnu|intel>` – select the compiler family (default: `gnu`).
* `--use-openmp` – enable OpenMP parallelization.
* `--use-radlib` – enable the optional **SLW** spectral model via radlib. Fetches
  the `lib/third_party/radlib` submodule and configures with `-DUSE_RADLIB=ON`.
* `--include-orion=PATH` – build an external ORION clone instead of the release GROOT selects.
* `--update-deps` – select again the newest ORION release that GROOT accepts (it replaces the one of `deps.lock`).
* `--include-finer=PATH` – build an external FiNeR clone instead of the submodule `lib/FiNeR`.

The script sets up environment variables for the chosen compilers and then invokes CMake. After a successful build, a `CMakePresets.json` file is written in the source root so that subsequent compilations can reuse the configuration.

<!-- The release path, build artifacts etc.--> 

**`compile` command**

Re‑runs CMake using the previously generated preset and rebuilds the project without clearing the build directory. This is useful during development when only the source has changed. Example usage:

```bash
./install.sh compile
```

### Build with CMake

If you prefer fine‑grained control, perform the configuration yourself. This is essentially what `install.sh` does under the hood.

```bash
mkdir build && cd build
cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_Fortran_COMPILER=gfortran \  # or ifx
    -DUSE_OPENMP=ON \                    # or OFF
    -DFETCHCONTENT_SOURCE_DIR_ORION=/path/to/ORION \   # optional: an external clone instead of the release
    -DFETCHCONTENT_SOURCE_DIR_FINER=/path/to/FiNeR     # optional
cmake --build . --parallel
```

The resulting artifacts are placed in `build/` by default. The static library is `lib/libGROOT.a` and the main executable is `bin/GROOT` (inside the build
directory unless you set `CMAKE_INSTALL_PREFIX`).

## CMake presets

The file `CMakePresets.json` produced by the install script records all of the cache variables that were used during configuration. You can build the project later simply with

```bash
cmake --preset default
cmake --build build
```

or using the `compile` command of the install script as described above.

## Optional components

### SLW spectral model (radlib)

The **SLW** (RC-SLW) spectral model is provided by the external
[radlib](https://github.com/BYUignite/radlib) library (MIT-licensed) and is
**off by default**. Enable it with:

```bash
./install.sh build --compilers=gnu --use-openmp --use-radlib
```

This fetches the `lib/third_party/radlib` submodule (if not already present)
and configures with `-DUSE_RADLIB=ON`. Equivalently, with a manual CMake build:

```bash
# --checkout is required: radlib is marked `update = none` in .gitmodules so it is
# skipped by the default (recursive) submodule update; --checkout overrides that.
git submodule update --init --checkout lib/third_party/radlib
cmake -B build -DUSE_RADLIB=ON ...
```

Notes:

* The ALBDF data files (`*.bin`) are **generated automatically** at build time
  from the shipped `*.txt` tables (via radlib's `writebin.x`); no extra step or
  `make install` of radlib is required.
* radlib is built as a **shared** library (`libradlib.so`); the GROOT executable
  locates it via its build-tree rpath. When relocating/installing the executable,
  make sure the shared library is reachable.
* Once built, select the model with `model = slw` in `[GROOT-Model]`
  (see [Input Parameters](../user/input/input-parameters.md#groot-model)).
* Building **without** `--use-radlib` is fully supported; only `model = slw`
  becomes unavailable.

## Library linking (advanced)

To use GROOT from an external Fortran program you can compile as follows:

```bash
gfortran -I/path/to/GROOT/include \
         -L/path/to/GROOT/lib \
         -lGROOT \
         your_program.f90 -o your_program
```

or, from a CMake project:

```cmake
find_package(GROOT REQUIRED)
add_executable(myapp main.f90)
target_link_libraries(myapp GROOT::GROOT)
```

Installation prefix and other details may be customised via standard CMake
variables such as `CMAKE_INSTALL_PREFIX` and by running `cmake --install`.

## Next steps

* **[Quick Start Tutorial](quick-start.md)** – build and verify the installation.

---