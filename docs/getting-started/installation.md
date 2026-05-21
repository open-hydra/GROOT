# Installation

This document describes how to obtain and build **GROOT**. The instructions describe the current `install.sh` script, CMake configuration and Git submodule layout.

!!! note
    GROOT has a dual nature: it is both a library and an executable. The installation process described here produces both the static library `libGROOT.a` and the main executable `bin/GROOT`. If you are only interested in using GROOT as a library, you can link against `libGROOT.a` without caring about the executable.

## Prerequisites

Before attempting to build GROOT make sure your system provides the following external tools and compilers:

- **CMake** – 3.23 or newer.
- **Fortran compiler** – either the GNU toolchain (`gfortran`) or Intel/oneAPI (`ifort`/`ifx`) are supported..
- **OpenMP** – needed for optional parallelization support.

### Git submodules

GROOT depends on several repositories that are included as Git submodules.

| Path                               | Repository URL                                        | Purpose                                  |
|------------------------------------|-------------------------------------------------------|------------------------------------------|
| `lib/ORION`                        | `https://github.com/MarcoGrossi92/ORION.git`         | I/O routines (TecIO, VTK, Plot3D, etc.)  |
| `lib/third_party/FiNeR`            | `https://github.com/szaghi/FiNeR.git`                | INI file parser                      |


## Build methods

First clone the repository with submodules:

```bash
git clone https://github.com/openhydra/GROOT.git
cd GROOT
# initialise submodules
git submodule update --init --recursive
```

To fully install GROOT, you may either use the bundled install script or invoke CMake manually. The script is convenient and is the preferred route for most users.

### Build with `install.sh` (recommended)

The script exposes three commands: `build`, `compile`, and `update`.  It also maintains a `CMakePresets.json` file that records the configuration used for
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
./install.sh build --compiler=gnu --use-openmp

# full configuration with Intel compilers and all optional features
./install.sh build --compiler=intel \
                  --use-openmp\
```

Options accepted by `build`

* `--compiler=<gnu|intel>` – select the compiler family (default: `gnu`).
* `--use-openmp` – enable OpenMP parallelization.
* `--include-orion=PATH` – use an external ORION tree instead of the submodule.
* `--include-finer=PATH` – same for FiNeR.

The script sets up environment variables for the chosen compilers and then invokes CMake. After a successful build, a `CMakePresets.json` file is written in the source root so that subsequent compilations can reuse the configuration.

<!-- The release path, build artifacts etc.--> 

**`compile` command**

Re‑runs CMake using the previously generated preset and rebuilds the project without clearing the build directory. This is useful during development when only the source has changed. Example usage:

```bash
./install.sh compile
```

**`update` command**

Synchronises the Git submodules. By default it checks out the commit recorded in `.gitmodules`; passing `--remote` will fetch the latest commit from each remote branch.

```bash
./install.sh update            # sync to recorded commit
./install.sh update --remote   # update to newest remote commit
```

### Build with CMake

If you prefer fine‑grained control, perform the configuration yourself. This is essentially what `install.sh` does under the hood.

```bash
mkdir build && cd build
cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_Fortran_COMPILER=gfortran \  # or ifx
    -DUSE_OPENMP=ON \                    # or OFF
    -DORION_PATH=/path/to/ORION \        # optional
    -DFINER_PATH=/path/to/FiNeR          # optional
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