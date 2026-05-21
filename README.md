<p align="center">
  <h1 align="center">GROOT</h1>
  <p align="center"><b>Grid based RadiatiOn fOrtran Toolbox</b></p>
</p>

<p align="center">
  <a href="https://github.com/open-hydra/GROOT/blob/main/LICENSE"><img src="https://img.shields.io/badge/license-GPLv3-blue.svg" alt="License: GPLv3"></a>
  <a href="https://open-hydra.github.io/GROOT/"><img src="https://img.shields.io/badge/docs-online-brightgreen.svg" alt="Documentation"></a>
  <img src="https://img.shields.io/badge/language-Fortran-734f96.svg" alt="Language: Fortran">
</p>

---

GROOT is an open-source radiative transfer equation solver on multi-block structured grids. Written in modern Fortran, it implements the discrete transfer method for gray and specral (WSGG) gases.

## Features

- **Discrete Transfer Method** - Heat flux and radiative source term calculation with full ray tracing
- **Weighted sum of gray gases approach** - WSGG models calibrated for high-pressure and temperature combustion products
- **Parallel execution**— shared-memory parallelism via OpenMP.
- **Flexible I/O** — solution output in Tecplot (ASCII and binary) and VTK formats.

## Quick Start

| Requirement | Details |
|---|---|
| **CMake** | ≥ 3.23 |
| **Fortran compiler** | GNU (`gfortran`) or Intel/oneAPI (`ifort` / `ifx`) |



### Build

```bash
git clone --recurse-submodules https://github.com/open-hydra/GROOT.git
cd GROOT

# Build with GNU compilers and OpenMP
./install.sh build --compiler=gnu --use-openmp

```

The executable is placed in `bin/GROOT`.

See the [Installation Guide](https://open-hydra.github.io/GROOT/getting-started/installation/) for all build options, CMake presets, and troubleshooting.

### Run the Finite Cylinder Test Case

```bash
cd test/finite_cylinder
./bin/cylinder
```

See the [Quick Start](https://open-hydra.github.io/GROOT/getting-started/quick-start/) for a full walkthrough.

## Dependencies

GROOT is built on top of several companion libraries, included as Git submodules:

| Library | Role |
|---|---|
| [ORION](https://github.com/MarcoGrossi92/ORION) | Multi-format I/O (Tecplot, VTK, Plot3D) |
| [FiNeR](https://github.com/szaghi/FiNeR) | INI configuration file parser |

Optional external libraries: **OpenMP**.

## Project Structure

```
GROOT/
├── src/
│   ├── app/           # Main application
│   └── lib/           # Solver library sources
├── lib/               # Git submodule dependencies
├── test/              # Verification & validation cases
│   └── ...
├── docs/              # MkDocs documentation source
├── cmake/             # CMake modules
├── install.sh         # Build helper script
└── CMakeLists.txt
```

## Documentation

Full documentation is available at **[open-hydra.github.io/GROOT](https://open-hydra.github.io/GROOT/)**, covering:

- Installation & quick start
- User guide & input file reference
- Verification & validation results
- Theory guide (governing equations, numerical methods, emissivity models)

## License

GROOT is free and open-source software released under the [GNU General Public License v3.0](LICENSE).
