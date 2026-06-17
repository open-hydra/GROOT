# Code Structure

This page documents the repository layout, the internal architecture
of the GROOT library, and the solver execution pipeline.  All diagrams
use [Mermaid](https://mermaid.js.org/) and render directly in the
documentation.

---

## Repository Layout

```
GROOT/
├── CMakeLists.txt          # Top-level CMake build
├── CMakePresets.json        # Developer presets (compilers, paths)
├── install.sh               # Build / compile / update helper
├── mkdocs.yml               # Documentation site configuration
│
├── src/
│   ├── app/                 # Executables
│   │   └── main.f90         # GROOT solver entry point
│   └── lib/                 # GROOT library (GROOTL)
│       ├── base/            # Fundamental types and global parameters
│       ├── config/          # Input parsing, registry, setup assignment
│       ├── driver/          # High-level solver orchestration
│       ├── io/              # File I/O (solution, walls)
│       └── physics/         # Radiation physics
│           ├── Lib_radproperties.f90  #   Spectral models (gray, WSGG)
│           ├── lib_geometry.f90       #   Face normals and cell volumes
│           ├── lib_math.f90           #   Utility math routines
│           ├── lib_raytracing.f90     #   DTM ray tracing
│           └── lib_rte.f90            #   RTE solve loop
│
├── lib/                     # External dependencies (git submodules)
│   ├── ORION/               # Structured-grid I/O library
│   └── third_party/
│       └── FiNeR/           # INI file parser
│
├── test/                    # Validation test suite
│   ├── homcube/             #   Homogeneous cube (3-D)
│   ├── cube_periodic/       #   Cube with periodic BCs
│   ├── finite_cylinder/     #   Axisymmetric cylinder
│   ├── furnace_source/      #   Axisymmetric furnace, non-uniform T
│   ├── quadrilateral/       #   Scalene trapezoid enclosure
│   └── furnace3d/           #   3-D box furnace (Selcuk 1985)
│
├── docs/                    # MkDocs documentation source
├── cmake/                   # CMake modules (compiler flags, OpenMP, etc.)
├── bin/                     # Built executables (GROOT)
└── build/                   # Build artefacts
```

---

## Dependency Graph

```mermaid
graph TD
    GROOT["<b>GROOT</b><br/>DTM radiation solver"]
    ORION["<b>ORION</b><br/>Structured-grid I/O"]
    FiNeR["<b>FiNeR</b><br/>INI parser"]

    GROOT --> ORION
    GROOT --> FiNeR

    style GROOT fill:#37474f,stroke:#cfd8dc,color:#fff
    style ORION fill:#1565c0,stroke:#90caf9,color:#fff
    style FiNeR fill:#6a1b9a,stroke:#ce93d8,color:#fff
```

Optional compile-time dependencies:

```mermaid
graph LR
    GROOT["GROOT"]
    OMP["OpenMP"]
    TecIO["TecIO"]

    GROOT -.->|USE_OPENMP| OMP
    GROOT -.->|USE_TECIO| TecIO

    style GROOT fill:#37474f,stroke:#cfd8dc,color:#fff
    style OMP fill:#00695c,stroke:#80cbc4,color:#fff
    style TecIO fill:#00695c,stroke:#80cbc4,color:#fff
```

---

## Library Architecture

The GROOT library (`libGROOTL`) is organised in five layers.  Lower
layers have no knowledge of higher layers.

```mermaid
graph TB
    subgraph driver ["<b>driver/</b> — Solver orchestration"]
        Procedures["Procedures_m<br/><i>GROOT_simulation_type</i>"]
        WSetup["Wrap_Setup"]
        WSolve["Wrap_Solve"]
        WPost["Wrap_Postprocess"]
        Alloc["Mod_Allocate_Data"]
        Metrics["Mod_Metrics"]
        Bound["Mod_Bound<br/><i>wall BCs</i>"]
        Phase["Mod_Phase<br/><i>species / WSGG bands</i>"]
    end

    subgraph config ["<b>config/</b> — Input & configuration"]
        ReadIni["Read_Ini"]
        AssignSetup["Assign_Setup"]
        ConfigTypes["Config_Types_m"]
    end

    subgraph physics ["<b>physics/</b> — Radiation physics"]
        RadProp["Lib_radproperties<br/><i>spectral models</i>"]
        Geom["Lib_geometry<br/><i>face normals · volumes</i>"]
        Math["Lib_math<br/><i>utility routines</i>"]
        Ray["Lib_raytracing<br/><i>DTM ray integration</i>"]
        RTE["Lib_rte<br/><i>RTE solve loop</i>"]
    end

    subgraph io ["<b>io/</b> — File I/O"]
        IOSol["IO_Solution"]
        IOWall["IO_Wall"]
    end

    subgraph base ["<b>base/</b> — Fundamental types"]
        AdvTypes["Advanced_Types_m<br/><i>GROOT_simulation_type</i>"]
        ModType["Mod_Type<br/><i>obj_grid, obj_block</i>"]
        Params["Mod_Parameters<br/><i>global scalars</i>"]
    end

    WSetup --> config & io & base
    WSolve --> RTE
    RTE --> Ray & RadProp
    Ray --> Geom & Math & Metrics
    WSetup --> Bound & Phase
    config --> base
    io --> base

    style driver fill:#263238,stroke:#90a4ae,color:#eceff1
    style config fill:#1a237e,stroke:#7986cb,color:#e8eaf6
    style physics fill:#4a148c,stroke:#ce93d8,color:#f3e5f5
    style io fill:#bf360c,stroke:#ff8a65,color:#fbe9e7
    style base fill:#3e2723,stroke:#a1887f,color:#efebe9
```

---

## Module Hierarchy

Each source file contains a Fortran module following a consistent
naming convention:

| Prefix | Role | Example |
|--------|------|---------|
| `Mod_*` | Module defining types, data, and procedure pointers | `Mod_Allocate_Data`, `Mod_Bound` |
| `Lib_*` | Library of pure computational routines | `Lib_rte`, `Lib_raytracing` |
| `*_m` | Fundamental type / parameter modules | `Advanced_Types_m`, `Config_Types_m` |
| `Wrap_*` | High-level driver wrappers | `Wrap_Setup`, `Wrap_Solve` |
| `IO_*` | File I/O routines | `IO_Solution` |
| `Read_*` | Input file parsing | `Read_Ini` |

All public symbols are prefixed with `GROOT_` to avoid namespace
collisions when GROOT is linked as a library.

---

## Solver Pipeline

The main program (`src/app/main.f90`) creates a `GROOT_simulation_type` object
and calls three phases: **setup**, **solve**, and **postprocess**.

```mermaid
sequenceDiagram
    participant Main as main.f90
    participant Setup as Wrap_Setup
    participant Solve as Wrap_Solve
    participant Post as Wrap_Postprocess

    Main->>Setup: GROOT_setup(simulation)
    activate Setup
    Setup->>Setup: Read_Inifile
    Setup->>Setup: Assign_Setup (spectral model, BC types)
    Setup->>Setup: Load_Gas_Field (ORION read)
    Setup->>Setup: Setup_Data_Structure
    Setup->>Setup: Setup_Metrics (volumes, face normals)
    Setup->>Setup: Setup_BC (wall temperatures, emissivities)
    Setup->>Setup: Setup_Rad_Properties (κ, Ib, WSGG weights)
    deactivate Setup

    Main->>Solve: GROOT_solve(simulation)
    activate Solve
    Solve->>Solve: RTE solve loop (DTM iterations)
    Solve->>Solve: Compute wall fluxes and source terms
    deactivate Solve

    Main->>Post: GROOT_postprocess(simulation)
    Note right of Post: Write rad-field.tec + rad-wall.tec
```

### DTM Solve Loop

Each call to the RTE solver performs the Discrete Transfer Method
iteration until the wall-flux residual falls below `tol`:

```mermaid
flowchart TB
    Start([Wrap_Solve entry]) --> Init[Initialise flux accumulators]
    Init --> RayLoop

    subgraph RayLoop ["DTM ray loop (Nr directions)"]
        direction TB
        SampleDir["Sample ray direction<br/>(uniform sphere)"] --> TraceRay
        TraceRay["Trace ray through cells<br/>(Lib_raytracing)"] --> AccumFlux
        AccumFlux["Accumulate wall flux Q<br/>and source term ∇·q"]
    end

    RayLoop --> Converge{Residual<br/>< tol?}
    Converge -->|No, iter < itermax| RayLoop
    Converge -->|Yes| End([Return])

    style Start fill:#37474f,stroke:#cfd8dc,color:#fff
    style End fill:#37474f,stroke:#cfd8dc,color:#fff
    style RayLoop fill:#4a148c,stroke:#ce93d8,color:#f3e5f5
```

---

## Data Structures

### Simulation container

The top-level type is `GROOT_simulation_type`, which wraps the
computational domain and configuration:

```mermaid
classDiagram
    class GROOT_simulation_type {
        +config : config_type
        +domain : obj_grid
        +OCP : ORION_data
        +ORP : ORION_data
    }
    class config_type {
        +nthreads : integer
        +Nr : integer
        +source : logical
        +cputime(2) : real64
    }
    class obj_grid {
        +Nb : integer
        +block(:) : obj_block
    }
    class obj_block {
        +ka(:,:,:,:) : real64
        +ks(:,:,:,:) : real64
        +Ib(:,:,:) : real64
        +source(:,:,:) : real64
        +G(:,:,:) : real64
        +face(6) : obj_face
        +metrics : d_metrics_type
    }

    GROOT_simulation_type --> config_type : config
    GROOT_simulation_type --> obj_grid : domain
    obj_grid "1" --> "*" obj_block : block
```

| Array | Shape | Content |
|-------|-------|---------|
| `ka` | `(Nx, Ny, Nz, 0:Ngg)` | Absorption coefficient per spectral band |
| `ks` | `(Nx, Ny, Nz, 0:Ngg)` | Scattering coefficient per spectral band |
| `Ib` | `(Nx, Ny, Nz)` | Blackbody intensity $\sigma T^4 / \pi$ |
| `source` | `(Nx, Ny, Nz)` | Radiative source term $\nabla \cdot \mathbf{q}_r$ |
| `G` | `(Nx, Ny, Nz)` | Irradiation (incident radiation) |

### Face data

Each of the six faces of a block stores wall quantities:

| Array | Shape | Content |
|-------|-------|---------|
| `face(f)%Q` | `(n1, n2, 0:Ngg)` | Spectral wall heat flux [W/m²] per band |
| `face(f)%Qtot` | `(n1, n2)` | Total wall heat flux [W/m²] |
| `face(f)%T` | `(n1, n2)` | Wall temperature [K] |
| `face(f)%eps` | `(n1, n2)` | Wall emissivity |
| `face(f)%bc` | `(n1, n2)` | BC type identifier |

---

## Build System

### CMake targets

```mermaid
graph LR
    GROOTL["<b>GROOTL</b><br/>(Fortran library)"]
    GROOT_EXE["<b>GROOT</b><br/>(executable)"]
    ORION_LIB["ORION"]
    FiNeR_LIB["FiNeR"]

    GROOT_EXE --> GROOTL
    GROOTL --> ORION_LIB
    GROOTL --> FiNeR_LIB

    style GROOTL fill:#37474f,stroke:#cfd8dc,color:#fff
    style GROOT_EXE fill:#1565c0,stroke:#90caf9,color:#fff
    style ORION_LIB fill:#6a1b9a,stroke:#ce93d8,color:#fff
    style FiNeR_LIB fill:#6a1b9a,stroke:#ce93d8,color:#fff
```

| Target | Type | Description |
|--------|------|-------------|
| `GROOTL` | Static library | Core solver + all physics |
| `GROOT` | Executable | Standalone solver (`src/app/main.f90`) |

### Build workflow

```bash
# Option A: install.sh (recommended for first build)
./install.sh build --compilers=gnu

# Option B: CMake presets (for iterative development)
./install.sh compile
# or equivalently:
cmake --preset default && cmake --build build
```

Key CMake options:

| Option | Default | Effect |
|--------|:-------:|--------|
| `USE_OPENMP` | OFF | Enable OpenMP threading |
| `USE_TECIO` | OFF | Enable TecIO binary output |

### Dependency paths

| Variable | Default | Library |
|----------|---------|---------|
| `ORION_PATH` | `lib/ORION/` | ORION I/O |
| `FINER_PATH` | `lib/third_party/FiNeR/` | FiNeR INI parser |

---

## Naming Conventions

| Convention | Example | Meaning |
|------------|---------|---------|
| `GROOT_` prefix | `GROOT_Config_Types_m` | Public Fortran module |
| `obj_` prefix | `obj_dtm`, `obj_rad_model` | Global configuration singleton |
| `Lib_` prefix | `Lib_rte`, `Lib_raytracing` | Computational routine library |
| `Mod_` prefix | `Mod_Allocate_Data` | Module with types or driver logic |
| `Wrap_` prefix | `Wrap_Solve` | Driver-level wrapper |
| `_m` suffix | `Config_Types_m` | Fundamental type module |
| `_type` suffix | `GROOT_simulation_type` | Derived type |
