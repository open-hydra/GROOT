# Running GROOT

This page describes the complete workflow for running a GROOT simulation: preparing the case directory, defining the mesh and boundary conditions, launching the solver, and inspecting the output.

---

## Simulation Workflow

A typical GROOT run follows three phases:

1. **Setup** — read the input file, load the gas-phase field and mesh, compute metrics and radiative properties.
2. **Solve** — run the Discrete Transfer Method (DTM) iteration loop.
3. **Postprocess** — write the radiative field and wall heat flux to output files.

All configuration is driven by a single INI file (`input.ini`).

---

## Case Directory Structure

Every GROOT case follows a standard layout:

```
my_case/
├── input.ini          ← solver configuration
├── INPUT/
│   ├── phase.txt      ← species list and molecular weights
│   └── bc.txt         ← boundary condition flags per face
└── OUTPUT/
    ├── field.*        ← (read) gas-phase field from MOSE (mesh + variables)
    ├── wall.tec       ← (read) wall temperatures from MOSE
    ├── rad-field.*    ← (written) radiative field solution (volume)
    └── rad-wall.*     ← (written) wall radiative flux (surface)
```

!!! note
    GROOT reads the computational mesh and gas-phase field from the MOSE solver output (`OUTPUT/field.*`). There is no separate mesh file: mesh node coordinates are embedded in the MOSE field file.

---

## Mesh

GROOT reads the computational mesh from `OUTPUT/field.*`, the same Tecplot-format multi-block file produced by MOSE. The format (ASCII or binary) is controlled by `sol-format` in `[MOSE-IO]`.

Each block must be a structured zone with node coordinates `(X, Y, Z)`. GROOT computes all mesh metrics (cell volumes, face areas and outward normals) internally from the node positions.

For **2-D axisymmetric** cases, provide a thin wedge mesh in the $x$–$y$ plane extruded by a small angle $\delta\theta$ in the $z$ direction (set `twoDax = .true.` and `ang = <degrees>` in `[GROOT-Discretization]`).

---

## Gas-Phase Field

The gas-phase field (temperature, density, species mass fractions) is read from `OUTPUT/field.*` — the MOSE solver output. The file format is set with `sol-format` in `[MOSE-IO]`.


---

## Species Definition

The file `INPUT/phase.txt` defines the species present in the gas mixture. It is read by GROOT to map CFD species to radiating species in the spectral database.

Format:

```
# Species    Mol. weight [g/mol]
N2           28.014
H2O          18.015
CO2          44.010
```


---

## Boundary Conditions

Wall boundary conditions (temperature and emissivity) are encoded in the `INPUT/bc.txt` file or set programmatically when using GROOT as a library. The BC type identifier stored per wall face determines how the wall interacts with rays.

| BC type | Description |
|---------|-------------|
| 200 | Azimuthal / periodic face (2-D axisymmetric) |
| 300 | Specular (mirror) — used for symmetry planes |
| 101/103 | Connection between blocks |
| Other | Wall |

---

## Running the Solver



```bash
export OMP_NUM_THREADS=4
ulimit -s unlimited
./bin/GROOT
```

!!! warning
    GROOT may require a large stack size for deep ray-tracing recursion. Always run `ulimit -s unlimited` (or `export KMP_STACKSIZE=100M` for Intel compilers) before launching the solver.

---

## Output

After a run completes, the `OUTPUT/` directory contains:

| File | Description |
|------|-------------|
| `rad-field.tec` | Volume solution |
| `rad-wall.tec` | Surface solution |

The output format is controlled by `sol-format` in `[GROOT-IO]`.

### Output Formats

| Value | Format |
|-------|--------|
| `tecplot ascii` | Tecplot ASCII (`.tec`) — default |
| `tecplot binary` | Tecplot binary (`.szplt`) — requires TecIO |
| `vtk ascii` | VTK ASCII (`.vtm`) |
| `vtk raw` | VTK binary (`.vtm`) |

### Output Variables

The radiative field file (`rad-field.*`) contains one zone per block. Variables depend on the spectral model:

=== "Gray / const"

    | Variable | Description |
    |----------|-------------|
    | `Ib` | Blackbody intensity $I_b = \sigma T^4 / \pi$ [W/m²/sr] |
    | `k` | Absorption coefficient [m⁻¹] |
    | `divq` | Radiative source term $\nabla \cdot \mathbf{q}_r$ [W/m³] |

=== "WSGG"

    | Variable | Description |
    |----------|-------------|
    | `Ib` | Blackbody intensity [W/m²/sr] |
    | `k(0)`…`k(4)` | Per-gray-gas absorption coefficients [m⁻¹] |
    | `a(0)`…`a(4)` | Per-gray-gas weights [–] |
    | `divq` | Radiative source term $\nabla \cdot \mathbf{q}_r$ [W/m³] |

The wall file (`rad-wall.*`) contains one zone per opaque boundary face:

| Variable | Description |
|----------|-------------|
| `emissivity` | Wall emissivity $\varepsilon$ [–] |
| `Q` | Net radiative wall heat flux [W/m²] |
| `a(0)`…`a(4)` | WSGG weights at the wall (WSGG models only) |

---

## Using GROOT as a Library

GROOT can be used as a static library (`libGROOTL.a`) from an external Fortran program. The main entry point is the `GROOT_simulation_type` derived type:

```fortran
use GROOT_Advanced_Types_m, only: GROOT_simulation_type
use GROOT_Wrap_Setup,        only: GROOT_setup
use GROOT_Wrap_Solve,        only: GROOT_solve
use GROOT_Wrap_Postprocess,  only: GROOT_postprocess

type(GROOT_simulation_type) :: simulation

call GROOT_setup(simulation)
call GROOT_solve(simulation)
call GROOT_postprocess(simulation)
```

Link against the library as described in the [Installation](../getting-started/installation.md#library-linking-advanced) page.

---
