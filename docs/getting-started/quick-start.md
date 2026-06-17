# Quick Start

This page walks you through running your first GROOT simulation — the **homogeneous cube** test, a classical benchmark for DTM radiation solvers in a 1 m³ box with uniform gas properties and black cold walls.

!!! note "Prerequisites"
    Make sure GROOT is built and the executable is available at `bin/GROOT`. See [Installation](installation.md) if you haven't built it yet.

---

## Test Case Overview

The homogeneous cube benchmark consists of a 1 m³ cubic enclosure filled with a uniform absorbing gas at 1000 K. All six walls are black and cold ($T_w = 0$ K). An analytical solution for the wall heat flux and the radiative source term is available, making this an ideal first verification test.

| Parameter | Value |
|-----------|-------|
| Domain | 1 × 1 × 1 m³ |
| Mesh | 30 × 30 × 30 cells |
| Gas temperature | 1000 K |
| Absorption coefficient $\kappa$ | 1 m⁻¹ |
| Wall condition | Black, cold ($T_w = 0$ K) |
| Rays | 256 |

---

## Directory Layout

The homcube test is pre-configured in the repository:

```
test/homcube/
├── homcube.f90        ← self-contained test program (no input.ini needed)
├── CMakeLists.txt     ← build configuration
├── data/              ← reference solutions for comparison
│   ├── ref_source*
│   └── ref_heatflux*
└── OUTPUT/            ← solver output (created on run)
    ├── source_centerline_k*.dat
    └── heatflux_face1_k*.dat
```

The test builds a 30 × 30 × 30 mesh and uniform field entirely in memory — no external mesh or field files are needed.

---

## Running the Test

After building GROOT, the `homcube` executable is placed in `test/homcube/bin/`. Run it from the test directory:

```bash
cd test/homcube
./bin/homcube
```

!!! tip "Parallel execution"
    If GROOT was built with OpenMP:
    ```bash
    export OMP_NUM_THREADS=4
    ulimit -s unlimited
    ./bin/homcube
    ```

During execution, GROOT prints the iteration history and CPU time to the terminal. A successful run ends with:

```
 Done. See OUTPUT/
```

---

## Checking the Results

After the run, the `OUTPUT/` directory contains:

| File | Content |
|------|---------|
| `source_centerline_k<k>.dat` | Radiative source term $\nabla \cdot \mathbf{q}_r$ along the z-axis centreline (one file per absorption coefficient) |
| `heatflux_face1_k<k>.dat` | Wall heat flux $Q$ on face 1 (x-min) vs y-coordinate (one file per absorption coefficient) |

Compare against the reference data in `data/`.
A successful run shows the GROOT solution matching the reference to within a few percent (mesh-dependent).

---

## Running the Full Test Suite

All validation cases can be built and run from the repository root after building with CMake:

```
test/
├── homcube/           ← homogeneous cube (3-D, uniform T)
├── cube_periodic/     ← cube with periodic boundary conditions
├── finite_cylinder/   ← axisymmetric cylinder (2-D wedge)
├── furnace_source/    ← axisymmetric furnace, non-uniform T
├── quadrilateral/     ← scalene trapezoid enclosure
└── furnace3d/         ← 3-D box furnace (Selcuk 1985)
```

Each test is a standalone executable with its own `bin/` and `OUTPUT/` directory. Run any test by executing its binary from within the test directory. Reference data are in `data/` for post-processing comparison.

---

## Next Steps

Now that GROOT is running, explore:

- **[User Guide](../user/index.md)** — running GROOT on your own meshes with `input.ini`
- **[Input Reference](../user/input/input-parameters.md)** — full parameter documentation
- **[Theory Guide](../theory/index.md)** — DTM algorithm and spectral models

---
