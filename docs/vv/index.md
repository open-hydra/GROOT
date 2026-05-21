# Verification & Validation

This section documents the test suite used to verify GROOT's implementation of
the Discrete Transfer Method and to validate its results against reference
solutions.

## Test suite

<div class="grid cards" markdown>

-   :material-cube-outline:{ .lg .middle } __Homogeneous Cube__

    ---

    3-D grey medium in a unit cube. Reference by angular quadrature (same Nr as GROOT). Tests heatflux and volumetric source for κ = 0.1, 1, 10 m⁻¹.

    [:octicons-arrow-right-24: Homogeneous cube](homcube.md)

-   :material-mirror:{ .lg .middle } __Periodic-BC Cube__

    ---

    Half-cube [0, 0.5]³ with symmetry planes. Verifies that the symmetry (periodic) boundary condition reproduces the full-cube solution.

    [:octicons-arrow-right-24: Periodic-BC cube](cube-periodic.md)

-   :material-axis-arrow:{ .lg .middle } __Finite Cylinder__

    ---

    Axisymmetric cylinder (R = 1 m, L = 2 m, T = 3500 K). Reference by analytical integration (scipy). Tests lateral wall flux and on-axis source.

    [:octicons-arrow-right-24: Finite cylinder](finite-cylinder.md)

-   :material-fire:{ .lg .middle } __Axisymmetric Furnace__

    ---

    Axisymmetric furnace with non-uniform temperature field from a polynomial flame profile. Reference from a high-resolution legacy GROOT run.

    [:octicons-arrow-right-24: Axisymmetric furnace](furnace.md)

-   :material-shape-outline:{ .lg .middle } __Scalene Trapezoid__

    ---

    Non-orthogonal quadrilateral enclosure extruded as a thin slab. Exercises bilinear mesh support and non-axis-aligned faces.

    [:octicons-arrow-right-24: Scalene trapezoid](quadrilateral.md)

</div>

## Test matrix

| Test | Geometry | T field | κ values | Reference type | Tolerance |
|------|----------|---------|----------|----------------|-----------|
| [Homogeneous cube](homcube.md) | 3-D cube 1×1×1 m | uniform 1000 K | 0.1, 1, 10 m⁻¹ | Angular quadrature | 5 % |
| [Periodic-BC cube](cube-periodic.md) | Half-cube 0.5³ | uniform 1000 K | 0.1, 1, 10 m⁻¹ | Same as homcube | 10 % |
| [Finite cylinder](finite-cylinder.md) | Cyl. R=1, L=2 m | uniform 3500 K | 0.1, 1, 10 m⁻¹ | Analytical (scipy) | 5 % |
| [Furnace](furnace.md) | Cyl. R=0.45, L=5 m | non-uniform | fixed | Legacy GROOT run | 15 % |
| [Quadrilateral](quadrilateral.md) | Trapezoid 2.2×1.2 m | uniform 3000 K | 0.1, 1, 10 m⁻¹ | Legacy GROOT run | 20 % |

## Running the tests

After building GROOT (`./install.sh build --compiler=gnu`), each test
executable is placed in `test/<name>/bin/`.  Run from the test directory:

```bash
cd test/homcube
./bin/homcube_test
python verify.py --k 1.0 --plot
```

The `verify.py` script reads `OUTPUT/` (solver output) and `data/` (reference),
prints the maximum relative error for each quantity, and exits with code 1 if
any error exceeds the tolerance.
