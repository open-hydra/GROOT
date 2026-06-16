# GROOT Input Parameters

## [MOSE-IO]

Controls the format of the gas-phase field file written by MOSE and read by GROOT (`OUTPUT/field.*`).

| Parameter | Default | Allowed | Required | Description |
|-----------|---------|---------|----------|-------------|
| sol-format | `tecplot ascii` | `tecplot ascii`, `tecplot binary`, `vtk ascii`, `vtk raw` | no | Format and encoding of the MOSE gas-phase field file |

## [GROOT-IO]

Controls the format of the radiative solution files written by GROOT (`OUTPUT/rad-field.*`, `OUTPUT/rad-wall.*`).

| Parameter | Default | Allowed | Required | Description |
|-----------|---------|---------|----------|-------------|
| sol-format | `tecplot ascii` | `tecplot ascii`, `tecplot binary`, `vtk ascii`, `vtk raw` | no | Format and encoding of the GROOT output files |

## [GROOT-Discretization]

| Parameter | Default | Allowed | Required | Description |
|-----------|---------|---------|----------|-------------|
| Nr | 256 | > 0 | yes | Number of ray directions (uniform sphere sampling) |
| tol | 1.0e-3 | > 0 | yes | Convergence tolerance on the wall-flux residual |
| itermax | 10 | >= 1 | no | Maximum number of DTM iterations |
| ang | 1.0 | > 0 | no | Wedge half-angle [deg] for 2-D axisymmetric cases (`twoDax = .true.`) |

## [GROOT-Options]

| Parameter | Default | Allowed | Required | Description |
|-----------|---------|---------|----------|-------------|
| source | `.false.` | `.true.`, `.false.` | no | Compute volumetric radiative source term $\nabla \cdot \mathbf{q}_r$ |
| twoDax | `.false.` | `.true.`, `.false.` | no | 2-D axisymmetric mode: domain is a thin wedge in the $x$–$r$ plane |
| optically_thin | `.false.` | `.true.`, `.false.` | no | Optically-thin approximation for the **volumetric source term only**: sets $G=0$, i.e. $\nabla\cdot\mathbf{q}_r \approx 4\pi\kappa I_b$. Wall heat flux is unaffected. |

## [GROOT-Model]

| Parameter | Default | Allowed | Required | Description |
|-----------|---------|---------|----------|-------------|
| model | `const` | `const`, `gray`, `wsgg-H2O`, `wsgg-H2OCO2`, `snbw`, `snb` | no | Spectral model for the gas-phase absorption coefficient |
| wall_emissivity | 1.0 | >= 0 | yes | Wall emissivity (uniform, applied to all opaque walls) |
| k | 1.0 | > 0 | no | Absorption coefficient [m⁻¹] — used only when `model = const` |
