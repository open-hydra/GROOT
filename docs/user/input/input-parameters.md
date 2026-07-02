# GROOT Input Parameters

## [MOSE-IO]

Controls the format of the gas-phase field file written by MOSE and read by GROOT (`OUTPUT/field.*`).

| Parameter | Default | Allowed | Required | Description |
|-----------|---------|---------|----------|-------------|
| sol-format | `tecplot ascii` | `tecplot ascii`, `tecplot binary`, `vtk ascii`, `vtk raw` | no | Format and encoding of the MOSE gas-phase field file |

!!! warning "VTK input not yet fully supported"
    Reading the gas-phase field in VTK format (`vtk ascii` / `vtk raw`) is not yet functional.
    Use `tecplot ascii` or `tecplot binary` for `[MOSE-IO] sol-format`.
    The wall file (`OUTPUT/wall.tec`) is always read as Tecplot ASCII regardless of this setting.

## [GROOT-IO]

Controls the format of the radiative solution files written by GROOT (`OUTPUT/rad-field.*`, `OUTPUT/rad-wall.*`).

| Parameter | Default | Allowed | Required | Description |
|-----------|---------|---------|----------|-------------|
| sol-format | `tecplot ascii` | `tecplot ascii`, `tecplot binary`, `vtk ascii`, `vtk raw` | no | Format and encoding of the GROOT output files |

!!! warning "VTK output experimental"
    VTK output (`vtk ascii` / `vtk raw`) for `rad-field` and `rad-wall` is implemented but not yet validated.
    Prefer `tecplot ascii` or `tecplot binary` for production runs.

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
| model | `const` | `const`, `gray`, `wsgg-H2O`, `wsgg-H2OCO2`, `snbw`, `snb`, `slw` | no | Spectral model for the gas-phase absorption coefficient |
| wall_emissivity | 1.0 | >= 0 | yes | Wall emissivity (uniform, applied to all opaque walls) |
| k | 1.0 | > 0 | no | Absorption coefficient [m⁻¹] — used only when `model = const` |
| slw_ngray | 4 | 1–25 | no | Number of gray gases for the SLW model — used only when `model = slw` |
| p_ref | 1.0e5 | > 0 | no | Reference pressure [Pa] for SLW wall emission — used only when `model = slw` |

!!! note "`model = slw` requires radlib"
    The SLW model is available only when GROOT is built with the optional
    **radlib** dependency (`./install.sh build --use-radlib`, or
    `-DUSE_RADLIB=ON`). A build without radlib will abort at runtime if
    `model = slw` is selected. See
    [Spectral Models → SLW](../../theory/spectral-models.md#spectral-line-weighted-sum-of-gray-gases-slw)
    and [Installation → Optional components](../../getting-started/installation.md#optional-components).
