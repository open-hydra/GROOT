# Spectral Models

This page describes the spectral models available in GROOT for computing
the **absorption coefficient** $\kappa$ and the **blackbody weighting**
needed to integrate the RTE over the full spectrum.

The spectral model is selected via the `model` parameter in
`[GROOT-Model]`.

---

## Gray Gas (constant $\kappa$)

**Model name:** `const`

The simplest model assumes a **spectrally constant** absorption coefficient:

$$
\kappa(\mathbf{r}, \eta) = \kappa_0 = \text{const}
$$

The RTE is solved once for the total (gray) intensity.  The absorption
coefficient is supplied directly by the user via the `k` parameter:

```ini
[GROOT-Model]
model = const
k     = 1.0   # [m⁻¹]
```

**Applicability:** benchmarking, scoping studies, and problems where a
known effective absorption coefficient is available.

---

## Planck-Mean Gray Gas

**Model name:** `gray`

The `gray` model computes a **temperature- and pressure-dependent**
absorption coefficient for each radiating species using tabulated
Planck-mean values:

$$
\kappa(\mathbf{r}) = \sum_s x_s\, \frac{p}{p_0}\, \kappa_{\mathrm{Planck},s}(T)
$$

where $x_s$ is the molar fraction of species $s$, $p_0 = 10^5$ Pa, and
$\kappa_{\mathrm{Planck},s}(T)$ is read from an internal polynomial
database.  Supported species: H₂O, CO₂, CO, CH₄.

A single (gray) RTE is solved with this $\kappa$.  The model is fast but
can overestimate the total heat flux because it weights the absorption
coefficient by the Planck function rather than by the actual spectral
distribution encountered along each ray.

```ini
[GROOT-Model]
model           = gray
wall_emissivity = 0.85
```

---

## Weighted-Sum-of-Gray-Gases (WSGG)

**Model names:** `wsgg-H2O`, `wsgg-H2OCO2`

The **WSGG** model approximates the non-gray absorption of real
participating gases by a weighted sum of $N_g$ gray gases plus one
transparent "window" band ($j = 0$):

$$
\varepsilon_g = \sum_{j=0}^{N_g} a_j(T)\, \left(1 - e^{-\kappa_j\, L}\right)
$$

where $a_j(T)$ are temperature-dependent weighting factors and $\kappa_j$
are the effective gray-gas absorption coefficients.  The total radiative
heat flux is assembled as

$$
q_r = \sum_{j=0}^{N_g} a_j\, q_{r,j}
$$

by solving $N_g + 1$ independent gray-gas RTEs.

### Available correlations

| Model name | Species | Correlation | Bands |
|------------|---------|-------------|-------|
| `wsgg-H2O` | H₂O only | Smith et al. (1982) | 4 gray gases + 1 clear window |
| `wsgg-H2OCO2` | H₂O + CO₂ | Fabiani et al. (JPP 2025) | 4 gray gases + 1 clear window |

The `wsgg-H2OCO2` correlation accounts for the H₂O/CO₂ molar ratio
$M_R = x_{\mathrm{H_2O}} / x_{\mathrm{CO_2}}$ via a set of discrete
reference values; the nearest $M_R$ is selected per cell.

```ini
[GROOT-Model]
model           = wsgg-H2OCO2
wall_emissivity = 0.85
```

---

## Statistical Narrow Band (SNB)

**Model names:** `snbw`, `snb`

Narrow-band models perform spectral integration over narrow bands, typically
spanning a few tens of cm⁻¹, to compute the average transmissivity of each
band.  This significantly reduces cost compared to line-by-line (LBL)
calculations while retaining the main spectral features.

!!! note "Transmissivity, not absorption coefficient"
    Unlike gray-gas or WSGG models, SNB models do not yield an absorption
    coefficient directly.  They yield the **column transmissivity** $\tau$ of
    each band, which is then used to integrate the RTE in its integral form
    (see [RTE integration](#rte-integration-for-snb) below).

The spectrum is divided into $N_b = 450$ bands of width
$\Delta\eta = 25\ \mathrm{cm}^{-1}$, covering $25$–$11250\ \mathrm{cm}^{-1}$.
The total intensity is assembled as

$$
\mathcal{I}(L) = \sum_{n=1}^{N_b} \overline{\mathcal{I}}_n(L)
$$

where $\overline{\mathcal{I}}_n$ is the band-averaged intensity, and the
Planck weights

$$
a_n(T) = \frac{\pi \overline{I}_{b,n}(T)}{\sigma T^4}\,\Delta\eta
$$

satisfy $\sum_n a_n \approx 1$.

### Spectral database

The SNB coefficients are tabulated from LBL calculations performed with
**HITEMP-2010** (H₂O, CO) and **CDSD-4000** (CO₂):

| Species | LBL database | SNB bands | $T$ range |
|---------|-------------|-----------|-----------|
| H₂O | HITEMP-2010 | 449 | 300–5000 K |
| CO₂ | CDSD-4000   | 323 | 300–5000 K |
| CO  | HITEMP-2010 | 194 | 300–5000 K |

Temperature interpolation is linear on a 100 K grid (48 points).
The SNB coefficients **do not depend explicitly on pressure** — pressure
enters only through the optical depth $u = xpL$.

The database path is baked in at CMake configure time via
`src/lib/physics/snb_datadir.f90.in`.

If a species is absent from the CFD field (e.g. no CO), its mole fraction
is set to zero automatically.

### Mean band parameters

For each band $n$ centred at $\eta_n = 25n\ \mathrm{cm}^{-1}$, the
relevant quantities are:

- $\bar{\kappa}_n$ — mean line intensity-to-spacing ratio [cm⁻¹ atm⁻¹]
- $\bar{\beta}_n$ — ratio of average Lorentz half-width to mean line spacing [-]

In GROOT these are combined with mole fractions and pressure to give a
per-cell effective absorption:

$$
\kappa_n^{\mathrm{eff}} =
  \sum_s x_s\,\frac{p}{10^5}\,\bar{K}_{s,n}(T) \times 100
  \quad [\mathrm{m}^{-1}]
$$

stored as `blk%ka(:,:,:,n)`, where $\bar{K}_{s,n}$ are the tabulated
coefficients (cm⁻¹ atm⁻¹), $p/10^5$ converts Pa to atm, and $\times 100$
converts cm⁻¹ to m⁻¹.

### SNBW — weak absorption limit

**Model name:** `snbw`

When $2\pi x p L \bar{\kappa}/\bar{\beta} \ll 1$ (weak absorption limit),
the Malkmus transmissivity (Eq. 3.15 of Fabiani 2025) simplifies to
**Beer–Lambert**:

$$
\tau_n = \exp\!\left(-\kappa_n^{\mathrm{eff}}\,\Delta s\right)
$$

SNBW shares the same `intRTE` path as WSGG and requires no additional
memory beyond `blk%ka`.

**Accuracy:** the error with respect to full SNB remains below 2% for
typical combustion-chamber conditions ($T > 2000\ \mathrm{K}$,
$p > 20\ \mathrm{bar}$) for all path lengths.  Larger errors occur at lower
pressures and temperatures, such as in the divergent section of a nozzle.

```ini
[GROOT-Model]
model           = snbw
wall_emissivity = 0.85
```

### SNB full — Malkmus model

**Model name:** `snb`

The classical SNB transmissivity for a homogeneous isothermal column of a
single species is given by the **Malkmus** formula:

$$
\tau_n = \exp\!\left[-\frac{\bar{\beta}_n}{\pi}\left(\sqrt{1 + \frac{2\pi\,x\,p\,L\,\bar{\kappa}_n}{\bar{\beta}_n}} - 1\right)\right]
$$

where $u = xpL$ is the optical path.

#### Multi-species mixtures

For mixtures with multiple radiating species, the total transmissivity is
computed as the **product** of the individual-species transmissivities:

$$
\tau_n^{\mathrm{mix}} = \prod_s \tau_{n,s}
$$

#### Non-homogeneous media: Curtis–Godson approximation

Real flames and rocket nozzles are neither homogeneous nor isothermal.
GROOT applies the **Curtis–Godson approximation**, which keeps the Malkmus
form but substitutes $\kappa$ and $\beta$ with path-averaged equivalent
parameters:

$$
\kappa_{\mathrm{eq}} = \frac{1}{u}\int_{s'}^{s} x_i\,p\,\kappa\,ds'',
\qquad
\beta_{\mathrm{eq}} = \frac{1}{u\,\kappa_{\mathrm{eq}}}\int_{s'}^{s} x_i\,p\,\kappa\,\beta\,ds''
$$

In GROOT, since each mesh cell is treated as locally homogeneous and
isothermal, the path integral reduces to a per-cell effective parameter
$\beta_{\mathrm{eff}}$ computed as the absorption-weighted average over
species:

$$
\beta_{\mathrm{eff},n} =
  \frac{\sum_s x_s\,\frac{p}{10^5}\,\bar{K}_{s,n}\,\beta_{s,n}}
       {\sum_s x_s\,\frac{p}{10^5}\,\bar{K}_{s,n}},
\qquad
\beta_{s,n} = 2\,\gamma_{s,n}\,D_{s,n}
$$

with $\gamma_{s,n}$ the Lorentz collision-broadening half-width and
$D_{s,n} = 1/\delta_{s,n}$ the inverse mean line spacing from the database.
$\beta_{\mathrm{eff},n}$ is stored as `blk%kb(i,j,k,0:Ngg)` (allocated only
for `model = snb`).

The per-cell Malkmus transmissivity then becomes:

$$
\tau_n = \exp\!\left[-\beta_{\mathrm{eff},n}\left(\sqrt{1 + \frac{2\,\kappa_n^{\mathrm{eff}}\,\Delta s}{\beta_{\mathrm{eff},n}}} - 1\right)\right]
$$

```ini
[GROOT-Model]
model           = snb
wall_emissivity = 0.85
```

### RTE integration for SNB

Because the SNB model yields transmissivity rather than an absorption
coefficient, the RTE is cast in its integral form.  The discrete form used
in GROOT (Eq. 4.14 of Fabiani 2025) is:

$$
\overline{\mathcal{I}}_n(L) =
  \overline{\mathcal{I}}_{w,n}\,\overline{\tau_n(0 \to L)}
  + \sum_{k=0}^{N}
    \Bigl(\overline{\tau_n(0 \to L_{k-1/2})} - \overline{\tau_n(0 \to L_{k+1/2})}\Bigr)
    \,\overline{\mathcal{I}}_{b,n}^k
$$

where the nodal transmissivities $\tau_n(0\to L_{k\pm 1/2})$ are obtained
as the cumulative product of per-cell transmissivities, and
$\overline{\mathcal{I}}_{w,n}$ is the wall emission weighted by
$a_n(T_{\mathrm{wall}})$.

### Computational cost

| Model | Extra memory | RTE cost vs gray |
|-------|-------------|-----------------|
| `snbw` | none beyond `blk%ka` | 450× |
| `snb`  | `blk%kb` (same size as `blk%ka`) | 450× |

The 450× factor is offset in practice by OpenMP parallelism over ray
directions.

The 450× factor is offset by OpenMP parallelism over ray directions.

---

## Spectral-Line Weighted-sum-of-gray-gases (SLW)

**Model name:** `slw`

!!! note "Optional — requires radlib"
    SLW is provided through the external **[radlib](https://github.com/BYUignite/radlib)**
    library (BYU, MIT-licensed) and is compiled only when GROOT is built with
    `--use-radlib` / `-DUSE_RADLIB=ON`. Without it, selecting `model = slw`
    aborts at runtime. See
    [Installation → Optional components](../getting-started/installation.md#optional-components).

The **Rank-Correlated SLW** (RC-SLW) model, like WSGG, represents the gas as a
set of $N_g$ gray gases plus one transparent window ($j=0$), producing weights
$a_j$ and absorption coefficients $\kappa_j$ with the **same structure as WSGG**:

$$
\varepsilon_g = \sum_{j=0}^{N_g} a_j\, \left(1 - e^{-\kappa_j\, L}\right)
$$

so GROOT reuses the existing `intRTE` gray-gas infrastructure unchanged. The
difference from WSGG lies in how $a_j$ and $\kappa_j$ are obtained: instead of
fixed polynomial correlations, RC-SLW reconstructs them **locally** for each
cell state from tabulated **Absorption-Line Blackbody Distribution Functions
(ALBDF)**, by dividing the absorption cross-section range into $N_g$ gray gases
via a rank-correlated quadrature. This makes SLW a more physically detailed
model than WSGG, at the cost of a per-cell ALBDF evaluation.

The number of gray gases is user-controlled via `slw_ngray` (radlib's $n_{GG}$,
**excluding** the transparent window, so $N_g$ = `slw_ngray`, exactly as WSGG
uses 4). More gray gases increase accuracy and cost; the solution is typically
converged by $N_g \approx 8$–$20$.

### Species and tabulated range

RC-SLW supports **H₂O, CO₂, CO** (looked up automatically — no `species =`
line needed). The ALBDF tables are valid over a limited range:

| Quantity | Valid range |
|----------|-------------|
| Temperature | 300–3000 K |
| Pressure | 0.1–50 atm |
| Gray gases (`slw_ngray`) | 1–25 |

Outside the pressure range radlib would abort, so GROOT **clamps the pressure**
to $[0.1, 50]$ atm and emits a one-time warning. Temperatures above 3000 K are
extrapolated by radlib and should be treated with caution.

### Wall emission

As for WSGG/SNB, the wall Planck fractions $a_j(T_\mathrm{wall})$ are computed
in a dedicated `Setup_SLWWall` routine at setup, using a representative
composition ($M_R = x_{\mathrm{H_2O}}/x_{\mathrm{CO_2}} = 2$, $x_\mathrm{CO}=0$)
at the reference pressure `p_ref`. Without this the wall would emit the full
$\sigma T^4$ in every gray gas.

```ini
[GROOT-Model]
model           = slw
slw_ngray       = 8         # number of gray gases (1–25)
p_ref           = 1.0e5     # reference pressure [Pa] for wall emission
wall_emissivity = 0.85
```

!!! note "SLW vs SNB/WSGG"
    RC-SLW and SNB/WSGG are independent approximations; they agree closely in
    the optically-thin limit but can differ by several percent as the optical
    thickness ($p \times L$) grows. SNB and `wsgg-H2OCO2` typically agree within
    ~2%; SLW may sit a few percent higher at elevated pressure–path.

---

## Soot

!!! warning "Not yet implemented"
    Soot radiation is **not implemented** in the current version.
    The gas-phase field may contain soot variables written by MOSE
    (controlled by `nsoot` in `[MOSE-Physics]`); GROOT reads and skips
    them correctly when indexing the field array, but the soot volume
    fraction $f_v$ is not used in any absorption coefficient calculation.
    The `fvs` argument of `compute_gasproperties` is reserved for a future
    $\kappa_\text{soot} = C_\text{abs} f_v T$ contribution.

---

## Optically Thin Approximation

When `optically_thin = .true.` (in `[GROOT-Options]`), the **volumetric
source term** $\nabla \cdot \mathbf{q}_r$ is computed without ray marching
by setting the incident irradiation $G = 0$:

$$
\nabla \cdot \mathbf{q}_r \approx 4\pi\,\kappa\, I_b
  \quad \text{(instead of } \kappa(4\pi I_b - G)\text{)}
$$

This is valid when $\kappa L \ll 1$ so that re-absorption is negligible
and each cell only emits without receiving radiation from other cells.

!!! note
    `optically_thin` affects **only `co_source`** (the volumetric source
    term).  The wall heat flux computed by `co_flux` always uses full DTM
    ray marching regardless of this flag.

---

## References

- M.F. Modest, *The Weighted-Sum-of-Gray-Gases Model for Arbitrary Solution
  Methods in Radiative Transfer*, Journal of Heat Transfer, 113(3):650–656,
  1991. (correlation used in `wsgg-H2O`)
- M. Fabiani, M.T. Migliorino, D. Bianchi, and F. Nasuti, *Spectral and
  Global Radiative Heat Transfer Models for Liquid Propellant Rocket Engines*,
  Journal of Propulsion and Power, 41:650–664, 2025.
  DOI: [10.2514/1.B39892](https://doi.org/10.2514/1.B39892)
  (correlation used in `wsgg-H2OCO2`)
- A. Soufiani and J. Taine, *High temperature gas radiative property
  parameters of statistical narrow-band model for H₂O, CO₂ and CO, and
  correlated-K model for H₂O and CO₂*, International Journal of Heat and
  Mass Transfer, 40(4):987–991, 1997.
- P. Rivière and A. Soufiani, *Updated band model parameters for H₂O, CO₂,
  CH₄ and CO radiation at high temperature*, International Journal of Heat
  and Mass Transfer, 55:3349–3358, 2012.
- M.F. Modest, *Radiative Heat Transfer*, 3rd edition, Academic Press,
  2013.
- V.P. Solovjov, D. Lemonnier, and B.W. Webb, *The Rank Correlated SLW model
  of gas radiation in non-uniform media*, Journal of Quantitative Spectroscopy
  and Radiative Transfer, 137:141–161, 2014. (RC-SLW model, `slw`)
- D.O. Lignell et al., *radlib* — radiation property library,
  [github.com/BYUignite/radlib](https://github.com/BYUignite/radlib)
  (implementation used by `slw`).
