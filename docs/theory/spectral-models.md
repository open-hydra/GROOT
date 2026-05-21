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

## Weighted-Sum-of-Gray-Gases (WSGG)

**Model names:** `wsgg`, `wsgg-H2O`, `wsgg-H2OCO2`

The **WSGG** model approximates the non-gray (spectrally dependent)
absorption of real participating gases by a weighted sum of $N_g$ gray
gases plus one transparent "window" band:

$$
\varepsilon_g = \sum_{j=0}^{N_g} a_j(T)\, \left(1 - e^{-\kappa_j(\eta_i, p_i) L}\right)
$$

where:

- $a_j(T)$ are temperature-dependent weighting factors (WSGG weights),
- $\kappa_j$ are the gray-gas absorption coefficients for band $j$,
- $L$ is the path length.

The WSGG model allows the total radiative heat transfer to be obtained as a
sum of $N_g + 1$ gray-gas RTE solutions, each solved with a different
$\kappa_j$, and assembled with the weights $a_j$:

$$
q_r = \sum_{j=0}^{N_g} a_j\, q_{r,j}
$$

### Available correlations

GROOT uses the WSGG correlations of **Smith et al. (1982)** for H₂O and
H₂O/CO₂ mixtures:

| Model name | Species | Bands |
|------------|---------|-------|
| `wsgg-H2O` | H₂O only | 4 gray gases + 1 clear window |
| `wsgg-H2OCO2` | H₂O + CO₂ | 4 gray gases + 1 clear window |

The mixture is automatically selected based on the `species` list in
`[GROOT-Model]`:

- If `CO2` is listed → `wsgg-H2OCO2`
- Otherwise → `wsgg-H2O`

```ini
[GROOT-Model]
model   = wsgg
species = H2O,CO2
wall_emissivity = 0.85
```

### Spectral database

The internal species database contains:

| Species | Role |
|---------|------|
| H₂O | Primary radiating species; included in both WSGG sets |
| CO₂ | Secondary radiating species; triggers H₂O/CO₂ WSGG |
| CO  | Minor radiating species |
| CH₄ | Minor radiating species |

Species are mapped from the CFD solver field via `phase.txt`. Only
species listed in `species` (in `[GROOT-Model]`) are treated as radiating.

---

## Soot

When `soot-field = .true.`, GROOT reads the soot volume fraction $f_v$
from the gas-phase field (two extra variables appended after the main
thermodynamic variables) and adds the soot contribution to $\kappa$:

$$
\kappa_{\text{soot}} = C_{\text{abs}} f_v T
$$

where $C_{\text{abs}}$ is a model constant.  The soot model can be
combined with any gas-phase spectral model.

---

## Optically Thin Approximation

When `optically_thin = .true.` (in `[GROOT-Options]`), GROOT ignores
the absorption term along each ray:

$$
I(L) \approx I(0) + \int_0^L \kappa\, I_b\, \mathrm{d}s
$$

This is valid only when $\kappa L \ll 1$ throughout the domain. The wall
heat flux is then dominated by direct emission rather than re-absorption.

---

## References

- T.F. Smith, Z.F. Shen, and J.N. Friedman, *Evaluation of coefficients
  for the weighted sum of gray gases model*, Journal of Heat Transfer,
  104(4):602–608, 1982.
- M.F. Modest, *Radiative Heat Transfer*, 3rd edition, Academic Press,
  2013.
