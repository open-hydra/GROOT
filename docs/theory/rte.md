# Radiative Transfer Equation

This page presents the governing equation solved by GROOT: the **radiative
transfer equation (RTE)** for a gray or spectrally participating medium.

---

## Radiative Intensity

The primary unknown is the **spectral radiative intensity**
$I_\eta(\mathbf{r}, \hat{\mathbf{s}})$ [W/(m²·sr·m⁻¹)], which describes the
energy per unit time, per unit area normal to $\hat{\mathbf{s}}$, per unit
solid angle around $\hat{\mathbf{s}}$, and per unit wavenumber interval, at
position $\mathbf{r}$ travelling in direction $\hat{\mathbf{s}}$.

For a **gray** medium (spectrally averaged) the subscript $\eta$ is dropped
and $I$ [W/(m²·sr)] denotes the total intensity.

---

## Governing Equation

In a participating medium with absorption coefficient $\kappa$ [m⁻¹],
scattering coefficient $\sigma_s$ [m⁻¹], and blackbody intensity
$I_b = \sigma T^4 / \pi$ [W/(m²·sr)], the steady-state RTE along direction
$\hat{\mathbf{s}}$ reads:

$$
\hat{\mathbf{s}} \cdot \nabla I(\mathbf{r}, \hat{\mathbf{s}})
= \kappa\, I_b(T) - (\kappa + \sigma_s)\, I(\mathbf{r}, \hat{\mathbf{s}})
  + \frac{\sigma_s}{4\pi} \int_{4\pi} \Phi(\hat{\mathbf{s}}', \hat{\mathbf{s}})\,
    I(\mathbf{r}, \hat{\mathbf{s}}')\, \mathrm{d}\Omega'
$$

The three terms on the right-hand side are:

1. **Emission** — energy emitted by the gas at the local temperature $T$.
2. **Extinction** — energy removed from direction $\hat{\mathbf{s}}$ by absorption and out-scattering.
3. **In-scattering** — energy redirected into $\hat{\mathbf{s}}$ from all other directions $\hat{\mathbf{s}}'$, weighted by the phase function $\Phi$.

In the current version of GROOT scattering is not implemented ($\sigma_s = 0$), so the RTE reduces to:

$$
\hat{\mathbf{s}} \cdot \nabla I
= \kappa\, I_b - \kappa\, I
$$

---

## Formal Solution

Along a ray path $s$ from an origin point $s = 0$ to $s = L$, the formal
solution of the non-scattering RTE is:

$$
I(L) = I(0)\, e^{-\kappa L}
  + \int_0^L \kappa\, I_b(s')\, e^{-\kappa(L - s')}\, \mathrm{d}s'
$$

The first term is the attenuated wall contribution; the second is the gas
emission collected along the path.

---

## Wall Boundary Condition

At an **opaque diffuse wall** with emissivity $\varepsilon$ and temperature
$T_w$, the outgoing intensity in any direction $\hat{\mathbf{s}}$ pointing
away from the wall is:

$$
I_w(\hat{\mathbf{s}}) = \varepsilon\, I_b(T_w)
  + \frac{1 - \varepsilon}{\pi} \int_{\hat{\mathbf{n}} \cdot \hat{\mathbf{s}}' < 0}
    I(\hat{\mathbf{s}}')\, |\hat{\mathbf{n}} \cdot \hat{\mathbf{s}}'|\, \mathrm{d}\Omega'
$$

where $\hat{\mathbf{n}}$ is the inward wall normal.  For a **black wall**
($\varepsilon = 1$) this reduces to $I_w = I_b(T_w)$.

A **specular wall** (used for symmetry planes, `bc = 3`) simply reflects the
incoming ray: $I_w(\hat{\mathbf{s}}) = I(\hat{\mathbf{s}}^*)$, where
$\hat{\mathbf{s}}^*$ is the mirror direction.

---

## Derived Quantities

### Radiative heat flux

The **net radiative heat flux vector** at a point is:

$$
\mathbf{q}_r = \int_{4\pi} I(\hat{\mathbf{s}})\, \hat{\mathbf{s}}\, \mathrm{d}\Omega
$$

The scalar **wall heat flux** $Q$ [W/m²] is the component of $\mathbf{q}_r$
normal to the wall surface, computed by GROOT and stored in `face(f)%Qtot`.

### Irradiation

The **irradiation** (incident radiation) is the total intensity integrated
over all directions:

$$
G = \int_{4\pi} I(\hat{\mathbf{s}})\, \mathrm{d}\Omega \quad [\text{W/m}^2]
$$

### Radiative source term

The **divergence of the radiative heat flux** (radiative source/sink per unit
volume) is:

$$
\nabla \cdot \mathbf{q}_r = \kappa\, (4\pi I_b - G)
  = 4\kappa\sigma T^4 - \kappa G \quad [\text{W/m}^3]
$$

This quantity (stored in `block%source`) appears as a source term in the
energy equation of the host CFD solver.

---
