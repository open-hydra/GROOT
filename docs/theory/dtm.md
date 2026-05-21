# Discrete Transfer Method

This page describes the **Discrete Transfer Method (DTM)** as implemented
in GROOT.  The DTM is a ray-tracing scheme that directly discretises the
formal solution of the RTE along a finite set of ray directions.

---

## Algorithm Overview

The DTM, originally proposed by Lockwood and Shah (1981), proceeds as
follows:

1. **Discretise the unit sphere** into $N_r$ discrete directions
   $\{\hat{\mathbf{s}}_m, w_m\}$, where $w_m$ is the solid-angle weight.
2. For each ray direction $\hat{\mathbf{s}}_m$:
   a. **Launch rays** from every wall face, travelling through the medium.
   b. **Integrate intensity** along each ray using the formal solution of the RTE.
   c. **Accumulate** the resulting intensity at the originating wall face.
3. **Compute wall heat flux** $Q$ and radiative source $\nabla \cdot \mathbf{q}_r$
   from the accumulated intensities.

```mermaid
flowchart LR
    A["Wall face<br/>(emission point)"] -->|"ray direction s_m"| B["Cell 1"]
    B --> C["Cell 2"]
    C --> D["..."]
    D --> E["Opposing wall<br/>(termination point)"]
    E -->|"reflected/transmitted intensity"| D
    D --> C
    C --> B
    B -->|"accumulated I_m"| A

    style A fill:#263238,stroke:#90a4ae,color:#eceff1
    style E fill:#263238,stroke:#90a4ae,color:#eceff1
```

---

## Ray Direction Sampling

GROOT samples ray directions uniformly over the full sphere $4\pi$ using
a deterministic quadrature.  The number of rays $N_r$ is set by the
parameter `Nr` in `[GROOT-Discretization]`.

For a **2-D axisymmetric** configuration (`twoDax = .true.`) the domain is
a thin wedge of angular extent $\delta\theta$ (set by `ang` in degrees).
The azimuthal integration is performed analytically, reducing the problem
to a 2-D integration in the meridional plane.

---

## Intensity Integration Along a Ray

Consider a ray travelling from point $P_0$ (a wall face centre) in
direction $\hat{\mathbf{s}}$ through a sequence of $N_c$ cells.  Cell $k$
has path length $\Delta s_k$, absorption coefficient $\kappa_k$, and
blackbody source $I_{b,k} = \sigma T_k^4 / \pi$.

The intensity leaving cell $k$ is obtained by integrating the RTE over
$\Delta s_k$:

$$
I_k = I_{k-1}\, e^{-\kappa_k \Delta s_k}
    + I_{b,k}\, (1 - e^{-\kappa_k \Delta s_k})
$$

Starting from the wall intensity $I_0 = I_w$ (determined by the wall BC),
this recurrence is marched from the terminating wall back to the emitting
wall $P_0$.

---

## Wall Boundary Conditions

### Opaque diffuse wall

The intensity leaving the wall in any outgoing direction is:

$$
I_w = \varepsilon\, I_b(T_w) + (1 - \varepsilon)\, q_{\text{inc}} / \pi
$$

where $q_{\text{inc}}$ is the hemispherically integrated incident flux
accumulated from all incoming rays.

### Specular wall (`bc = 300`)

The ray is reflected according to Snell's law:

$$
\hat{\mathbf{s}}' = \hat{\mathbf{s}} - 2 (\hat{\mathbf{s}} \cdot \hat{\mathbf{n}})\, \hat{\mathbf{n}}
$$

No energy is emitted or absorbed. Used for symmetry planes (e.g. the
$z$-faces of a 2-D slab, or the axial faces of a 2-D axisymmetric wedge).

### Azimuthal / periodic wall (`bc = 200`)

The ray is re-entered at the corresponding periodic face.  Used for the
azimuthal faces of the 2-D axisymmetric wedge.

---

## Wall Heat Flux

The net heat flux at a wall face $f$ is assembled from all ray
contributions:

$$
Q_f = \sum_{m=1}^{N_r} I_m(\hat{\mathbf{s}}_m) \cdot |\hat{\mathbf{n}}_f \cdot \hat{\mathbf{s}}_m| \cdot w_m
    - \varepsilon_f \sigma T_{w,f}^4
$$

Stored as `face(f)%Qtot` [W/m²]; positive values indicate energy
deposited on the wall (net absorption).

---

## Radiative Source Term

The volumetric source term in cell $c$ is:

$$
\left(\nabla \cdot \mathbf{q}_r\right)_c
= 4 \kappa_c \sigma T_c^4
  - \kappa_c \sum_{m=1}^{N_r} I_{m,c} \cdot w_m
$$

Stored as `block%source` [W/m³].

---

## Convergence

The DTM iteration converges when the relative change in wall heat flux
between successive sweeps falls below the tolerance `tol` in
`[GROOT-Discretization]`.  For non-scattering problems (`scatter = .false.`)
a single pass is sufficient (set `itermax = 1`) because there is no
coupling between ray directions.

---

## References

- F.C. Lockwood and N.G. Shah, *A new radiation solution method for
  incorporation in general combustion prediction procedures*, 18th
  Symposium (Int.) on Combustion, pp. 1405–1414, 1981.
