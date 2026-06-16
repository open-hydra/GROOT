# Theoretical Guide

This section provides the theoretical foundations of the physical models and
numerical methods implemented in GROOT.

<div class="grid cards" markdown>

-   :material-math-integral:{ .lg .middle } __Radiative Transfer Equation__

    ---

    The governing integro-differential equation for radiative intensity,
    emission, absorption, and scattering in a participating medium

    [:octicons-arrow-right-24: Radiative transfer equation](rte.md)

-   :material-ray-start:{ .lg .middle } __Discrete Transfer Method__

    ---

    Ray-tracing discretization of the RTE: ray sampling, intensity
    integration along a ray, wall boundary conditions, and convergence

    [:octicons-arrow-right-24: Discrete transfer method](dtm.md)

-   :material-waves:{ .lg .middle } __Spectral Models__

    ---

    Gray gas (constant $\kappa$), weighted-sum-of-gray-gases (WSGG)
    approximations and SNB models for non-gray participating media

    [:octicons-arrow-right-24: Spectral models](spectral-models.md)

</div>

---

## Overview

GROOT solves the radiative transfer equation (RTE) by the Discrete
Transfer Method (DTM), a ray-tracing approach in which a finite number
of rays are traced from each wall face through the participating medium.
The numerical pipeline can be summarised as follows:

| Stage | Method | Page |
|-------|--------|------|
| **Governing equation** | Radiative transfer equation (RTE) | [RTE](rte.md) |
| **Spatial discretization** | Discrete Transfer Method (DTM) | [DTM](dtm.md) |
| **Spectral integration** | Gray gas, WSGG or SNB | [Spectral Models](spectral-models.md) |

---
