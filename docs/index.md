---
title: GROOT
hide:
  - navigation
  - toc
---

<div class="groot-hero" markdown>

<p class="groot-hero__eyebrow">HYDRA CFD SUITE</p>

<h1 class="groot-hero__title">GROOT</h1>

<p class="groot-hero__lead">Open-source solver of the radiative transfer equation on multi-block structured grids, by
the discrete transfer method, for gray and spectral gases.</p>

<p class="groot-hero__cmd"><span>$</span> ./install.sh build --compilers=gnu --use-openmp</p>

[Get started](getting-started/index.md){ .md-button .md-button--primary }
[Theory](theory/index.md){ .md-button .groot-hero__ghost }

<p class="groot-hero__name"><strong>G</strong>rid based <strong>R</strong>adiati<strong>O</strong>n
f<strong>O</strong>rtran <strong>T</strong>oolbox</p>

</div>

<p class="groot-kicker">The solver</p>

## Radiative heat transfer on your mesh { .groot-heading }

<p class="groot-sub" markdown>GROOT gives the radiative heat flux at the walls and the source term in the gas. Start with
the [Quick Start](getting-started/quick-start.md): a shipped case runs in a few minutes.</p>

<div class="groot-tools" markdown>

[:material-white-balance-sunny: <b>Discrete transfer method</b> <i>Steady radiative transfer by ray tracing, with the heat flux and the source term as results.</i>](theory/dtm.md){ .groot-tool }
[:material-cube-outline: <b>1D, 2D and 3D</b> <i>Multi-block structured grids, with periodic boundaries.</i>](user/using.md){ .groot-tool }
[:material-chart-bell-curve: <b>Spectral gas models</b> <i>Gray gas, weighted sum of gray gases (WSGG) and SNB, for the wavelength-dependent emission and absorption.</i>](theory/spectral-models.md){ .groot-tool }
[:material-chip: <b>OpenMP</b> <i>Shared-memory parallel runs, on a workstation or on a cluster node.</i>](getting-started/installation.md){ .groot-tool }
[:material-link-variant: <b>Coupling</b> <i>Used as a library by the other Hydra solvers, to give them the radiative flux and source.</i>](development/structure.md){ .groot-tool }
[:material-file-chart-outline: <b>Inputs and outputs</b> <i>Solution files in Tecplot (ASCII and binary) and VTK.</i>](user/using.md#output){ .groot-tool }

</div>

<p class="groot-kicker"></p>

## Hydra CFD Suite { .groot-heading }

<p class="groot-sub">Hydra couples solvers of gases, real fluids, solids, particles and radiation in one
simulation. GROOT is its radiation solver.</p>

<div class="hydra" markdown>

<div class="hydra__case" markdown>

<p class="hydra__label">Your case</p>

<p class="hydra__file" markdown>:material-grid: structured mesh</p>

<p class="hydra__file" markdown>:material-database: [databases](https://open-hydra.github.io/ATLAS/databases/)</p>

</div>

<div class="hydra__atlas" markdown>

<p class="hydra__name">ATLAS</p>

<p class="hydra__tools" markdown>[GPB](https://open-hydra.github.io/ATLAS/user-guide/gpb/) [BCB](https://open-hydra.github.io/ATLAS/user-guide/bcb/) [ICB](https://open-hydra.github.io/ATLAS/user-guide/icb/) [STB](https://open-hydra.github.io/ATLAS/user-guide/stb/) [MDB](https://open-hydra.github.io/ATLAS/user-guide/mdb/) [KAnT](https://open-hydra.github.io/ATLAS/user-guide/kant/)</p>

</div>

<div class="hydra__link">
<span class="hydra__link-line"></span>
<span class="hydra__link-bottom">input files for each run</span>
</div>

<div class="hydra__suite" markdown>

<p class="hydra__label" markdown>[Hydra](https://github.com/open-hydra) · the coupled simulation</p>

<div class="hydra__solvers" markdown>

[<b>MOSE</b> compressible flow <i>ideal-gas</i>](https://github.com/open-hydra/MOSE){ .hydra__solver .hydra__solver--gas }
[<b>GROOT</b> radiative transfer <i>no ATLAS phase</i>](https://github.com/open-hydra/GROOT){ .hydra__solver .hydra__solver--radiation .is-here }
[<b>ARES</b> real-fluid flow <i>real-fluid</i>](https://github.com/open-hydra/ARES){ .hydra__solver .hydra__solver--fluid }
[<b>FUSS</b> heat conduction in solids <i>solid</i>](https://github.com/open-hydra/FUSS){ .hydra__solver .hydra__solver--solid }
[<b>IGLOO</b> Lagrangian particles <i>dispersed</i>](https://github.com/open-hydra/IGLOO){ .hydra__solver .hydra__solver--particle }
[<b>ICE</b> Eulerian particles <i>dispersed</i>](https://github.com/open-hydra/ICE){ .hydra__solver .hydra__solver--particle }

</div>

</div>

<div class="hydra__libs" markdown>

<p class="hydra__label">Built on</p>

[<b>ORION</b> Tecplot, VTK and PLOT3D I/O](https://MarcoGrossi92.github.io/ORION/){ .hydra__lib }
[<b>FLINT</b> thermochemistry](https://MarcoGrossi92.github.io/FLINT/){ .hydra__lib }
[<b>OSLO</b> ODE integrators](https://github.com/MarcoGrossi92/OSLO){ .hydra__lib }
[<b>FiNeR</b> INI file parser](https://github.com/szaghi/FiNeR){ .hydra__lib }

</div>

</div>

<div class="groot-links" markdown>

[Quick Start](getting-started/quick-start.md)
[User Guide](user/index.md)
[Input Parameters](user/input/input-parameters.md)
[Theory](theory/index.md)
[Verification](vv/index.md)
[Developer Guide](development/index.md)

</div>
