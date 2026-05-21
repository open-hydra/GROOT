#!/usr/bin/env python3
"""
Analytical reference solution for the axisymmetric finite cylinder test.

Geometry: radius R=1 m, full length 2 m (half-length L=1 m).
          Axis = x; z in [-L, L] → output coordinate x = L+z in [0, 2].
Medium:   uniform T=3500 K, absorption coefficient k [m⁻¹], no scattering.
Walls:    black (eps=1), cold (Tw=0 K) on endcaps and lateral wall.

Replaces very_old_groot/test/algorithm/finite_cylinder/data/sol_exact.m.

Usage:
  python3 cylinder.py           # compute for k = 0.1, 1.0, 10.0
  python3 cylinder.py 1.0       # single k value

Outputs (written to this script's directory):
  ref_qwall_k{k}.dat            lateral wall heat flux   [W/m²]  (100 pts)
  ref_source_axis_k{k}.dat      axis source term         [W/m³]  (100 pts)
"""

import sys
import os
import numpy as np
from scipy.integrate import dblquad, quad

# ── Physical and geometry constants ─────────────────────────────────────────
sigma = 5.67037e-8   # Stefan-Boltzmann  [W/m²/K⁴]
T     = 3500.0       # Gas temperature   [K]
R     = 1.0          # Cylinder radius   [m]
L     = 1.0          # Half-length       [m]  (full cylinder: 0 to 2L)
npts  = 100          # Output points (matches old sol_exact.m)

Ib = sigma * T**4 / np.pi   # Blackbody intensity [W/m²/sr]

# ── Lateral wall heat flux at r=R ───────────────────────────────────────────
def compute_qwall(k_abs):
    """
    Compute heat flux on the lateral surface (r=R) at npts axial locations.

    Three contributions for each point z ∈ (-L, L):
      q1 : far  endcap (at z = +L), theta in [themax(phi), pi]
      q2 : near endcap (at z = -L), theta in [0, themin(phi)]
      q3 : lateral wall,            theta in [themin(phi), themax(phi)]

    dblquad(func(theta, phi), phi_lo, phi_hi, theta_lo, theta_hi)
    """
    r = R
    # Exclude exact endpoints to avoid 0/0 in path lengths
    z_arr = np.linspace(-L, L, npts + 2)[1:-1]
    q     = np.zeros(npts)

    for idx, zi in enumerate(z_arr):
        d_up   = L - zi    # distance to far  endcap (> 0)
        d_down = L + zi    # distance to near endcap (> 0)

        def chord(phi):
            c2 = max(R**2 - r**2 * np.sin(phi)**2, 0.0)
            return r * np.cos(phi) + np.sqrt(c2)

        def themin_f(phi):
            return np.arctan(chord(phi) / d_down)

        def themax_f(phi):
            return np.pi - np.arctan(chord(phi) / d_up)

        # Far endcap: cos(the) < 0  →  path = d_up / (-cos(the))
        def fun1(the, phi):
            return (1.0 - np.exp(-k_abs * d_up / (-np.cos(the)))) \
                   * np.sin(the)**2 * np.cos(phi)

        # Near endcap: cos(the) > 0  →  path = d_down / cos(the)
        def fun2(the, phi):
            return (1.0 - np.exp(-k_abs * d_down / np.cos(the))) \
                   * np.sin(the)**2 * np.cos(phi)

        # Lateral wall: path = chord(phi) / sin(the)
        def fun3(the, phi):
            path = chord(phi) / np.sin(the)
            return (1.0 - np.exp(-k_abs * path)) * np.sin(the)**2 * np.cos(phi)

        q1, _ = dblquad(fun1, -np.pi/2, np.pi/2,
                        themax_f, np.pi,    epsabs=1e-5, epsrel=1e-4)
        q2, _ = dblquad(fun2, -np.pi/2, np.pi/2,
                        0.0,      themin_f, epsabs=1e-5, epsrel=1e-4)
        q3, _ = dblquad(fun3, -np.pi/2, np.pi/2,
                        themin_f, themax_f, epsabs=1e-5, epsrel=1e-4)

        q[idx] = Ib * (q1 + q2 + q3)

    return L + z_arr, q


# ── Axis source term at r=0 ─────────────────────────────────────────────────
def compute_source_axis(k_abs):
    """
    Compute the radiative source S = k*(4*sigma*T^4 - G) on the axis (r=0).

    At r=0 the theta limits are phi-independent → 1-D quad over theta only,
    then multiply by the trivial phi integral = 2*pi.  Much faster than dblquad.
    """
    z_arr = np.linspace(-L, L, npts + 2)[1:-1]
    src   = np.zeros(npts)

    for idx, zi in enumerate(z_arr):
        d_up   = L - zi
        d_down = L + zi

        # Theta limits at r=0 (phi-independent)
        themin_v = np.arctan(R / d_down)
        themax_v = np.pi - np.arctan(R / d_up)

        # Irradiance integrand (sin(the) factor only; phi gives 2*pi)
        def g1(the): return (1.0 - np.exp(-k_abs * d_up   / (-np.cos(the)))) * np.sin(the)
        def g2(the): return (1.0 - np.exp(-k_abs * d_down /   np.cos(the)))  * np.sin(the)
        def g3(the): return (1.0 - np.exp(-k_abs * R      /   np.sin(the)))  * np.sin(the)

        G1, _ = quad(g1, themax_v, np.pi,     limit=200)
        G2, _ = quad(g2, 0.0,      themin_v,  limit=200)
        G3, _ = quad(g3, themin_v, themax_v,  limit=200)

        G      = Ib * 2.0 * np.pi * (G1 + G2 + G3)
        src[idx] = k_abs * (4.0 * sigma * T**4 - G)

    return L + z_arr, src


# ── I/O ─────────────────────────────────────────────────────────────────────
def write_dat(fname, x, y, header):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(script_dir, fname)
    with open(path, 'w') as f:
        f.write(header + '\n')
        for xi, yi in zip(x, y):
            f.write(f'{xi:20.8E} {yi:20.8E}\n')
    print(f'  Written {path}')


# ── Main ─────────────────────────────────────────────────────────────────────
if __name__ == '__main__':
    k_vals = [0.1, 1.0, 10.0] if len(sys.argv) < 2 else [float(sys.argv[1])]

    for k in k_vals:
        kstr = f'{k:g}'
        x, q = compute_qwall(k)
        write_dat(f'ref_qwall_k{kstr}.dat',
                  x, q, '# x [m]              qwall [W/m2]')
        x, s = compute_source_axis(k)
        write_dat(f'ref_source_axis_k{kstr}.dat',
                  x, s, '# x [m]              source_axis [W/m3]')
