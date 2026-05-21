#!/usr/bin/env python3
"""
Homogeneous cube reference solution.

Geometry : L x L x L cube, black cold walls (T_wall = 0, eps = 1).
Medium   : uniform temperature T, absorption coeff k (no scattering).
Method   : Monte-Carlo-style angular quadrature (same as DTM angle set).

Generates
  source_k{k}_T_{T}.dat   -- z-axis centerline volumetric source [W/m³] vs z
  heatflux_k{k}_T{T}.dat  -- face-1 (x=0) wall heat flux [W/m²] vs z at y=L/2

Usage
  python cube.py                         # default: L=1, T=1000, Nr=256, Nz=30
  python cube.py L T Nr Nz              # custom parameters
"""

import numpy as np
import os
import sys

SIGMA = 5.670374419e-8   # Stefan-Boltzmann constant [W/m²/K⁴]


def ray_cube_pathlength(x0, y0, z0, sx, sy, sz, L):
    """Pathlength from (x0,y0,z0) in direction (sx,sy,sz) to the first face
    of the cube [0,L]^3.  Returns 0.0 if the starting point is on a boundary
    and the ray immediately exits (e.g. corner point with outward component);
    in that case I = Ib*(1-exp(0)) = 0 so the flux contribution is zero."""
    t_min = np.inf
    eps   = 1e-12 * L

    faces = [
        (0.0, 0, x0, 1, y0, 2, z0),   # x = 0
        (L,   0, x0, 1, y0, 2, z0),   # x = L
        (0.0, 1, y0, 0, x0, 2, z0),   # y = 0
        (L,   1, y0, 0, x0, 2, z0),   # y = L
        (0.0, 2, z0, 0, x0, 1, y0),   # z = 0
        (L,   2, z0, 0, x0, 1, y0),   # z = L
    ]
    s_vec = (sx, sy, sz)
    p_vec = (x0, y0, z0)

    for fval, nidx, n0, a1, pa1, a2, pa2 in faces:
        sn = s_vec[nidx]
        if abs(sn) < eps:
            continue
        t = (fval - n0) / sn
        if t <= eps:
            continue
        p1 = pa1 + s_vec[a1] * t
        p2 = pa2 + s_vec[a2] * t
        if -eps <= p1 <= L + eps and -eps <= p2 <= L + eps:
            if t < t_min:
                t_min = t

    return t_min if t_min < np.inf else 0.0


def compute_source(L, k, T, Nr, Nz, filename):
    """Centerline volumetric radiative source term div(q) [W/m³] vs z.

    Angular quadrature identical to DTM source integration:
      nphi = 2*sqrt(Nr/2),  nthe = sqrt(Nr/2)  (full-sphere, same as lib_rte)
      the in [-pi/2, pi/2], phi in [0, 2*pi]
      direction: s = (sin(the), cos(the)*cos(phi), cos(the)*sin(phi))
    """
    Ib   = SIGMA * T**4 / np.pi
    nphi = int(2 * np.sqrt(Nr / 2.0))
    nthe = int(    np.sqrt(Nr / 2.0))
    dphi = 2.0 * np.pi / nphi
    dthe = np.sin(np.pi / nthe)       # matches GROOT lib_rte.f90 angle set

    x0 = L / 2.0
    y0 = L / 2.0
    z0vec = np.linspace(0.01 * L, (1.0 - 0.01) * L, Nz)

    rows = []
    for z0 in z0vec:
        G = 0.0
        for iphi in range(1, nphi + 1):
            phi = np.pi / nphi + (iphi - 1) * dphi
            for ithe in range(1, nthe + 1):
                the = np.pi / (2.0 * nthe) + (ithe - 1) * dthe - np.pi / 2.0
                sx = np.sin(the)
                sy = np.cos(the) * np.cos(phi)
                sz = np.cos(the) * np.sin(phi)
                t  = ray_cube_pathlength(x0, y0, z0, sx, sy, sz, L)
                I  = Ib * (1.0 - np.exp(-k * t))
                # solid angle weight: cos(the) = sin(pi/2 - the)
                G += I * np.cos(the) * dthe * dphi
        div_q = k * (4.0 * SIGMA * T**4 - G)
        rows.append((z0, div_q))

    with open(filename, 'w') as fid:
        fid.write('# z [m]              source [W/m3]\n')
        for z, s in rows:
            fid.write(f'{z:20.8E} {s:20.8E}\n')
    print(f'Written {filename}  ({len(rows)} points)')


def compute_heatflux(L, k, T, Nr, Nz, filename):
    """Wall heat flux on face x=0 at y=L/2 vs z [W/m²].

    Angular quadrature identical to DTM wall integration:
      nphi = 2*sqrt(Nr),  nthe = 0.5*sqrt(Nr)  (hemisphere, same as lib_rte)
      the in [0, pi/2] measured from the face plane (sx = sin(the))
      direction: s = (sin(the), cos(the)*cos(phi), cos(the)*sin(phi))
    """
    Ib     = SIGMA * T**4 / np.pi
    nphi   = int(2.0 * np.sqrt(Nr))
    nthe   = int(0.5 * np.sqrt(Nr))
    dtheta = (np.pi / 2.0) / nthe
    dphi   = 2.0 * np.pi / nphi

    x0 = 0.0
    y0 = L / 2.0
    z0vec = np.linspace(0.0, L, Nz)

    rows = []
    for z0 in z0vec:
        Q = 0.0
        for iphi in range(1, nphi + 1):
            phi = np.pi / nphi + (iphi - 1) * dphi
            for ithe in range(1, nthe + 1):
                the = np.pi / (4.0 * nthe) + (ithe - 1) * dtheta
                sx  = np.sin(the)               # > 0: ray enters domain from x=0
                sy  = np.cos(the) * np.cos(phi)
                sz  = np.cos(the) * np.sin(phi)
                t   = ray_cube_pathlength(x0, y0, z0, sx, sy, sz, L)
                I = Ib * (1.0 - np.exp(-k * t))
                # Lambert weight: cos(theta_from_normal) = sin(the)
                Q += I * np.cos(the) * np.sin(the) * dtheta * dphi
        rows.append((z0, Q))

    with open(filename, 'w') as fid:
        fid.write('# z [m]              Q_wall [W/m2]\n')
        for z, q in rows:
            fid.write(f'{z:20.8E} {q:20.8E}\n')
    print(f'Written {filename}  ({len(rows)} points)')


def k_label(k):
    """Format k for filename: 1.0->'1', 0.1->'0.1', 10.0->'10'."""
    s = f'{k:g}'
    return s


if __name__ == '__main__':
    args = sys.argv[1:]
    L  = float(args[0]) if len(args) > 0 else 1.0
    T  = float(args[1]) if len(args) > 1 else 1000.0
    Nr = int(args[2])   if len(args) > 2 else 256
    Nz = int(args[3])   if len(args) > 3 else 30

    os.makedirs('data', exist_ok=True)

    T_int = int(T)
    for k in [0.1, 1.0, 10.0]:
        kl = k_label(k)
        compute_source(
            L, k, T, Nr, Nz,
            f'data/ref_source_k{kl}_T_{T_int}.dat'
        )
        compute_heatflux(
            L, k, T, Nr, Nz,
            f'data/ref_heatflux_k{kl}_T{T_int}.dat'
        )
