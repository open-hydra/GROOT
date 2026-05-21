#!/usr/bin/env python3
"""
Verification script for the scalene-trapezoid (quadrilateral) test.

Reference data in data/ comes from a legacy GROOT version and stores
heat flux normalised by sigma*T^4 (dimensionless).  GROOT outputs
physical heat flux in W/m².

The reference is indicative only (the old code reportedly segfaulted on
this case), so a generous 20 % tolerance is used.

Checks
------
  heatflux_face3.dat   bottom-wall (face 3) heat flux vs x  → ref_k{k}.dat

Usage
-----
  python verify.py              # k=1.0, no plot window
  python verify.py --k 0.1
  python verify.py --k 10 --plot
"""

import argparse
import sys
import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

# ── Matplotlib: transparent SVG, theme-aware text ──────────────────────────
mpl.rcParams.update({
    "figure.facecolor": "none",
    "axes.facecolor":   "none",
    "savefig.facecolor":"none",
    "svg.fonttype":     "none",
})

# ── Physical constants ─────────────────────────────────────────────────────
SIGMA = 5.670374419e-8   # Stefan-Boltzmann  [W m⁻² K⁻⁴]
T_GAS = 3000.0           # Uniform gas temperature [K]
# Reference is normalised by sigma*T^4
NORM  = SIGMA * T_GAS**4

# ── CLI ────────────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(description="Verify quadrilateral test")
parser.add_argument("--k",    type=float, default=1.0,
                    choices=[0.1, 1.0, 10.0],
                    help="Absorption coefficient used in the run (default: 1.0)")
parser.add_argument("--plot", action="store_true",
                    help="Show interactive plot window")
args = parser.parse_args()

k     = args.k
k_lbl = f"{k:g}"
output_dir  = "OUTPUT"
folder_name = os.path.basename(os.getcwd())

# ── Load data ──────────────────────────────────────────────────────────────
# Reference: plain space-separated, no header, Q is dimensionless (Q/sigma*T^4)
ref_raw = np.loadtxt(f"data/ref_k{k_lbl}.dat")
x_ref   = ref_raw[:,0]
q_ref   = ref_raw[:,1] * NORM       # convert to W/m²

# GROOT output: comment header, Q in W/m²
out_raw = np.loadtxt(f"{output_dir}/heatflux_face3.dat", comments="#")
x_out   = out_raw[:,0]
q_out   = out_raw[:,1]

# Interpolate GROOT output at reference x positions
q_out_at_ref = np.interp(x_ref, x_out, q_out)

# ── Plot ───────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(7, 4))
ax.plot(x_ref, q_ref,           "k--", label="Reference (legacy, normalised)")
ax.plot(x_out, q_out,           "r-",  label="GROOT")
ax.set_xlabel("x  [m]")
ax.set_ylabel("Wall heat flux  [W m⁻²]")
ax.set_title(f"Scalene trapezoid — face-3 heat flux  (k = {k_lbl} m⁻¹,  T = {T_GAS:.0f} K)")
ax.legend()
ax.grid(True, alpha=0.3)

plt.tight_layout()
svg_path = os.path.join(output_dir, f"{folder_name}_k{k_lbl}.svg")
plt.savefig(svg_path, bbox_inches="tight", transparent=True)
print(f"Plot saved: {svg_path}")

if args.plot:
    plt.show()

# ── Error check ────────────────────────────────────────────────────────────
# Generous threshold: reference is from legacy code (indicative only).
THRESHOLD = 20.0   # max relative error [%]

err = np.max(np.abs((q_out_at_ref - q_ref) / (np.abs(q_ref) + 1e-30))) * 100.0

print(f"Max relative error — face-3 heat flux : {err:.2f} %  (threshold {THRESHOLD} %)")

if err > THRESHOLD:
    print(f"FAIL: heat flux error {err:.2f}% exceeds {THRESHOLD}%")
    sys.exit(1)

print("PASS")
