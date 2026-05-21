#!/usr/bin/env python3
"""
Verification script for the symmetry (periodic BC) cube test.

The half-cube [0, 0.5]^3 with symmetry BCs on three faces must reproduce
the radiation field of the corresponding octant of the full 1x1x1 cube.
Reference data is therefore shared with the homcube test (data/ symlinks to
homcube/data or files are copied there).

Checks
------
  heatflux_face1.dat   face-1 (x=0) wall heat flux vs z     → ref_heatflux_k{k}_T1000.dat
  source_corner.dat    source at (x≈0.483, y≈0.483) vs z    → ref_source_k{k}_T_1000.dat
                       (corner ≈ mirror of centerline, same by symmetry)

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

# ── CLI ────────────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(description="Verify periodic-BC cube test")
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
def load_dat(path):
    return np.loadtxt(path, comments="#")

ref_hf  = load_dat(f"data/ref_heatflux_k{k_lbl}_T1000.dat")
ref_src = load_dat(f"data/ref_source_k{k_lbl}_T_1000.dat")
out_hf  = load_dat(f"{output_dir}/heatflux_face1.dat")
out_src = load_dat(f"{output_dir}/source_corner.dat")

z_ref_hf,  q_ref  = ref_hf[:,0],  ref_hf[:,1]
z_ref_src, s_ref  = ref_src[:,0], ref_src[:,1]
z_hf,      q_out  = out_hf[:,0],  out_hf[:,1]
z_src,     s_out  = out_src[:,0], out_src[:,1]

# The half-cube domain is [0, 0.5] m in x, y, z.
# Keep only the ref points that fall within the half-domain for comparison.
mask_hf  = z_ref_hf  <= 0.5 + 1e-9
mask_src = z_ref_src <= 0.5 + 1e-9
q_ref_i = np.interp(z_hf,  z_ref_hf[mask_hf],   q_ref[mask_hf])
s_ref_i = np.interp(z_src, z_ref_src[mask_src],  s_ref[mask_src])

# ── Plot ───────────────────────────────────────────────────────────────────
fig, axs = plt.subplots(1, 2, figsize=(10, 4))
fig.suptitle(f"Periodic-BC cube  —  k = {k_lbl} m⁻¹,  T = 1000 K")

axs[0].plot(z_ref_hf[mask_hf],  q_ref[mask_hf], "k--", label="Reference (homcube)")
axs[0].plot(z_hf,               q_out,           "r-",  label="GROOT")
axs[0].set_xlabel("z  [m]")
axs[0].set_ylabel("Wall heat flux  [W m⁻²]")
axs[0].set_title("Face-1 heat flux")
axs[0].legend()
axs[0].grid(True, alpha=0.3)

axs[1].plot(z_ref_src[mask_src], s_ref[mask_src], "k--", label="Reference (homcube CL)")
axs[1].plot(z_src,               s_out,            "r-",  label="GROOT (corner)")
axs[1].set_xlabel("z  [m]")
axs[1].set_ylabel("Source term  [W m⁻³]")
axs[1].set_title("Near-corner source  (x = y ≈ 0.483 m)")
axs[1].legend()
axs[1].grid(True, alpha=0.3)

plt.tight_layout()
svg_path = os.path.join(output_dir, f"{folder_name}_k{k_lbl}.svg")
plt.savefig(svg_path, bbox_inches="tight", transparent=True)
print(f"Plot saved: {svg_path}")

if args.plot:
    plt.show()

# ── Error check ────────────────────────────────────────────────────────────
# Slightly looser threshold: the corner is not exactly the mirror of the
# centreline due to angular-quadrature asymmetry.
THRESHOLD = 10.0   # max relative error [%]

err_hf  = np.max(np.abs((q_out - q_ref_i) / (np.abs(q_ref_i) + 1e-30))) * 100.0
err_src = np.max(np.abs((s_out - s_ref_i) / (np.abs(s_ref_i) + 1e-30))) * 100.0

print(f"Max relative error — heat flux : {err_hf:.2f} %  (threshold {THRESHOLD} %)")
print(f"Max relative error — source    : {err_src:.2f} %  (threshold {THRESHOLD} %)")

failed = False
if err_hf > THRESHOLD:
    print(f"FAIL: heat flux error {err_hf:.2f}% exceeds {THRESHOLD}%")
    failed = True
if err_src > THRESHOLD:
    print(f"FAIL: source error {err_src:.2f}% exceeds {THRESHOLD}%")
    failed = True

if failed:
    sys.exit(1)

print("PASS")
