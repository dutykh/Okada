#!/usr/bin/env python3
# Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology, Abu Dhabi, UAE)
"""Figure of the README: the surface displacements computed by OkadaSol.m for a sample thrust fault.

Run okada_readme_data.m first (GNU Octave or MATLAB), then, from this folder,
    python3 plot_readme_figure.py
which writes okada-displacement.png. The style is the stamped copy manuscript_style.py beside it.
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize
from PIL import Image

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import manuscript_style as hs  # noqa: E402

# The fault of okada_readme_data.m: its surface projection and the section of panel (b).
LENGTH, WIDTH, DIP = 60.0, 40.0, np.radians(13.0)       # km, km, radians
SECTION = 30.0                                           # km, the middle of the fault along strike

x = np.loadtxt(HERE / "x.csv") / 1e3
y = np.loadtxt(HERE / "y.csv") / 1e3
v = np.loadtxt(HERE / "v.csv", delimiter=",")
w = np.loadtxt(HERE / "w.csv", delimiter=",")             # w[i, j]: at (x[j], y[i])
level = np.ceil(100 * np.abs(w).max()) / 100
top_edge = WIDTH * np.cos(DIP)

hs.configure(family="libertine")
plt.rcParams["savefig.dpi"] = 260
plt.rcParams["figure.facecolor"] = hs.PAPER
plt.rcParams["savefig.facecolor"] = hs.PAPER

fig = plt.figure(figsize=(hs.WIDTH, 3.0))
grid = fig.add_gridspec(1, 3, width_ratios=[1.0, 0.04, 1.05])
left = fig.add_subplot(grid[0, 0])
bar_axis = fig.add_subplot(grid[0, 1])
right = fig.add_subplot(grid[0, 2])

mesh = left.pcolormesh(x, y, w, shading="gouraud", cmap=hs.SIGNED_COLORMAP,
                       norm=Normalize(-level, level), rasterized=True)
left.plot([0, LENGTH, LENGTH, 0, 0], [0, 0, top_edge, top_edge, 0], color=hs.INK, lw=0.8,
          linestyle=(0, (5.0, 2.0)))
left.axvline(SECTION, color=hs.SLATE, lw=0.6, linestyle=(0, (1.6, 1.6)))
left.set_aspect("equal")
left.set_xlim(x.min(), x.max())
left.set_ylim(y.min(), y.max())
left.set_xlabel(r"$x$ (km)")
left.set_ylabel(r"$y$ (km)")
hs.panel(left, "(a) vertical displacement")
bar = fig.colorbar(mesh, cax=bar_axis, label=r"$w$ (m)")
bar.outline.set_linewidth(0.5)
bar.outline.set_edgecolor(hs.INK)
bar.ax.tick_params(width=0.5, length=2.5)

j = int(np.argmin(np.abs(x - SECTION)))
right.axvspan(0.0, top_edge, color=hs.SAND, zorder=0, gid="background projection of the fault")
right.axhline(0.0, color=hs.SLATE, lw=0.5, linestyle=(0, (1.6, 1.6)))
for n, (field, label) in enumerate([(w[:, j], r"vertical, $w$"), (v[:, j], r"horizontal across the strike, $v$")]):
    style = hs.series(n)
    style["marker"] = None
    right.plot(y, field, lw=1.1, label=label, **style)
right.set_xlim(y.min(), y.max())
right.set_xlabel(r"$y$ (km)")
right.set_ylabel("displacement (m)")
hs.panel(right, rf"(b) section $x = {SECTION:g}$ km")
hs.legend_outside(right, where="below", ncol=1)

out = HERE / "okada-displacement.png"
hs.save(fig, out)
# The background is opaque: without its alpha channel, and compressed a second time without
# loss, the file is smaller.
Image.open(out).convert("RGB").save(out, optimize=True)
print("written", out, f"(colour scale from -{level:g} to {level:g} m)")
