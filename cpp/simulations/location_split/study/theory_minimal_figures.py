"""Figures for the minimal theory experiments (location_split_theory), simulated mean against the closed forms.

Usage: python theory_minimal_figures.py <directory with the csv files of location_split_theory> <figure directory>
"""
import os
import sys

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt

DATA = sys.argv[1] if len(sys.argv) > 1 else "study_output/theory_minimal"
OUT = sys.argv[2] if len(sys.argv) > 2 else os.path.join(DATA, "figures")
os.makedirs(OUT, exist_ok=True)

INK, PAPER, MUTED, GRID = "#16233A", "#F6F5F1", "#5B6578", "#DADDE3"
A_COLORS = {0.05: "#9DB4D3", 0.5: "#6C8EBF", 1.0: "#2F6297", 2.0: INK}
ROOM_COLORS = {200: INK, 50: "#9DB4D3", 25: "#6C8EBF", 10: "#2F6297", 5: "#D4692A"}
plt.rcParams.update({"font.size": 14, "font.family": "DejaVu Sans", "axes.spines.top": False,
                     "axes.spines.right": False, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
                     "xtick.color": MUTED, "ytick.color": MUTED, "text.color": INK, "axes.titlesize": 15,
                     "axes.titleweight": "bold", "axes.titlelocation": "left", "axes.grid": True,
                     "grid.color": GRID, "grid.linewidth": 0.8, "figure.facecolor": PAPER,
                     "axes.facecolor": PAPER, "savefig.facecolor": PAPER, "legend.frameon": False,
                     "figure.dpi": 160, "lines.linewidth": 2.2})

# 1. One step: Proposition 1 (left) and Proposition 2 (right), both relative to the linear value.
single = pd.read_csv(os.path.join(DATA, "one_step_single.csv"))
clus = pd.read_csv(os.path.join(DATA, "one_step_clustering.csv"))
fig, axes = plt.subplots(1, 2, figsize=(13, 5))
ax = axes[0]
for a, g in single.groupby("a"):
    ax.plot(g.susceptible, g.closed_form / a, color=A_COLORS[a], label=f"a = {a:g}")
    ax.errorbar(g.susceptible, g.observed_per_room / a, yerr=2 * g.se / a, fmt="o", ms=6, color=A_COLORS[a],
                capsize=3, lw=1.2)
ax.set_xscale("log")
ax.set(xlabel="susceptibles s in the room", ylabel="new infections per step / a",
       title="Prop. 1: one infectious agent, one step", ylim=(0, 1.15))
ax.legend(fontsize=11, title="lines: g(s, 1); dots: ABM", title_fontsize=11, loc="lower right")
ax = axes[1]
for a, g in clus.groupby("a"):
    for arrangement, ls, marker in [("split", "-", "o"), ("clustered", "--", "s")]:
        h = g[g.arrangement == arrangement]
        ax.plot(h.x, h.closed_form / (a * h.x), color=A_COLORS[a], ls=ls,
                label=f"a = {a:g}, {arrangement}" if a in (0.05, 1.0) else None)
        ax.errorbar(h.x, h.observed_per_group / (a * h.x), yerr=2 * h.se / (a * h.x), fmt=marker, ms=6,
                    color=A_COLORS[a], capsize=3, lw=1.2)
ax.set(xlabel="x infectious and x susceptible agents", ylabel="new infections per step / (a · x)",
       title="Prop. 2: one room (dashed) vs x rooms (solid)", ylim=(0, 1.15))
ax.legend(fontsize=11, loc="lower left")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "minimal_1_one_step.png"))
plt.close(fig)

# 2. Contact time: cumulative new infections per school of 200 with 10 infectious agents.
ct = pd.read_csv(os.path.join(DATA, "contact_time.csv"))
fig, axes = plt.subplots(1, 2, figsize=(13, 5))
for ax, a in zip(axes, sorted(ct.a.unique())):
    d = ct[ct.a == a]
    for (arrangement, size), g in d.groupby(["arrangement", "room_size"], sort=False):
        g = g.sort_values("hour")
        color = ROOM_COLORS[size]
        if arrangement == "redrawn":
            label = f"redrawn every 8 h, rooms of {size}"
            ax.plot(g.hour, g.expected_abm_cumulative, color=color, ls=":", lw=2.4, label=label)
        else:
            label = "one room of 200" if arrangement == "one room" else f"fixed rooms of {size}"
            ax.plot(g.hour, g.closed_form, color=color, label=label)
        pts = g[g.hour % 16 == 0]
        ax.plot(pts.hour, pts.observed_cumulative, "o", ms=5, color=color)
    ax.set(xlabel="hours together", ylabel="new infections per school",
           title=f"a = {a:g} per hour" + (" (as in the city ABM)" if a < 0.1 else ""))
    ax.set_xlim(0, 160)
axes[0].legend(fontsize=10.5, loc="upper left")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "minimal_2_contact_time.png"))
plt.close(fig)
print("figures written to", OUT)
