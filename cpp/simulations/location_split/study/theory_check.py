"""Compare the school arms with the one-step theory of the note "State Transitions in the MEmilio-ABM".

In the ABM, a susceptible in a room with N persons present, i of them infectious, is infected within one
hour with probability p = 1 - exp(-a i / (N - 1)), a = Lambda * dt (the note's Eq. (9) with N = s + i in its
homogeneous setting; here N also counts exposed and recovered persons). The expected number of new
infections in the room is g = s p. The runs write, per school room and hour, N, s, i and the new infections
(--write_school_rooms). This script

1. fits the single parameter a by maximum likelihood over all room hours of all school arms,
2. compares observed and predicted new infections (one-step check),
3. shows that the new infections per infectious pupil hour depend only on the susceptible share
   s / (N - 1) of the room, for every room size (to first order a * s / (N - 1), Proposition 1 and Sec. 5.3),
4. follows the susceptible share that infectious pupils see over the days, for fixed and redrawn rooms,
5. compares the outcomes of fixed and daily redrawn rooms, paired with the no-split arm.

Writes results/theory_check.csv and results/figures/theory_*.png.
"""
import glob
import math
import os

import h5py
import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt

STUDY_DIR = os.environ.get("STUDY_DIR", "study_output")
OUT_DIR = os.environ.get("STUDY_OUT", os.path.join(STUDY_DIR, "results"))
FIG = os.path.join(OUT_DIR, "figures")
os.makedirs(FIG, exist_ok=True)
N_PERSONS = 10000
SIZES = [50, 25, 10, 5]
BASE = "sw1_no_split"

INK, PAPER, MUTED, GRID = "#16233A", "#F6F5F1", "#5B6578", "#DADDE3"
FIXED, REDRAWN = "#2F6297", "#D4692A"
SIZE_COLORS = {0: INK, 50: "#9DB4D3", 25: "#6C8EBF", 10: "#2F6297", 5: "#D4692A"}
plt.rcParams.update({"font.size": 14, "font.family": "DejaVu Sans", "axes.spines.top": False,
                     "axes.spines.right": False, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
                     "xtick.color": MUTED, "ytick.color": MUTED, "text.color": INK, "axes.titlesize": 15,
                     "axes.titleweight": "bold", "axes.titlelocation": "left", "axes.grid": True,
                     "grid.color": GRID, "grid.linewidth": 0.8, "figure.facecolor": PAPER,
                     "axes.facecolor": PAPER, "savefig.facecolor": PAPER, "legend.frameon": False,
                     "figure.dpi": 160, "lines.linewidth": 2.4})


def arm_name(size, redrawn=False):
    if size == 0:
        return BASE
    return f"sw1_{'redrawn_' if redrawn else ''}rooms_{size}"


def arms():
    found = [(0, False)] + [(s, r) for r in (False, True) for s in SIZES]
    return [(s, r) for s, r in found if os.path.isdir(os.path.join(STUDY_DIR, arm_name(s, r), "runs"))]


def room_table(arm):
    files = sorted(glob.glob(os.path.join(STUDY_DIR, arm, "runs", "run_*_school_rooms.csv")))
    if not files:
        raise SystemExit(f"{arm}: no school room tables; run it with --write_school_rooms")
    t = pd.concat(pd.read_csv(f) for f in files)
    return t.groupby(["Day", "Persons", "Susceptible", "Infectious"], as_index=False)[["Room_Hours", "New_Infections"]].sum(), len(files)


# ---------------------------------------------------------------------------------------------------------
# 1. Fit a = Lambda * dt over all room hours with at least one susceptible and one infectious present.
# ---------------------------------------------------------------------------------------------------------
tables = {(s, r): room_table(arm_name(s, r))[0] for s, r in arms()}
cells = pd.concat(tables.values())
cells = cells[(cells.Susceptible > 0) & (cells.Infectious > 0) & (cells.Persons > 1)]
cells = cells.groupby(["Persons", "Susceptible", "Infectious"], as_index=False)[["Room_Hours", "New_Infections"]].sum()
x = (cells.Infectious / (cells.Persons - 1)).values  # i / (N - 1)
trials = (cells.Susceptible * cells.Room_Hours).values  # susceptible hours
events = cells.New_Infections.values


def fit_a(x, trials, events):
    """Maximum likelihood a for P(infection of one susceptible in one hour) = 1 - exp(-a x)."""
    def log_likelihood(a):
        p = 1 - np.exp(-a * x)
        return np.sum(events * np.log(np.maximum(p, 1e-300)) - (trials - events) * a * x)

    grid = np.exp(np.linspace(math.log(1e-4), math.log(10), 400))
    a = grid[np.argmax([log_likelihood(v) for v in grid])]
    for _ in range(60):  # golden section refinement around the grid maximum
        lo, hi = a / 1.05, a * 1.05
        m1, m2 = lo + (hi - lo) * 0.382, lo + (hi - lo) * 0.618
        a = (lo + m2) / 2 if log_likelihood(m1) > log_likelihood(m2) else (m1 + hi) / 2
    return a


a = fit_a(x, trials, events)
# The same fit per arm: a lower a in small rooms is a one-step difference that the homogeneous theory does not
# have (e.g. age mixing, viral load over time); the rest of the room effect has to come from the depletion.
a_per_arm = {}
for key, t in tables.items():
    c = t[(t.Susceptible > 0) & (t.Infectious > 0) & (t.Persons > 1)]
    a_per_arm[key] = fit_a((c.Infectious / (c.Persons - 1)).values, (c.Susceptible * c.Room_Hours).values,
                           c.New_Infections.values)
    print(f"  a for {arm_name(*key)}: {a_per_arm[key]:.4f} ({(a_per_arm[key] / a_per_arm.get((0, False), a) - 1) * 100:+.1f}% vs no split)")
cells["Predicted"] = cells.Susceptible * (1 - np.exp(-a * x)) * cells.Room_Hours
print(f"fitted a = Lambda*dt = {a:.4f} per hour; z = a i/(N-1) <= {a:.3f} for one infectious, so p ~ z (Sec. 5)")
print(f"room hours with s>0, i>0: {int(cells.Room_Hours.sum())}, new infections observed {int(events.sum())}, "
      f"predicted {cells.Predicted.sum():.0f}")

# ---------------------------------------------------------------------------------------------------------
# 2. Observed vs predicted, per arm, binned by the predicted new infections per room hour.
# 3. New infections per infectious pupil hour against the susceptible share s/(N-1).
# ---------------------------------------------------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(13, 5))
share_bins = np.linspace(0, 1, 11)
for (size, redrawn), t in tables.items():
    if redrawn:
        continue
    c = t[(t.Susceptible > 0) & (t.Infectious > 0) & (t.Persons > 1)].copy()
    cx = c.Infectious / (c.Persons - 1)
    c["Pred"] = c.Susceptible * (1 - np.exp(-a * cx)) * c.Room_Hours
    c["PredPerHour"] = c.Pred / c.Room_Hours
    bins = np.quantile(c.PredPerHour, np.linspace(0, 1, 9))
    c["bin"] = np.digitize(c.PredPerHour, np.unique(bins)[1:-1])
    g = c.groupby("bin")[["Pred", "New_Infections", "Room_Hours"]].sum()
    label = "no split" if size == 0 else f"fixed rooms of {size}"
    axes[0].plot(g.Pred / g.Room_Hours, g.New_Infections / g.Room_Hours, "o-", ms=6, color=SIZE_COLORS[size], label=label)

    c["share"] = c.Susceptible / (c.Persons - 1)
    c["bin2"] = np.digitize(c.share, share_bins[1:-1])
    h = c.groupby("bin2").agg(new=("New_Infections", "sum"),
                              inf_hours=("Room_Hours", lambda r: (r * c.loc[r.index, "Infectious"]).sum()),
                              share=("share", lambda v: np.average(v, weights=c.loc[v.index, "Room_Hours"] * c.loc[v.index, "Infectious"])))
    h = h[h.inf_hours > 200]
    axes[1].plot(h.share, h.new / h.inf_hours, "o-", ms=6, color=SIZE_COLORS[size], label=label)
lim = axes[0].get_xlim()[1]
axes[0].plot([0, lim], [0, lim], color=MUTED, lw=1, ls="--", label="theory, Eq. (9)")
axes[0].set(xlabel="predicted new infections per room hour", ylabel="observed",
            title="One step: ABM rooms close to g(s, i)")
axes[1].plot([0, 1], [0, a], color=MUTED, lw=1, ls="--", label=f"a · s/(N−1), a = {a:.3f}")
axes[1].set(xlabel="susceptible share of the room, s/(N−1)", ylabel="new infections per infectious pupil hour",
            title="Output per infectious pupil vs susceptible share")
axes[0].legend(fontsize=11); axes[1].legend(fontsize=11)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "theory_1_one_step.png")); plt.close(fig)

# ---------------------------------------------------------------------------------------------------------
# 4. Susceptible share seen by infectious pupils, per day (weighted by infectious pupil hours).
# ---------------------------------------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(9, 5))
seen = {}
for (size, redrawn), t in tables.items():
    if size not in (0, 5, 10):
        continue
    c = t[(t.Infectious > 0) & (t.Persons > 1)]
    w = c.Infectious * c.Room_Hours
    daily = (c.Susceptible / (c.Persons - 1) * w).groupby(c.Day).sum() / w.groupby(c.Day).sum()
    seen[(size, redrawn)] = daily
    smooth = daily.rolling(7, min_periods=1).mean()
    label = "no split" if size == 0 else f"{'redrawn' if redrawn else 'fixed'} rooms of {size}"
    ax.plot(smooth.index, smooth.values, color=SIZE_COLORS[size], ls="--" if redrawn else "-", label=label)
ax.set(xlabel="day", ylabel="susceptible share around infectious pupils",
       title="Fixed small rooms run out of susceptibles", ylim=(0, 1))
ax.legend(fontsize=11)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "theory_2_depletion.png")); plt.close(fig)

# ---------------------------------------------------------------------------------------------------------
# 5. Outcomes, paired with the no-split arm: school and all infections by day 30, delay to 20 %, peak.
# ---------------------------------------------------------------------------------------------------------
def events_of(arm, run):
    e = pd.read_csv(os.path.join(STUDY_DIR, arm, "runs", f"run_{run}_infection_events.csv"))
    e = e[e.Location_Type != "Initial"]
    return (e.Timestep - 1) // 24, e.Location_Type == "School"


def curve(arm, run):
    with h5py.File(os.path.join(STUDY_DIR, arm, "total_infections", f"results_run{run}.h5"), "r") as h:
        return h["0/Time"][()], h["0/Total"][()].ravel()


def day_reaching(t, y, fraction):
    above = np.nonzero(y >= fraction * N_PERSONS)[0]
    return t[above[0]] if len(above) else np.nan


def peak(t, y):
    daily = np.diff(np.interp(np.arange(0, int(t[-1]) + 1), t, y))
    return pd.Series(daily).rolling(7, center=True, min_periods=1).mean().max()


def ci(values):
    v = np.asarray(values, float)
    m, se = np.nanmean(v), np.nanstd(v, ddof=1) / math.sqrt(np.sum(~np.isnan(v)))
    return m, m - 1.96 * se, m + 1.96 * se


runs = min(len(glob.glob(os.path.join(STUDY_DIR, arm_name(s, r), "total_infections", "results_run*.h5"))) for s, r in arms())
base = []
for run in range(runs):
    day, school = events_of(BASE, run)
    t, y = curve(BASE, run)
    base.append(dict(school30=(school & (day < 30)).sum(), all30=(day < 30).sum(), d20=day_reaching(t, y, 0.2),
                     peak=peak(t, y), final=y[-1]))
rows = []
for size, redrawn in arms():
    if size == 0:
        continue
    arm = arm_name(size, redrawn)
    school30, all30, delay, peaks, final = [], [], [], [], []
    for run in range(runs):
        day, school = events_of(arm, run)
        t, y = curve(arm, run)
        b = base[run]
        school30.append(((school & (day < 30)).sum() / b["school30"] - 1) * 100)
        all30.append(((day < 30).sum() / b["all30"] - 1) * 100)
        delay.append(day_reaching(t, y, 0.2) - b["d20"])
        peaks.append((peak(t, y) / b["peak"] - 1) * 100)
        final.append((y[-1] - b["final"]) / N_PERSONS * 100)
    row = dict(arm=arm, room_size=size, rooms="redrawn" if redrawn else "fixed", runs=runs)
    for name, values in [("school_infections_day30_pct", school30), ("all_infections_day30_pct", all30),
                         ("delay_20pct_days", delay), ("peak_change_pct", peaks), ("final_size_diff_pp", final)]:
        row[name], row[name + "_ci_low"], row[name + "_ci_high"] = [round(v, 2) for v in ci(values)]
    row["fitted_a_arm"] = round(a_per_arm[(size, redrawn)], 5)
    row["fitted_a_vs_no_split_pct"] = round((a_per_arm[(size, redrawn)] / a_per_arm[(0, False)] - 1) * 100, 1)
    sh = seen.get((size, redrawn))
    if sh is not None:
        row["susceptible_share_seen_days0_30"] = round(float(sh[sh.index < 30].mean()), 3)
    rows.append(row)
result = pd.DataFrame(rows).sort_values(["rooms", "room_size"], ascending=[True, False])
result["fitted_a"] = round(a, 5)
result.to_csv(os.path.join(OUT_DIR, "theory_check.csv"), index=False)
print(result[["arm", "school_infections_day30_pct", "all_infections_day30_pct", "delay_20pct_days",
              "peak_change_pct", "final_size_diff_pp"]].to_string(index=False))

fig, axes = plt.subplots(1, 3, figsize=(14, 4.6))
pos = {50: 0, 25: 1, 10: 2, 5: 3}
for rooms, color, offset in [("fixed", FIXED, -0.08), ("redrawn", REDRAWN, 0.08)]:
    r = result[result.rooms == rooms]
    xs = np.array([pos[s] for s in r.room_size]) + offset
    for ax, col in zip(axes, ["school_infections_day30_pct", "delay_20pct_days", "peak_change_pct"]):
        ax.errorbar(xs, r[col], yerr=[r[col] - r[col + "_ci_low"], r[col + "_ci_high"] - r[col]], marker="o", ms=7,
                    capsize=3, color=color, label=f"{rooms} rooms")
for ax, title, ylabel in zip(axes, ["School infections by day 30", "Delay to 20 % infected", "Peak of daily incidence"],
                             ["change vs no split (%)", "days later than no split", "change vs no split (%)"]):
    ax.axhline(0, color=MUTED, lw=0.8)
    ax.set(title=title, ylabel=ylabel, xlabel="pupils per room")
    ax.set_xticks(range(4)); ax.set_xticklabels(["50", "25", "10", "5"])
axes[0].legend(fontsize=11)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "theory_3_fixed_vs_redrawn.png")); plt.close(fig)
print("figures written to", FIG)
