"""Figures for the location split presentation, from the finished arms in $STUDY_DIR."""
import os
STUDY_DIR = os.environ.get("STUDY_DIR", "study_output")
OUT_DIR = os.environ.get("STUDY_OUT", os.path.join(STUDY_DIR, "results"))
BINARY = os.environ.get("SIM_BINARY", "build/bin/location_split_simulation")
import glob, json, os
import h5py, numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = STUDY_DIR
OUT = os.path.join(OUT_DIR, "figures")
os.makedirs(OUT, exist_ok=True)
N = 10000
table = pd.read_csv(os.path.join(OUT_DIR, "paired_results_day90.csv"))
# Deck palette: ink, paper, blue ramp for room sizes, orange for the smallest rooms.
INK, PAPER, MUTED, GRID = "#16233A", "#F6F5F1", "#5B6578", "#DADDE3"
COLORS = {"no split": INK, "50": "#9DB4D3", "25": "#6C8EBF", "10": "#2F6297", "5": "#D4692A", "3": "#8A3B12"}
plt.rcParams.update({"font.size": 15, "font.family": "DejaVu Sans", "axes.spines.top": False, "axes.spines.right": False,
                     "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED,
                     "text.color": INK, "axes.titlesize": 17, "axes.titleweight": "bold", "axes.titlelocation": "left",
                     "axes.grid": True, "axes.grid.axis": "y", "grid.color": GRID, "grid.linewidth": 0.8,
                     "figure.facecolor": PAPER, "axes.facecolor": PAPER, "savefig.facecolor": PAPER,
                     "legend.frameon": False, "figure.dpi": 160, "lines.linewidth": 2.6})


def has(arm):
    return arm in set(table["arm"])


def median_curve(arm):
    curves = []
    for path in glob.glob(f"{BASE}/{arm}/total_infections/results_run*.h5"):
        with h5py.File(path, "r") as h:
            t, y = h["0/Time"][()], h["0/Total"][()].ravel()
        curves.append(np.interp(np.arange(0, 91), t, y))
    curves = np.array(curves) / N * 100
    return np.median(curves, axis=0), np.percentile(curves, 25, axis=0), np.percentile(curves, 75, axis=0)


def label_of(arm):
    return "no split" if arm.endswith("no_split") else arm.split("_")[-1]


# 1. Epidemic curves of the school split (school weight 1).
arms = [a for a in ["sw1_no_split", "sw1_rooms_50", "sw1_rooms_25", "sw1_rooms_10", "sw1_rooms_5"] if has(a)]
fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))
for arm in arms:
    m, lo, hi = median_curve(arm)
    lab = label_of(arm)
    axes[0].plot(np.arange(91), m, color=COLORS[lab], lw=2, label=("rooms of " + lab) if lab != "no split" else lab)
    if lab == "no split":
        axes[0].fill_between(np.arange(91), lo, hi, color=COLORS[lab], alpha=0.12)
    daily = np.diff(m) * N / 100
    axes[1].plot(np.arange(1, 91), pd.Series(daily).rolling(7, center=True, min_periods=1).mean(), color=COLORS[lab], lw=2)
axes[0].set(xlabel="day", ylabel="infected (% of city)", title="Cumulative infections")
axes[1].set(xlabel="day", ylabel="new infections per day (7-day mean)", title="Daily incidence")
axes[0].legend(frameon=False)
fig.tight_layout(); fig.savefig(f"{OUT}/01_school_split_curves.png"); plt.close(fig)

# 2. Three outcomes vs room size, per school share group.
groups = [("sw1", "school share 8%"), ("share15", "school share 15%"), ("share25", "school share 25%")]
fig, axes = plt.subplots(1, 3, figsize=(13, 4.2))
for g, (prefix, name) in enumerate(groups):
    rows = table[table["arm"].str.startswith(prefix + "_")]
    if rows.empty:
        continue
    sizes = []
    for _, r in rows.iterrows():
        sizes.append(0 if r["arm"].endswith("no_split") else int(r["arm"].split("_")[-1]))
    rows = rows.assign(rooms=sizes).sort_values("rooms", key=lambda s: s.replace(0, 10**6), ascending=False)
    position = {0: 0, 50: 1, 25: 2, 10: 3, 5: 4}
    x = np.array([position[s] for s in rows["rooms"]]) + (g - 1) * 0.08
    gc = ["#2F6297", "#D4692A", "#8A3B12"][g]
    axes[0].errorbar(x, rows["diff_pp"], yerr=[rows["diff_pp"] - rows["ci95_low"], rows["ci95_high"] - rows["diff_pp"]],
                     marker="o", ms=8, capsize=4, color=gc, label=name)
    axes[1].plot(x, rows["delay_20pct_days"], marker="o", ms=8, color=gc, label=name)
    axes[2].plot(x, rows["peak_change_pct"], marker="o", ms=8, color=gc, label=name)
for ax, title, ylabel in zip(axes, ["Infected by day 90", "Delay to 20% infected", "Peak of daily incidence"],
                             ["change vs no split (pp, 95% CI)", "days later than no split", "change vs no split (%)"]):
    ax.axhline(0, color="grey", lw=0.8)
    ax.set(title=title, ylabel=ylabel, xlabel="pupils per room")
    ax.set_xticks(range(5)); ax.set_xticklabels(["no split", "50", "25", "10", "5"])
if len(table[table["arm"].str.startswith("share")]):
    axes[2].legend(fontsize=12)
fig.tight_layout(); fig.savefig(f"{OUT}/02_school_split_effects.png"); plt.close(fig)

# 3. Where infections happen.
settings = ["Home", "SocialEvent", "Work", "School", "BasicsShop"]
names = {"Home": "home", "SocialEvent": "social events", "Work": "work", "School": "school", "BasicsShop": "shops"}
setcol = [INK, "#D4692A", "#2F6297", "#6C8EBF", "#B9C0CC"]
order = [a for a in ["sw1_no_split", "sw1_rooms_25", "sw1_rooms_5", "work_rooms_5", "work_rooms_3", "event_rooms_5",
                     "event_rooms_3", "share25_no_split", "share25_rooms_5"] if has(a)]
pretty = {"sw1_no_split": "no split", "sw1_rooms_25": "school 25", "sw1_rooms_5": "school 5", "work_rooms_5": "work 5",
          "work_rooms_3": "work 3", "event_rooms_5": "events 5", "event_rooms_3": "events 3",
          "share25_no_split": "no split\n(school 25%)", "share25_rooms_5": "school 5\n(school 25%)"}
fig, ax = plt.subplots(figsize=(12, 4.4))
bottom = np.zeros(len(order))
for s, c in zip(settings, setcol):
    vals = np.array([table.loc[table["arm"] == a, f"share_{s}_pct"].values[0] for a in order])
    ax.bar(range(len(order)), vals, bottom=bottom, color=c, label=names[s])
    for i, v in enumerate(vals):
        if v > 5:
            ax.text(i, bottom[i] + v / 2, f"{v:.0f}", ha="center", va="center", color="white", fontsize=9)
    bottom += vals
ax.set_xticks(range(len(order))); ax.set_xticklabels([pretty[a] for a in order], fontsize=10)
ax.set(ylabel="share of infections (%)", title="")
ax.legend(frameon=False, ncol=5, loc="upper center", bbox_to_anchor=(0.5, -0.18))
fig.tight_layout(); fig.savefig(f"{OUT}/03_infection_settings.png"); plt.close(fig)

# 4. Workplace and social event splits.
we = [a for a in ["work_rooms_5", "work_rooms_3", "event_rooms_5", "event_rooms_3", "sw1_rooms_5"] if has(a)]
if we:
    rows = table.set_index("arm").loc[we]
    lab = [pretty[a] for a in we]
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2))
    axes[0].bar(lab, rows["diff_pp"], yerr=[rows["diff_pp"] - rows["ci95_low"], rows["ci95_high"] - rows["diff_pp"]],
                color="#3b6ea5", capsize=3)
    axes[1].bar(lab, rows["delay_20pct_days"], color="#e09f3e")
    axes[2].bar(lab, rows["peak_change_pct"], color="#c8553d")
    for ax, title, ylabel in zip(axes, ["Infected by day 90", "Delay to 20% infected", "Peak of daily incidence"],
                                 ["change vs no split (pp)", "days", "change (%)"]):
        ax.axhline(0, color="grey", lw=0.8); ax.set(title=title, ylabel=ylabel); ax.tick_params(axis="x", labelsize=9)
    fig.tight_layout(); fig.savefig(f"{OUT}/04_work_event_split.png"); plt.close(fig)

# 5. Calibration.
calib = json.load(open(os.environ.get("CALIBRATION_RESULT", "calibration_output/fit_result.json")))
first, last = calib["history"][0]["shares"], calib["history"][-1]["shares"]
keys = ["home", "social_event", "work"]
fig, ax = plt.subplots(figsize=(7.5, 4.2))
x = np.arange(3)
ax.bar(x - 0.27, [calib["target_shares"][k] * 100 for k in keys], 0.27, label="ComCor target", color="#B9C0CC")
ax.bar(x, [first[k] * 100 for k in keys], 0.27, label="model before", color="#6C8EBF")
ax.bar(x + 0.27, [last[k] * 100 for k in keys], 0.27, label="model after fit", color="#D4692A")
ax.set_xticks(x); ax.set_xticklabels(["home", "social events", "work"])
ax.set(ylabel="share of adult infections (%)", title="")
ax.legend(frameon=False)
fig.tight_layout(); fig.savefig(f"{OUT}/05_calibration.png"); plt.close(fig)
print("figures written to", OUT)
