"""Paired analysis of all finished arms: each arm against the no-split arm of its group, run i vs run i."""
import os
STUDY_DIR = os.environ.get("STUDY_DIR", "study_output")
OUT_DIR = os.environ.get("STUDY_OUT", os.path.join(STUDY_DIR, "results"))
BINARY = os.environ.get("SIM_BINARY", "build/bin/location_split_simulation")
import glob, os
import h5py, numpy as np, pandas as pd

BASE = STUDY_DIR
OUT = OUT_DIR
N = 10000
SETTINGS = ['Home', 'School', 'Work', 'SocialEvent', 'BasicsShop']


def baseline_of(name):
    for prefix in ('share10_', 'share15_', 'share25_'):
        if name.startswith(prefix):
            return prefix + 'no_split'
    return 'sw1_no_split'


def curve(arm, run):
    with h5py.File(f'{BASE}/{arm}/total_infections/results_run{run}.h5', 'r') as h:
        return h['0/Time'][()], h['0/Total'][()].ravel()


def day_reaching(t, y, fraction):
    above = np.nonzero(y >= fraction * N)[0]
    return t[above[0]] if len(above) else np.nan


def peak_7d(t, y):
    daily = np.diff(np.interp(np.arange(0, int(t[-1]) + 1), t, y))
    smooth = pd.Series(daily).rolling(7, center=True, min_periods=1).mean().values
    return smooth.max(), np.argmax(smooth) + 1


def setting_shares(arm, runs):
    counts = pd.Series(0.0, index=SETTINGS)
    for r in runs:
        e = pd.read_csv(f'{BASE}/{arm}/runs/run_{r}_infection_events.csv')
        counts = counts.add(e['Location_Type'].value_counts().reindex(SETTINGS).fillna(0), fill_value=0)
    return counts / counts.sum()


def num_runs(arm):
    return len(glob.glob(f'{BASE}/{arm}/total_infections/results_run*.h5'))


rows = []
arms = sorted(os.path.basename(os.path.dirname(p)) for p in glob.glob(f'{BASE}/*/summary.txt'))
for arm in arms:
    base = baseline_of(arm)
    if base not in arms:
        continue
    runs = range(min(num_runs(arm), num_runs(base)))
    final, final_b, d20, d20_b, peak, peak_b = [], [], [], [], [], []
    for r in runs:
        t, y = curve(arm, r)
        tb, yb = curve(base, r)
        final.append(y[-1]); final_b.append(yb[-1])
        d20.append(day_reaching(t, y, 0.2)); d20_b.append(day_reaching(tb, yb, 0.2))
        peak.append(peak_7d(t, y)); peak_b.append(peak_7d(tb, yb))
    final, final_b = np.array(final), np.array(final_b)
    diff = (final - final_b) / N * 100
    se = diff.std(ddof=1) / np.sqrt(len(diff)) if arm != base else 0.0
    delay = np.array(d20) - np.array(d20_b)
    shares = setting_shares(arm, runs)
    rows.append(dict(
        arm=arm, baseline=base, runs=len(runs),
        infected_day90_pct=round(final.mean() / N * 100, 1),
        diff_pp=round(diff.mean(), 2), ci95_low=round(diff.mean() - 1.96 * se, 2), ci95_high=round(diff.mean() + 1.96 * se, 2),
        relative_pct=round((final.mean() / final_b.mean() - 1) * 100, 1),
        runs_lower=int((diff < 0).sum()), runs_higher=int((diff > 0).sum()),
        day_20pct=round(np.nanmean(d20), 1), delay_20pct_days=round(np.nanmean(delay), 2),
        peak_7d_per_day=round(np.mean([p[0] for p in peak]), 1), peak_day=round(np.mean([p[1] for p in peak]), 1),
        peak_change_pct=round((np.mean([p[0] for p in peak]) / np.mean([p[0] for p in peak_b]) - 1) * 100, 1),
        **{f'share_{s}_pct': round(shares[s] * 100, 1) for s in SETTINGS}))

os.makedirs(OUT, exist_ok=True)
pd.DataFrame(rows).to_csv(f'{OUT}/paired_results_day90.csv', index=False)
