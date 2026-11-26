"""Find school weight w and k so that schools cause a target share of all infections over the epidemic,
with the early growth rate (weeks 2-3) of the calibrated model. 10 runs, 90 days per evaluation."""
import os
STUDY_DIR = os.environ.get("STUDY_DIR", "study_output")
OUT_DIR = os.environ.get("STUDY_OUT", os.path.join(STUDY_DIR, "results"))
BINARY = os.environ.get("SIM_BINARY", "build/bin/location_split_simulation")
import glob, json, math, os, subprocess, sys
import numpy as np, pandas as pd
B = BINARY; D = os.path.join(STUDY_DIR, 'scan')
os.makedirs(D, exist_ok=True)

def evaluate(w, k):
    out = f'{D}/w{w:.4f}_k{k:.4f}'
    if not os.path.exists(out + '/summary.txt'):
        subprocess.run([B, '--runs', '10', '--days', '90', '--n_persons', '10000', '--initial_infections', '10', '--seed', '42',
                        *(['--threads', os.environ['THREADS']] if os.environ.get('THREADS') else []), '--infection_k', f'{k:.4f}',
                        '--contact_weights', f'home=21.616,work=0.625,social_event=4.8,shop=0.1584,school={1.6 * w:.4f}', '--output_dir', out],
                       check=True, stdout=subprocess.DEVNULL)
    e = pd.concat(pd.read_csv(f) for f in glob.glob(out + '/runs/*events.csv'))
    e = e[e.Location_Type != 'Initial']
    share = (e.Location_Type == 'School').mean()
    day = (e.Timestep - 1) // 24
    growth = math.log(((day >= 14) & (day < 21)).sum() / ((day >= 7) & (day < 14)).sum()) / 7
    return share, growth

log = open(f'{D}/scan_log.txt', 'a')
s0, g0 = evaluate(1.0, 2.11)
print(f'reference w=1 k=2.11: share {s0:.3f} growth {g0:.4f}', file=log, flush=True)
result = {}
for target in [0.15, 0.25]:
    odds = lambda s: s / (1 - s)
    w, k = odds(target) / odds(s0), 2.11
    k_slope, prev = 0.1, None
    for it in range(5):
        s, g = evaluate(w, k)
        print(f'target {target}: it {it} w={w:.3f} k={k:.3f} share {s:.3f} growth {g:.4f}', file=log, flush=True)
        if abs(s - target) < 0.01 and abs(g - g0) < 0.005:
            break
        if prev is not None and abs(math.log(k / prev[1])) > 1e-6:
            k_slope = float(np.clip((g - prev[3]) / math.log(k / prev[1]), 0.03, 0.5))
        prev = (w, k, s, g)
        w *= odds(target) / odds(s)
        k *= math.exp(np.clip((g0 - g) / k_slope, -0.3, 0.3))
    result[str(target)] = {'school_weight': 1.6 * w, 'school_weight_relative': w, 'k': k, 'share': s, 'growth': g}
json.dump({'reference': {'share': s0, 'growth': g0}, 'targets': result}, open(f'{D}/scan_result.json', 'w'), indent=2)
print('done', file=log, flush=True)
