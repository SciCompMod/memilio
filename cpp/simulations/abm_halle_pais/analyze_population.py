#############################################################################
# Copyright (C) 2020-2026 MEmilio
#
# Authors: Sascha Korf
#
# Contact: Martin J. Kuehn <Martin.Kuehn@DLR.de>
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#############################################################################
"""Plausibility check of the synthetic Halle population used by the Halle ABM.

Reads halle_population_data.csv (columns PersonID, home_id, age, Viertel, Home,
school_id, event_id, shopping_id, work_id), prints a summary and writes plots:

  age_distribution.png   persons per year of age, compared with the official
                         population of Halle (county 15002) if available
  household_size.png     distribution of household sizes
  viertel.png            persons per district (Viertel); -99 = outside Halle
  location_sizes.png     persons per school, workplace, shop and event location

Every person carries a school_id and a work_id, whatever their age. The model
only assigns schools to age group 5-14 and workplaces to 15-59, so the location
sizes here are counted for those ages only, as in halle_model.cpp.

Usage (from the repository root):
    python3 cpp/simulations/abm_halle_pais/analyze_population.py \
        --population data/Germany/halle/halle_population_data.csv \
        --county-population data/pydata/Germany/county_current_population.json \
        --out halle_population_plots
"""
import argparse
import json
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HALLE_COUNTY_ID = 15002
OUTSIDE_HALLE = -99

# Age index of the population file -> [lower, upper) age in years. The last
# group is open; 100 is used as its upper bound for the per-year density.
AGE_BINS = [(0, 5), (5, 7), (7, 16), (16, 19), (19, 22), (22, 36), (36, 61),
            (61, 66), (66, 76), (76, 81), (81, 100)]
SCHOOL_AGE_INDICES = [1, 2]        # model age group 5-14
WORK_AGE_INDICES = [3, 4, 5, 6]    # model age groups 15-34 and 35-59

# Columns of county_current_population.json -> [lower, upper) age in years.
OFFICIAL_BINS = {'<3 years': (0, 3), '3-5 years': (3, 6), '6-14 years': (6, 15),
                 '15-17 years': (15, 18), '18-24 years': (18, 25),
                 '25-29 years': (25, 30), '30-39 years': (30, 40),
                 '40-49 years': (40, 50), '50-64 years': (50, 65),
                 '65-74 years': (65, 75), '>74 years': (75, 100)}

# Validated categorical palette, fixed order.
BLUE, ORANGE, AQUA = '#2a78d6', '#eb6834', '#1baf7a'
INK, MUTED, GRID = '#0b0b0b', '#52514e', '#e4e3df'


def style():
    plt.rcParams.update({
        'figure.facecolor': '#fcfcfb', 'axes.facecolor': '#fcfcfb',
        'axes.edgecolor': MUTED, 'axes.labelcolor': INK, 'text.color': INK,
        'xtick.color': MUTED, 'ytick.color': MUTED,
        'axes.spines.top': False, 'axes.spines.right': False,
        'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.8,
        'axes.axisbelow': True, 'font.size': 11, 'axes.titlesize': 13,
        'axes.titleweight': 'bold', 'legend.frameon': False,
        'savefig.dpi': 150, 'savefig.bbox': 'tight'})


def per_year(counts, bins):
    """Turn counts over age bins into persons per year of age, as step data."""
    edges = [b[0] for b in bins] + [bins[-1][1]]
    density = [c / (hi - lo) for c, (lo, hi) in zip(counts, bins)]
    return edges, density + [density[-1]]


def read_official(path):
    if not path or not os.path.exists(path):
        return None
    with open(path) as handle:
        rows = [r for r in json.load(handle) if r.get('ID_County') == HALLE_COUNTY_ID]
    return rows[0] if rows else None


def plot_age(df, official, out):
    fig, ax = plt.subplots(figsize=(9, 5))
    groups = [('Synthetic, Halle residents', df[df.Viertel != OUTSIDE_HALLE], BLUE),
              ('Synthetic, from outside Halle (Viertel -99)', df[df.Viertel == OUTSIDE_HALLE], ORANGE)]
    for label, sub, color in groups:
        counts = sub.age.value_counts().reindex(range(len(AGE_BINS)), fill_value=0)
        edges, dens = per_year(counts.values, AGE_BINS)
        ax.step(edges, dens, where='post', color=color, linewidth=2, label=label)
    if official is not None:
        bins = list(OFFICIAL_BINS.values())
        edges, dens = per_year([official[k] for k in OFFICIAL_BINS], bins)
        ax.step(edges, dens, where='post', color=AQUA, linewidth=2,
                label=f'Official Halle (15002), {official["Population"]:,} total')
    ax.set_xlim(0, 100)
    ax.set_xlabel('Age in years (last group open, drawn up to 100)')
    ax.set_ylabel('Persons per year of age')
    ax.set_title('Age structure')
    ax.legend(loc='upper right')
    fig.savefig(os.path.join(out, 'age_distribution.png'))
    plt.close(fig)


def plot_households(df, out):
    sizes = df[df.Viertel != OUTSIDE_HALLE].groupby('home_id').size()
    share = sizes.value_counts(normalize=True).sort_index() * 100
    fig, ax = plt.subplots(figsize=(7, 4.5))
    bars = ax.bar(share.index, share.values, width=0.8, color=BLUE,
                  edgecolor='#fcfcfb', linewidth=2)
    ax.bar_label(bars, labels=[f'{v:.1f}%' for v in share.values], color=MUTED, padding=3)
    ax.set_xticks(share.index)
    ax.set_xlabel('Persons per household')
    ax.set_ylabel('Share of households (%)')
    ax.set_title(f'Household sizes, Halle residents (mean {sizes.mean():.2f})')
    ax.grid(axis='x', visible=False)
    fig.savefig(os.path.join(out, 'household_size.png'))
    plt.close(fig)


def plot_viertel(df, out):
    counts = df.Viertel.value_counts().sort_values()
    colors = [ORANGE if v == OUTSIDE_HALLE else BLUE for v in counts.index]
    labels = ['outside (-99)' if v == OUTSIDE_HALLE else str(v) for v in counts.index]
    fig, ax = plt.subplots(figsize=(8, 0.22 * len(counts) + 1.5))
    ax.barh(labels, counts.values, height=0.75, color=colors)
    ax.set_xlabel('Persons')
    ax.set_title('Persons per district (Viertel)')
    ax.grid(axis='y', visible=False)
    ax.tick_params(axis='y', labelsize=8)
    fig.savefig(os.path.join(out, 'viertel.png'))
    plt.close(fig)


def location_sizes(df):
    return {
        'School (ages 5-15)': df[df.age.isin(SCHOOL_AGE_INDICES)].groupby('school_id').size(),
        'Work (ages 16-60)': df[df.age.isin(WORK_AGE_INDICES)].groupby('work_id').size(),
        'Basics shop (all)': df.groupby('shopping_id').size(),
        'Social event (all)': df.groupby('event_id').size(),
    }


def plot_locations(sizes, out):
    fig, axes = plt.subplots(2, 2, figsize=(11, 7.5))
    for ax, (name, s) in zip(axes.flat, sizes.items()):
        bins = np.logspace(0, np.log10(s.max() + 1), 30)
        ax.hist(s.values, bins=bins, color=BLUE, edgecolor='#fcfcfb', linewidth=1)
        ax.set_xscale('log')
        ax.set_title(name)
        ax.set_xlabel('Persons assigned (log scale)')
        ax.set_ylabel('Locations')
        ax.text(0.98, 0.95, f'n = {len(s):,}\nmedian = {s.median():,.0f}\nmax = {s.max():,}',
                transform=ax.transAxes, ha='right', va='top', color=MUTED)
        ax.grid(axis='x', visible=False)
    fig.suptitle('Location sizes', fontweight='bold')
    fig.tight_layout()
    fig.savefig(os.path.join(out, 'location_sizes.png'))
    plt.close(fig)


def summary(df, sizes, official):
    resident = df[df.Viertel != OUTSIDE_HALLE]
    hh = resident.groupby('home_id').size()
    print(f'Persons:                 {len(df):,}')
    print(f'  Halle residents:       {len(resident):,}')
    print(f'  outside Halle (-99):   {len(df) - len(resident):,} '
          f'(age indices {sorted(int(a) for a in df[df.Viertel == OUTSIDE_HALLE].age.unique())})')
    if official is not None:
        print(f'  official Halle (15002): {official["Population"]:,} '
              f'-> residents / official = {len(resident) / official["Population"]:.3f}')
    print(f'Households (residents):  {len(hh):,}, mean size {hh.mean():.2f}, '
          f'single {100 * (hh == 1).mean():.1f}%, max {hh.max()}')
    print(f'Districts (Viertel):     {resident.Viertel.nunique()}')
    for name, s in sizes.items():
        print(f'{name:24s} {len(s):6,} locations, median {s.median():7,.0f}, '
              f'p95 {s.quantile(0.95):7,.0f}, max {s.max():7,}')


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--population', default='data/Germany/halle/halle_population_data.csv')
    parser.add_argument('--county-population',
                        default='data/pydata/Germany/county_current_population.json',
                        help='Official population per county; skipped if missing.')
    parser.add_argument('--out', default='halle_population_plots')
    args = parser.parse_args()

    df = pd.read_csv(args.population)
    official = read_official(args.county_population)
    os.makedirs(args.out, exist_ok=True)
    style()

    sizes = location_sizes(df)
    summary(df, sizes, official)
    plot_age(df, official, args.out)
    plot_households(df, args.out)
    plot_viertel(df, args.out)
    plot_locations(sizes, args.out)
    print(f'Plots written to {args.out}/')


if __name__ == '__main__':
    main()
