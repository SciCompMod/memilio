#!/usr/bin/env python3
"""Infection burden per household (or workplace) over time, per split level.

Generalizes comparative_temporal_heatmap.py from the panvXabmSim branch, which compared exactly two
initialization methods. Rows are arms, columns are time points, and every cell of a panel is one
household coloured by the share of its members that have been infected by that time, averaged over
the runs. Households are sorted by size and separated by thin lines so that the size structure stays
visible, which is what the branch encoded with per-size colours.

Usage:
    python plot_household_heatmap.py --scenario 'no split=out/no_split' \
        --scenario 'half split=out/half_split' --days 5 10 15 21
    python plot_household_heatmap.py ... --use-workplaces
"""

import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import scenario_io as sio


def _group_membership(scenario, use_workplaces):
    """@return (person -> group id, group id -> size), for households or workplaces."""
    if use_workplaces:
        frame = sio.load_workplaces(scenario)
        column = "Work_ID"
    else:
        frame = sio.load_households(scenario)
        column = "Household_ID"
    membership = dict(zip(frame["Person_ID"], frame[column]))
    sizes = frame.groupby(column).size()
    return membership, sizes


def _infection_rates(scenario, membership, sizes, day_cutoffs, use_workplaces):
    """Mean share of infected members per group and cutoff, shaped (cutoffs, groups)."""
    events = sio.load_infection_events(scenario)
    events = events[events["Person_ID"].isin(membership)]
    group_order = list(sizes.index)
    group_index = {group: i for i, group in enumerate(group_order)}
    size_vector = sizes.values.astype(float)

    runs = sorted(events["run"].unique())
    accumulated = np.zeros((len(day_cutoffs), len(group_order)))

    for run in runs:
        run_events = events[events["run"] == run]
        # Seeded infections have Timestep 0 and no location, so they count from the start.
        hours = run_events["Timestep"].values
        groups = np.array([group_index[membership[person]]
                           for person in run_events["Person_ID"].values])
        for cutoff_index, day in enumerate(day_cutoffs):
            mask = hours <= day * 24
            counts = np.bincount(groups[mask], minlength=len(group_order))
            accumulated[cutoff_index] += counts / size_vector

    return np.clip(accumulated / max(1, len(runs)), 0.0, 1.0), group_order


def _layout(sizes, group_order):
    """Grid shape and the row at which each household size block starts."""
    order = sorted(range(len(group_order)), key=lambda i: (sizes.values[i], group_order[i]))
    columns = int(np.ceil(np.sqrt(len(order) * 1.6)))
    rows = int(np.ceil(len(order) / columns))

    boundaries = []
    previous = None
    for position, index in enumerate(order):
        size = sizes.values[index]
        if previous is not None and size != previous:
            boundaries.append(position / columns)
        previous = size
    return order, rows, columns, boundaries


def plot(scenarios, args):
    label_word = "workplace" if args.use_workplaces else "household"
    baseline = sio.baseline_of(scenarios)

    # Padding cells of the last grid row are NaN and must not read as "never infected".
    colormap = matplotlib.colormaps["inferno"].copy()
    colormap.set_bad(alpha=0.0)

    membership, sizes = _group_membership(baseline, args.use_workplaces)
    order, rows, columns, boundaries = _layout(sizes, list(sizes.index))

    figure, axes = plt.subplots(len(scenarios), len(args.days), squeeze=False,
                                figsize=(2.9 * len(args.days) + 1.6, 2.9 * len(scenarios) + 1.2))

    for row, scenario in enumerate(scenarios):
        # Each arm has its own city realization, so the grouping is read per arm.
        arm_membership, arm_sizes = _group_membership(scenario, args.use_workplaces)
        arm_order, arm_rows, arm_columns, arm_boundaries = _layout(arm_sizes, list(arm_sizes.index))
        rates, _ = _infection_rates(scenario, arm_membership, arm_sizes, args.days,
                                    args.use_workplaces)

        for column, day in enumerate(args.days):
            axis = axes[row][column]
            grid = np.full(arm_rows * arm_columns, np.nan)
            grid[:len(arm_order)] = rates[column][arm_order]
            image = axis.imshow(grid.reshape(arm_rows, arm_columns), cmap=colormap,
                                vmin=0, vmax=1, aspect="auto", interpolation="nearest")
            for boundary in arm_boundaries:
                axis.axhline(boundary - 0.5, color="#5FD3F3", linewidth=0.8, alpha=0.9)

            axis.set_xticks([])
            axis.set_yticks([])
            if row == 0:
                axis.set_title(f"day {day:g}", fontsize=11)
            if column == 0:
                name = scenario.label + ("\n(baseline)" if scenario.is_baseline else "")
                axis.set_ylabel(name, fontsize=10)
            mean_rate = np.nanmean(rates[column])
            axis.set_xlabel(f"mean {mean_rate * 100:.0f}%", fontsize=9)

    figure.colorbar(image, ax=axes.ravel().tolist(),
                    label=f"share of {label_word} members ever infected",
                    fraction=0.022, pad=0.015)
    figure.suptitle(
        f"Infection burden per {label_word} over time  "
        f"(mean over runs; {label_word}s sorted by size, blue lines separate the sizes)",
        fontsize=12)
    return figure


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    sio.add_scenario_arguments(parser)
    parser.add_argument("--days", nargs="+", type=float, default=None,
                        help="Days to show. Default: four evenly spaced points.")
    parser.add_argument("--use-workplaces", action="store_true",
                        help="Group by workplace instead of household.")
    args = parser.parse_args()

    scenarios = sio.parse_scenarios(args)
    print("Scenarios:")
    print(sio.describe(scenarios))

    if args.days is None:
        total = sio.simulated_days(scenarios[0])
        args.days = [round(total * fraction, 1) for fraction in (0.25, 0.5, 0.75, 1.0)]

    figure = plot(scenarios, args)
    name = "workplace_heatmap" if args.use_workplaces else "household_heatmap"
    sio.savefig(figure, args, name)
    plt.close(figure)


if __name__ == "__main__":
    main()
