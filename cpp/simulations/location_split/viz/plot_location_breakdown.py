#!/usr/bin/env python3
"""Where infections happen, per split level.

Generalizes location_infection_pie_chart.py from the panvXabmSim branch. A pie chart per arm is
hard to compare once there are more than two arms, so the default is a share plot plus absolute
counts with the spread over runs; --pies reproduces the original pie layout as well.

Panels:
  1. share of infections per location type, one stacked bar per arm
  2. absolute infections per location type, arms grouped, whiskers = min/max over runs
  3. difference to the baseline arm per location type

Usage:
    python plot_location_breakdown.py \
        --scenario 'no split=out/no_split' --scenario 'half split=out/half_split'
"""

import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import scenario_io as sio


def _counts(scenario, locations):
    """@return (mean, min, max) infections per location type over the runs."""
    per_run = sio.events_per_run_by(scenario, "Location_Type", locations)
    return per_run.mean().values, per_run.min().values, per_run.max().values


def plot(scenarios, args):
    locations = args.locations
    baseline = sio.baseline_of(scenarios)

    statistics = {s.label: _counts(s, locations) for s in scenarios}
    baseline_mean = statistics[baseline.label][0]

    figure, axes = plt.subplots(3, 1, figsize=(11, 13))
    ax_share, ax_absolute, ax_difference = axes

    labels = [s.label for s in scenarios]
    positions = np.arange(len(scenarios))

    # 1. composition
    bottom = np.zeros(len(scenarios))
    for index, location in enumerate(locations):
        means = np.array([statistics[label][0][index] for label in labels])
        totals = np.array([statistics[label][0].sum() for label in labels])
        share = 100 * means / np.where(totals == 0, 1, totals)
        ax_share.bar(positions, share, bottom=bottom, width=0.62,
                     color=sio.LOCATION_COLORS.get(location, "#777777"),
                     edgecolor="white", linewidth=0.7, label=location)
        for x, (value, base) in enumerate(zip(share, bottom)):
            if value >= 4:
                ax_share.text(x, base + value / 2, f"{value:.0f}%", ha="center", va="center",
                              fontsize=9, color="white", fontweight="bold")
        bottom += share

    ax_share.set_xticks(positions)
    ax_share.set_xticklabels(labels)
    ax_share.set_ylabel("share of infections [%]")
    ax_share.set_ylim(0, 100)
    ax_share.set_title("Composition of transmission by location type")
    ax_share.legend(frameon=False, ncol=len(locations), loc="upper center",
                    bbox_to_anchor=(0.5, -0.08))

    # 2. absolute counts
    width = 0.8 / len(scenarios)
    group = np.arange(len(locations))
    for index, scenario in enumerate(scenarios):
        mean, low, high = statistics[scenario.label]
        offset = (index - (len(scenarios) - 1) / 2) * width
        ax_absolute.bar(group + offset, mean, width=width * 0.92, color=scenario.color,
                        label=scenario.label + (" (baseline)" if scenario.is_baseline else ""),
                        edgecolor="white", linewidth=0.5)
        ax_absolute.errorbar(group + offset, mean, yerr=[mean - low, high - mean],
                             fmt="none", ecolor="#333333", elinewidth=0.9, capsize=2.5)

    ax_absolute.set_xticks(group)
    ax_absolute.set_xticklabels(locations)
    ax_absolute.set_ylabel("infections per run")
    ax_absolute.set_title("Infections per location type (mean over runs, whiskers = min/max)")
    ax_absolute.legend(frameon=False)
    ax_absolute.grid(axis="y", alpha=0.25)

    # 3. difference to baseline
    for index, scenario in enumerate(scenarios):
        if scenario.is_baseline:
            continue
        mean = statistics[scenario.label][0]
        offset = (index - (len(scenarios) - 1) / 2) * width
        ax_difference.bar(group + offset, mean - baseline_mean, width=width * 0.92,
                          color=scenario.color, label=scenario.label,
                          edgecolor="white", linewidth=0.5)

    ax_difference.axhline(0, color="#222222", linewidth=1.0)
    ax_difference.set_xticks(group)
    ax_difference.set_xticklabels(locations)
    ax_difference.set_ylabel(f"difference to '{baseline.label}'")
    ax_difference.set_title("Change in infections per location type relative to the baseline arm")
    ax_difference.legend(frameon=False)
    ax_difference.grid(axis="y", alpha=0.25)

    figure.tight_layout()
    return figure


def plot_pies(scenarios, args):
    """The original per-scenario pie chart layout, one pie per arm."""
    locations = args.locations
    columns = min(len(scenarios), 4)
    rows = int(np.ceil(len(scenarios) / columns))
    figure, axes = plt.subplots(rows, columns, figsize=(4.2 * columns, 4.6 * rows), squeeze=False)

    for index, scenario in enumerate(scenarios):
        axis = axes[index // columns][index % columns]
        mean = sio.events_per_run_by(scenario, "Location_Type", locations).mean()
        axis.pie(mean.values, labels=locations, autopct="%1.1f%%", startangle=90,
                 colors=[sio.LOCATION_COLORS.get(name, "#777777") for name in locations],
                 textprops={"fontsize": 9},
                 wedgeprops={"edgecolor": "white", "linewidth": 1.0})
        title = scenario.label + (" (baseline)" if scenario.is_baseline else "")
        axis.set_title(f"{title}\n{mean.sum():.0f} infections per run", fontsize=11)

    for index in range(len(scenarios), rows * columns):
        axes[index // columns][index % columns].axis("off")

    figure.suptitle("Infections by location type", fontsize=14)
    figure.tight_layout()
    return figure


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    sio.add_scenario_arguments(parser)
    parser.add_argument("--locations", nargs="+", default=sio.TRANSMISSION_LOCATIONS,
                        help="Location types to show.")
    parser.add_argument("--pies", action="store_true", help="Also write the pie chart layout.")
    args = parser.parse_args()

    scenarios = sio.parse_scenarios(args)
    print("Scenarios:")
    print(sio.describe(scenarios))

    figure = plot(scenarios, args)
    sio.savefig(figure, args, "location_breakdown")
    plt.close(figure)

    if args.pies:
        figure = plot_pies(scenarios, args)
        sio.savefig(figure, args, "location_pies")
        plt.close(figure)


if __name__ == "__main__":
    main()
