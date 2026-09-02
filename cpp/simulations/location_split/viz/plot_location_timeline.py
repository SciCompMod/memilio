#!/usr/bin/env python3
"""When and where infections happen, per split level.

Generalizes temporal_infection_heatmap.py / time_loc_type.py from the panvXabmSim branch into an
arm-by-arm heatmap: rows are location types, columns are time, colour is the mean number of new
infections. All arms share one colour scale so they can be read against each other, and --differences
adds one row of heatmaps showing each arm minus the baseline on a diverging scale.

Usage:
    python plot_location_timeline.py --scenario 'no split=out/no_split' \
        --scenario 'half split=out/half_split' --differences
"""

import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import scenario_io as sio


def _grid_per_day(scenario, locations, num_days):
    """Mean new infections per day and location type, shaped (locations, days)."""
    events = sio.load_infection_events(scenario)
    events = events[events["Location_Type"].isin(locations)].copy()
    events["day"] = (events["Timestep"] // 24).astype(int)
    runs = max(1, sio.num_runs(scenario))

    grid = np.zeros((len(locations), num_days))
    index = {name: i for i, name in enumerate(locations)}
    for (location, day), count in events.groupby(["Location_Type", "day"]).size().items():
        if 0 <= day < num_days:
            grid[index[location], day] = count
    return grid / runs


def plot(scenarios, args):
    locations = args.locations
    num_days = int(np.ceil(sio.simulated_days(scenarios[0])))
    grids = {s.label: _grid_per_day(s, locations, num_days) for s in scenarios}

    baseline = sio.baseline_of(scenarios)
    peak = max(grid.max() for grid in grids.values()) or 1.0

    rows = 2 if args.differences else 1
    figure, axes = plt.subplots(rows, len(scenarios), squeeze=False,
                                figsize=(3.6 * len(scenarios) + 1.8, 2.7 * rows + 1.6))

    def draw(axis, grid, cmap, vmin, vmax, title, column, show_xlabel):
        image = axis.imshow(grid, aspect="auto", origin="lower", cmap=cmap, vmin=vmin, vmax=vmax,
                            extent=[0, num_days, -0.5, len(locations) - 0.5])
        axis.set_title(title, fontsize=10)
        axis.set_yticks(range(len(locations)))
        axis.set_yticklabels(locations)
        axis.tick_params(axis="y", labelleft=(column == 0))
        if show_xlabel:
            axis.set_xlabel("day")
        return image

    for column, scenario in enumerate(scenarios):
        title = scenario.label + (" (baseline)" if scenario.is_baseline else "")
        image = draw(axes[0][column], grids[scenario.label], "magma", 0, peak, title, column,
                     show_xlabel=not args.differences)
    figure.colorbar(image, ax=axes[0].tolist(), label="new infections per day",
                    fraction=0.025, pad=0.015)

    if args.differences:
        differences = {s.label: grids[s.label] - grids[baseline.label] for s in scenarios}
        span = max(abs(grid).max() for grid in differences.values()) or 1.0
        for column, scenario in enumerate(scenarios):
            title = "baseline" if scenario.is_baseline else f"{scenario.label} - {baseline.label}"
            image = draw(axes[1][column], differences[scenario.label], "RdBu_r", -span, span,
                         title, column, show_xlabel=True)
        figure.colorbar(image, ax=axes[1].tolist(), label=f"difference to '{baseline.label}'",
                        fraction=0.025, pad=0.015)

    figure.suptitle("New infections by location type over time", fontsize=13)
    return figure


def plot_stacked(scenarios, args):
    """Stacked area of new infections per day, one panel per arm."""
    locations = args.locations
    num_days = int(np.ceil(sio.simulated_days(scenarios[0])))
    days = np.arange(num_days)

    figure, axes = plt.subplots(len(scenarios), 1, squeeze=False, sharex=True, sharey=True,
                                figsize=(10, 2.6 * len(scenarios) + 1.0))
    peak = 0.0
    for row, scenario in enumerate(scenarios):
        axis = axes[row][0]
        grid = _grid_per_day(scenario, locations, num_days)
        peak = max(peak, grid.sum(axis=0).max())
        axis.stackplot(days, grid,
                       colors=[sio.LOCATION_COLORS.get(name, "#777777") for name in locations],
                       labels=locations, edgecolor="white", linewidth=0.3)
        title = scenario.label + (" (baseline)" if scenario.is_baseline else "")
        axis.set_title(title, fontsize=11, loc="left")
        axis.set_ylabel("new infections\nper day")
        axis.grid(alpha=0.2)

    axes[0][0].legend(frameon=False, ncol=len(locations), loc="upper left", fontsize=9)
    axes[-1][0].set_xlabel("day")
    axes[0][0].set_ylim(0, peak * 1.35)
    figure.suptitle("Daily incidence by location type", fontsize=13)
    figure.tight_layout()
    return figure


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    sio.add_scenario_arguments(parser)
    parser.add_argument("--locations", nargs="+", default=sio.TRANSMISSION_LOCATIONS)
    parser.add_argument("--differences", action="store_true",
                        help="Add a row of difference heatmaps against the baseline arm.")
    parser.add_argument("--no-stacked", action="store_true", help="Skip the stacked area figure.")
    args = parser.parse_args()

    scenarios = sio.parse_scenarios(args)
    print("Scenarios:")
    print(sio.describe(scenarios))

    figure = plot(scenarios, args)
    sio.savefig(figure, args, "location_timeline")
    plt.close(figure)

    if not args.no_stacked:
        figure = plot_stacked(scenarios, args)
        sio.savefig(figure, args, "location_timeline_stacked")
        plt.close(figure)


if __name__ == "__main__":
    main()
