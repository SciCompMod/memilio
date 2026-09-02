#!/usr/bin/env python3
"""Run the whole comparison: every figure plus a summary table.

Replaces compare_all_scenarios.py / combined_scenario_comparison.py from the panvXabmSim branch.
Those drove a fixed set of event types; this one drives an arbitrary number of split levels.

Usage:
    python compare_scenarios.py \
        --scenario 'no split=out/no_split' \
        --scenario 'quarter split=out/quarter_split' \
        --scenario 'half split=out/half_split' \
        --scenario 'full split=out/full_split' \
        --output figures
"""

import argparse
import copy
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import plot_epidemic_curves
import plot_household_heatmap
import plot_location_breakdown
import plot_location_timeline
import plot_transmission_tree
import scenario_io as sio


def summary_table(scenarios, args) -> pd.DataFrame:
    """Per arm: outbreak size, timing, and where transmission happened."""
    baseline = sio.baseline_of(scenarios)
    _, baseline_final = sio.load_total_infections(baseline, "p50")

    rows = []
    for scenario in scenarios:
        time, median = sio.load_total_infections(scenario, "p50")
        curves = sio.per_run_curves(scenario)
        population = sio.population_of(scenario)

        steps_per_day = max(1, int(round(1.0 / (time[1] - time[0]))))
        day_index = np.arange(0, len(median), steps_per_day)
        daily = np.diff(median[day_index])

        counts = sio.events_per_run_by(scenario, "Location_Type", sio.TRANSMISSION_LOCATIONS)
        location_means = counts.mean()
        location_total = location_means.sum() or 1.0

        row = {
            "scenario": scenario.label,
            "split_fraction": scenario.split_fraction,
            "baseline": scenario.is_baseline,
            "runs": sio.num_runs(scenario),
            "population": population,
            "final_size_median": float(median[-1]),
            "final_size_p25": float(np.percentile(curves[:, -1], 25)),
            "final_size_p75": float(np.percentile(curves[:, -1], 75)),
            "attack_rate_pct": 100.0 * float(median[-1]) / population,
            "vs_baseline_abs": float(median[-1]) - float(baseline_final[-1]),
            "vs_baseline_pct": (100.0 * (float(median[-1]) - float(baseline_final[-1]))
                                / max(1.0, float(baseline_final[-1]))),
            "peak_daily_incidence": float(daily.max()) if len(daily) else float("nan"),
            "peak_day": float(np.argmax(daily) + 1) if len(daily) else float("nan"),
        }
        for location in sio.TRANSMISSION_LOCATIONS:
            row[f"share_{location}_pct"] = 100.0 * float(location_means[location]) / location_total
        rows.append(row)

    return pd.DataFrame(rows)


def _namespace(args, **extra):
    namespace = copy.copy(args)
    for key, value in extra.items():
        setattr(namespace, key, value)
    return namespace


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    sio.add_scenario_arguments(parser)
    parser.add_argument("--run", type=int, default=0,
                        help="Run used for the transmission tree.")
    parser.add_argument("--days", nargs="+", type=float, default=None,
                        help="Time points for the household heatmap.")
    parser.add_argument("--skip-trees", action="store_true",
                        help="Skip the transmission trees; they need --write_person_locations.")
    args = parser.parse_args()

    scenarios = sio.parse_scenarios(args)
    print("Scenarios:")
    print(sio.describe(scenarios))
    os.makedirs(args.output, exist_ok=True)

    print("\nSummary statistics")
    table = summary_table(scenarios, args)
    path = os.path.join(args.output, "summary_statistics.csv")
    table.to_csv(path, index=False)
    print(f"  wrote {path}\n")
    columns = ["scenario", "final_size_median", "attack_rate_pct", "vs_baseline_abs",
               "vs_baseline_pct", "peak_daily_incidence", "peak_day"]
    print(table[columns].to_string(index=False, float_format=lambda v: f"{v:.1f}"))

    print("\nEpidemic curves")
    figure = plot_epidemic_curves.plot(scenarios, _namespace(args, show90=False, show_runs=False))
    sio.savefig(figure, args, "epidemic_curves")
    plt.close(figure)

    print("Location breakdown")
    breakdown_args = _namespace(args, locations=sio.TRANSMISSION_LOCATIONS, pies=False)
    figure = plot_location_breakdown.plot(scenarios, breakdown_args)
    sio.savefig(figure, args, "location_breakdown")
    plt.close(figure)

    print("Location timeline")
    timeline_args = _namespace(args, locations=sio.TRANSMISSION_LOCATIONS, differences=True)
    figure = plot_location_timeline.plot(scenarios, timeline_args)
    sio.savefig(figure, args, "location_timeline")
    plt.close(figure)
    figure = plot_location_timeline.plot_stacked(scenarios, timeline_args)
    sio.savefig(figure, args, "location_timeline_stacked")
    plt.close(figure)

    print("Household burden")
    days = args.days
    if days is None:
        total = sio.simulated_days(scenarios[0])
        days = [round(total * fraction, 1) for fraction in (0.25, 0.5, 0.75, 1.0)]
    figure = plot_household_heatmap.plot(scenarios, _namespace(args, days=days,
                                                               use_workplaces=False))
    sio.savefig(figure, args, "household_heatmap")
    plt.close(figure)

    if not args.skip_trees:
        print("Transmission trees")
        try:
            tree_args = _namespace(args, locations=sio.TRANSMISSION_LOCATIONS)
            trees = {s.label: plot_transmission_tree.build_tree(s, args.run) for s in scenarios}
            figure = plot_transmission_tree.plot_trees(trees, scenarios, tree_args)
            sio.savefig(figure, args, "transmission_trees")
            plt.close(figure)
            figure = plot_transmission_tree.plot_generations(trees, scenarios, tree_args)
            sio.savefig(figure, args, "transmission_generations")
            plt.close(figure)
        except FileNotFoundError as error:
            print(f"  skipped: {error}")

    print(f"\nDone. Figures in {args.output}")


if __name__ == "__main__":
    main()
