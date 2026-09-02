#!/usr/bin/env python3
"""Epidemic curves across split levels.

Replaces infection_cumulated_comp.py / simple_multi_panel.py / multi_seed_comparison.py from the
panvXabmSim branch, which each compared exactly two hardcoded arms. Here the number of arms is
arbitrary and they are ordered along the split axis.

Three panels:
  1. cumulative infections, median with an interquartile band per arm
  2. new infections per day, median per arm
  3. difference to the baseline arm, so that the effect of splitting is read off directly

Usage:
    python plot_epidemic_curves.py \
        --scenario 'no split=out/no_split' --scenario 'half split=out/half_split' \
        --output figures
"""

import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import scenario_io as sio


def _cumulative(scenario, percentile):
    time, total = sio.load_total_infections(scenario, percentile)
    return time, total


def _daily_incidence(time, cumulative):
    """New infections per day from a cumulative curve sampled hourly."""
    steps_per_day = max(1, int(round(1.0 / (time[1] - time[0]))))
    day_edges = np.arange(0, len(cumulative), steps_per_day)
    days = time[day_edges]
    values = cumulative[day_edges]
    return days[1:], np.diff(values)


def plot(scenarios, args):
    figure, axes = plt.subplots(3, 1, figsize=(11, 13), sharex=True)
    ax_cumulative, ax_incidence, ax_difference = axes

    baseline = sio.baseline_of(scenarios)
    baseline_time, baseline_cumulative = _cumulative(baseline, "p50")
    population = sio.population_of(baseline)

    for scenario in scenarios:
        time, median = _cumulative(scenario, "p50")
        _, lower = _cumulative(scenario, "p25")
        _, upper = _cumulative(scenario, "p75")

        style = dict(color=scenario.color, linewidth=2.4 if scenario.is_baseline else 1.9,
                     linestyle="--" if scenario.is_baseline else "-")
        label = scenario.label + (" (baseline)" if scenario.is_baseline else "")

        ax_cumulative.plot(time, median, label=label, **style)
        ax_cumulative.fill_between(time, lower, upper, color=scenario.color, alpha=0.14, linewidth=0)

        if args.show90:
            _, p05 = _cumulative(scenario, "p05")
            _, p95 = _cumulative(scenario, "p95")
            ax_cumulative.fill_between(time, p05, p95, color=scenario.color, alpha=0.07, linewidth=0)

        if args.show_runs:
            for curve in sio.per_run_curves(scenario):
                ax_cumulative.plot(time, curve, color=scenario.color, alpha=0.18, linewidth=0.6)

        days, incidence = _daily_incidence(time, median)
        ax_incidence.plot(days, incidence, **style)

        if len(median) == len(baseline_cumulative):
            ax_difference.plot(time, median - baseline_cumulative, **style)

    ax_cumulative.set_ylabel("cumulative infections")
    ax_cumulative.set_title(
        f"Cumulative infections by split level  (median, IQR band; population {population})")
    ax_cumulative.legend(frameon=False, loc="upper left")
    ax_cumulative.grid(alpha=0.25)

    secondary = ax_cumulative.secondary_yaxis(
        "right", functions=(lambda v: 100 * v / population, lambda v: v * population / 100))
    secondary.set_ylabel("attack rate [%]")

    ax_incidence.set_ylabel("new infections per day")
    ax_incidence.set_title("Daily incidence (median)")
    ax_incidence.grid(alpha=0.25)

    ax_difference.axhline(0, color="#222222", linewidth=1.0, linestyle="--")
    ax_difference.set_ylabel(f"difference to '{baseline.label}'")
    ax_difference.set_xlabel("day")
    ax_difference.set_title("Cumulative infections relative to the baseline arm")
    ax_difference.grid(alpha=0.25)

    figure.tight_layout()
    return figure


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    sio.add_scenario_arguments(parser)
    parser.add_argument("--show90", action="store_true", help="Also shade the 5-95 percentile band.")
    parser.add_argument("--show-runs", action="store_true", help="Draw every single run as a thin line.")
    args = parser.parse_args()

    scenarios = sio.parse_scenarios(args)
    print("Scenarios:")
    print(sio.describe(scenarios))

    figure = plot(scenarios, args)
    sio.savefig(figure, args, "epidemic_curves")
    plt.close(figure)


if __name__ == "__main__":
    main()
