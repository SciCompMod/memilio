#!/usr/bin/env python3
"""Transmission trees and generation structure, per split level.

Reconstructs who infected whom and compares the resulting chains across arms. This replaces
infection_timeline.py / comparative_infection_trees.py from the panvXabmSim branch, which laid out
two hardcoded arms side by side with a heuristic infector search over contact records.

The reconstruction here uses the two things the simulation writes: every infection with the exact
Location it happened at, and where every Person was at every hour. The infector of a Person infected
during hour t at Location L must have been at L during hour t and already infectious, which usually
narrows the candidates down to a handful; among those the one whose own infection time puts it
furthest into its infectious window is chosen.

Needs the simulation to have been run with --write_person_locations.

Usage:
    python plot_transmission_tree.py --scenario 'no split=out/no_split' \
        --scenario 'half split=out/half_split' --run 0
"""

import argparse
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import scenario_io as sio

# Hours after its own infection during which a Person is treated as a plausible infector. Roughly
# the exposed period through the end of the symptomatic period for the parameters of this model.
INFECTIOUS_WINDOW = (24, 14 * 24)


def build_tree(scenario, run):
    """@return (nodes, parent) where nodes maps person -> dict(time, location_type, generation)."""
    events = sio.load_infection_events(scenario)
    events = events[events["run"] == run]
    if events.empty:
        raise ValueError(f"No infection events for run {run} of '{scenario.label}'.")

    infection_hour = dict(zip(events["Person_ID"], events["Timestep"]))
    seeds = set(events.loc[events["Location_Type"] == "Initial", "Person_ID"])

    locations = sio.load_person_locations(scenario, run)
    locations = locations[locations["Person_ID"].isin(infection_hour)]
    # (hour, location) -> the infected persons present. Only ever-infected persons can be infectors.
    present = defaultdict(list)
    for hour, location, person in zip(locations["Timestep"], locations["Location_ID"],
                                      locations["Person_ID"]):
        present[(hour, location)].append(person)

    nodes = {}
    parent = {}
    for person in seeds:
        nodes[person] = {"time": 0, "location_type": "Initial", "generation": 0}

    transmissions = events[events["Location_Type"] != "Initial"].sort_values("Timestep")
    unresolved = 0

    for person, hour, location, location_type in zip(
            transmissions["Person_ID"], transmissions["Timestep"],
            transmissions["Location_ID"], transmissions["Location_Type"]):
        # The infection happened during the step before it was observed, at the Location the
        # Person occupied then; see derive_infection_events in multi_run_simulator.cpp.
        candidates = []
        for other in present.get((hour - 1, int(location)), ()):
            if other == person:
                continue
            other_hour = infection_hour.get(other)
            if other_hour is None or other_hour >= hour:
                continue
            age = hour - other_hour
            score = 1.0 if INFECTIOUS_WINDOW[0] <= age <= INFECTIOUS_WINDOW[1] else 0.0
            candidates.append((score, -abs(age - 5 * 24), other))

        generation = None
        if candidates:
            infector = max(candidates)[2]
            parent[person] = infector
            if infector in nodes:
                generation = nodes[infector]["generation"] + 1
        else:
            unresolved += 1

        nodes[person] = {"time": hour, "location_type": location_type,
                         "generation": generation if generation is not None else -1}

    # Persons whose infector was resolved only later get their generation filled in afterwards.
    for person in sorted(nodes, key=lambda p: nodes[p]["time"]):
        if nodes[person]["generation"] == -1 and person in parent:
            infector = parent[person]
            if infector in nodes and nodes[infector]["generation"] >= 0:
                nodes[person]["generation"] = nodes[infector]["generation"] + 1

    resolved = sum(1 for node in nodes.values() if node["generation"] >= 0)
    print(f"  {scenario.label}: {len(nodes)} infections, {resolved} placed in the tree, "
          f"{unresolved} without an identifiable infector")
    return nodes, parent


def plot_trees(trees, scenarios, args):
    figure, axes = plt.subplots(len(scenarios), 1, squeeze=False, sharex=True,
                                figsize=(11, 3.1 * len(scenarios) + 1.0))

    for row, scenario in enumerate(scenarios):
        axis = axes[row][0]
        nodes, parent = trees[scenario.label]
        placed = {p: n for p, n in nodes.items() if n["generation"] >= 0}

        for person, node in placed.items():
            infector = parent.get(person)
            if infector in placed:
                axis.plot([placed[infector]["time"] / 24.0, node["time"] / 24.0],
                          [placed[infector]["generation"], node["generation"]],
                          color="#999999", linewidth=0.35, alpha=0.25, zorder=1)

        for location in ["Initial"] + args.locations:
            subset = [n for n in placed.values() if n["location_type"] == location]
            if not subset:
                continue
            axis.scatter([n["time"] / 24.0 for n in subset], [n["generation"] for n in subset],
                         s=9, alpha=0.65, zorder=2, label=location,
                         color=sio.LOCATION_COLORS.get(location, "#777777"))

        title = scenario.label + (" (baseline)" if scenario.is_baseline else "")
        axis.set_title(f"{title} - run {args.run}", fontsize=11, loc="left")
        axis.set_ylabel("generation")
        axis.grid(alpha=0.2)

    axes[0][0].legend(frameon=False, ncol=6, fontsize=9, loc="upper left")
    axes[-1][0].set_xlabel("day of infection")
    figure.suptitle("Transmission chains: infection generation over time, coloured by where it happened",
                    fontsize=12)
    figure.tight_layout()
    return figure


def plot_generations(trees, scenarios, args):
    """Generation distribution per arm. This is the part that compares directly across arms."""
    figure, axes = plt.subplots(1, 2, figsize=(13, 4.6))
    ax_histogram, ax_mean = axes

    def generations_of(scenario):
        nodes, _ = trees[scenario.label]
        return [n["generation"] for n in nodes.values() if n["generation"] >= 0]

    max_generation = 0
    for scenario in scenarios:
        generations = generations_of(scenario)
        max_generation = max(max_generation, max(generations) if generations else 0)

    bins = np.arange(-0.5, max_generation + 1.5)
    width = 0.8 / len(scenarios)
    for index, scenario in enumerate(scenarios):
        generations = generations_of(scenario)
        counts, _ = np.histogram(generations, bins=bins)
        offset = (index - (len(scenarios) - 1) / 2) * width
        ax_histogram.bar(np.arange(max_generation + 1) + offset, counts, width=width * 0.92,
                         color=scenario.color, edgecolor="white", linewidth=0.4,
                         label=scenario.label + (" (baseline)" if scenario.is_baseline else ""))

    ax_histogram.set_xlabel("generation")
    ax_histogram.set_ylabel("infections")
    ax_histogram.set_title("How deep the transmission chains go")
    ax_histogram.legend(frameon=False)
    ax_histogram.grid(axis="y", alpha=0.25)

    # Mean offspring per generation: the effective reproduction number as the tree sees it.
    for scenario in scenarios:
        nodes, parent = trees[scenario.label]
        offspring = defaultdict(int)
        for person, infector in parent.items():
            if infector in nodes and nodes[infector]["generation"] >= 0:
                offspring[infector] += 1
        by_generation = defaultdict(list)
        for person, node in nodes.items():
            if node["generation"] >= 0:
                by_generation[node["generation"]].append(offspring.get(person, 0))
        generations = sorted(by_generation)
        means = [np.mean(by_generation[g]) for g in generations]
        ax_mean.plot(generations, means, marker="o", markersize=4,
                     color=scenario.color, linewidth=2.0 if scenario.is_baseline else 1.6,
                     linestyle="--" if scenario.is_baseline else "-",
                     label=scenario.label)

    ax_mean.axhline(1.0, color="#888888", linewidth=0.9, linestyle=":")
    ax_mean.set_xlabel("generation")
    ax_mean.set_ylabel("mean secondary infections")
    ax_mean.set_title("Effective reproduction by generation")
    ax_mean.legend(frameon=False)
    ax_mean.grid(alpha=0.25)

    figure.tight_layout()
    return figure


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    sio.add_scenario_arguments(parser)
    parser.add_argument("--run", type=int, default=0, help="Which run to reconstruct.")
    parser.add_argument("--locations", nargs="+", default=sio.TRANSMISSION_LOCATIONS)
    args = parser.parse_args()

    scenarios = sio.parse_scenarios(args)
    print("Scenarios:")
    print(sio.describe(scenarios))

    print("Reconstructing transmission trees...")
    trees = {s.label: build_tree(s, args.run) for s in scenarios}

    figure = plot_trees(trees, scenarios, args)
    sio.savefig(figure, args, "transmission_trees")
    plt.close(figure)

    figure = plot_generations(trees, scenarios, args)
    sio.savefig(figure, args, "transmission_generations")
    plt.close(figure)


if __name__ == "__main__":
    main()
