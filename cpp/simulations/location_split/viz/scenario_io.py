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
"""Shared loading and styling for the location split comparison plots.

A *scenario* is one output directory of the location_split_simulation, i.e. one arm of the
comparison. Arms are ordinal: they differ in how much of the split location is handled outside the
ABM, from no split (the baseline) through quarter and half up to a full split. Every plot in this
directory therefore takes an arbitrary number of arms, orders them by that split fraction and draws
the baseline as the reference the others are compared against.

How the split fraction of an arm is determined, in order of precedence:

1. ``scenario.json`` in the results directory::

       {"label": "half split", "split_fraction": 0.5, "split_location_type": "SocialEvent"}

2. the label given on the command line or the directory name, parsed for a known word
   (``no``/``quarter``/``third``/``half``/``three_quarter``/``full``) or a number
   (``split_0.75``, ``0.75``).

3. the position on the command line, if neither of the above says anything.

The split mechanism itself is not implemented in the simulation yet. Point 1 is the seam it will
plug into: once an arm is actually produced by splitting, it writes its own scenario.json and
nothing in these scripts has to change.
"""

import argparse
import glob
import json
import os
import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional

import h5py
import numpy as np
import pandas as pd

# Order matters: it is the column order in the h5 files written by the simulation.
INFECTION_STATES = [
    "Susceptible",
    "Exposed",
    "InfectedNoSymptoms",
    "InfectedSymptoms",
    "InfectedSevere",
    "InfectedCritical",
    "Recovered",
    "Dead",
]

LOCATION_TYPES = [
    "Home",
    "School",
    "Work",
    "SocialEvent",
    "BasicsShop",
    "Hospital",
    "ICU",
    "Car",
    "PublicTransport",
    "TransportWithoutContact",
    "Cemetery",
]

# The location types that actually carry transmission in this model. Used to keep the plots free of
# empty categories.
TRANSMISSION_LOCATIONS = ["Home", "Work", "School", "SocialEvent", "BasicsShop"]

LOCATION_COLORS = {
    "Home": "#4C72B0",
    "Work": "#DD8452",
    "School": "#55A868",
    "SocialEvent": "#C44E52",
    "BasicsShop": "#8172B3",
    "Hospital": "#937860",
    "ICU": "#DA8BC3",
    "Initial": "#999999",
}

AGE_GROUPS = ["0-4", "5-14", "15-34", "35-59", "60-79", "80+"]

# Words that may appear in a scenario label, and the share of the split location they stand for.
_SPLIT_WORDS = {
    "no": 0.0,
    "none": 0.0,
    "baseline": 0.0,
    "quarter": 0.25,
    "third": 1.0 / 3.0,
    "half": 0.5,
    "three_quarter": 0.75,
    "threequarter": 0.75,
    "full": 1.0,
    "complete": 1.0,
}


@dataclass
class Scenario:
    """One arm of the comparison."""

    label: str
    path: str
    split_fraction: Optional[float] = None
    split_location_type: Optional[str] = None
    color: str = "black"
    is_baseline: bool = False
    _events: Optional[pd.DataFrame] = field(default=None, repr=False)

    def file(self, *parts: str) -> str:
        return os.path.join(self.path, *parts)


def split_fraction_from_label(label: str) -> Optional[float]:
    """Guess the split fraction from a label such as ``half_split`` or ``split_0.75``."""
    normalized = label.lower().replace("-", "_").replace(" ", "_")

    match = re.search(r"(\d+(?:\.\d+)?)\s*%", normalized)
    if match:
        return float(match.group(1)) / 100.0
    match = re.search(r"(?:^|_)(\d*\.\d+|[01])(?:_|$)", normalized)
    if match:
        return float(match.group(1))

    for word, fraction in _SPLIT_WORDS.items():
        if re.search(rf"(?:^|_){word}(?:_|$)", normalized):
            return fraction
    return None


def _read_scenario_json(path: str) -> dict:
    config_file = os.path.join(path, "scenario.json")
    if not os.path.exists(config_file):
        return {}
    with open(config_file) as handle:
        return json.load(handle)


def add_scenario_arguments(parser: argparse.ArgumentParser) -> None:
    """Add the arguments every comparison plot understands."""
    parser.add_argument(
        "--scenario",
        action="append",
        required=True,
        metavar="[LABEL=]PATH",
        help="One arm of the comparison. Repeat for every arm. Without LABEL the directory name "
        "is used. Example: --scenario 'no split=out/no_split' --scenario 'half split=out/half'",
    )
    parser.add_argument("--baseline", default=None,
                        help="Label of the arm the others are compared against. "
                             "Default: the arm with the smallest split fraction.")
    parser.add_argument("-o", "--output", default="viz_out", help="Directory for the figures.")
    parser.add_argument("--dpi", type=int, default=150)
    parser.add_argument("--format", default="png", choices=["png", "pdf", "svg"])


def parse_scenarios(args: argparse.Namespace) -> List[Scenario]:
    """Turn the --scenario arguments into an ordered list of arms."""
    scenarios: List[Scenario] = []

    for position, spec in enumerate(args.scenario):
        if "=" in spec:
            label, path = spec.split("=", 1)
        else:
            label, path = os.path.basename(os.path.normpath(spec)), spec
        label, path = label.strip(), os.path.expanduser(path.strip())

        if not os.path.isdir(path):
            raise FileNotFoundError(f"Scenario directory does not exist: {path}")

        config = _read_scenario_json(path)
        label = config.get("label", label)
        fraction = config.get("split_fraction", split_fraction_from_label(label))

        scenarios.append(Scenario(label=label, path=path, split_fraction=fraction,
                                  split_location_type=config.get("split_location_type")))
        scenarios[-1]._position = position  # type: ignore[attr-defined]

    known = [s.split_fraction for s in scenarios if s.split_fraction is not None]
    if len(known) != len(scenarios):
        missing = [s.label for s in scenarios if s.split_fraction is None]
        print(f"Note: no split fraction for {missing}; keeping the command line order.")
        for position, scenario in enumerate(scenarios):
            scenario.split_fraction = float(position)

    scenarios.sort(key=lambda s: (s.split_fraction, s.label))
    _assign_roles(scenarios, getattr(args, "baseline", None))
    return scenarios


def _assign_roles(scenarios: List[Scenario], baseline_label: Optional[str]) -> None:
    """Mark the baseline arm and give every arm a colour along the split axis."""
    import matplotlib

    baseline_index = 0
    if baseline_label is not None:
        matches = [i for i, s in enumerate(scenarios) if s.label == baseline_label]
        if not matches:
            raise ValueError(f"Unknown baseline '{baseline_label}'. "
                             f"Known arms: {[s.label for s in scenarios]}")
        baseline_index = matches[0]

    scenarios[baseline_index].is_baseline = True

    # The baseline is the reference and gets a neutral colour, the split arms a sequential ramp so
    # that "more split" reads as "further along the colour axis" even with six of them.
    others = [i for i in range(len(scenarios)) if i != baseline_index]
    ramp = matplotlib.colormaps["viridis"]
    for rank, index in enumerate(others):
        position = 0.12 + 0.72 * (rank / max(1, len(others) - 1)) if len(others) > 1 else 0.45
        scenarios[index].color = ramp(position)
    scenarios[baseline_index].color = "#222222"


def baseline_of(scenarios: List[Scenario]) -> Scenario:
    for scenario in scenarios:
        if scenario.is_baseline:
            return scenario
    return scenarios[0]


def describe(scenarios: List[Scenario]) -> str:
    rows = []
    for scenario in scenarios:
        marker = " (baseline)" if scenario.is_baseline else ""
        rows.append(f"  {scenario.label}{marker}: split={scenario.split_fraction:g}  {scenario.path}")
    return "\n".join(rows)


# --------------------------------------------------------------------------------------------
# Loading
# --------------------------------------------------------------------------------------------

def load_h5(scenario: Scenario, quantity: str, percentile: str = "p50") -> Dict[str, np.ndarray]:
    """Load one percentile file of one quantity.

    @param quantity One of total_infections, infection_state_per_age_group,
                    infection_per_location_type_per_age_group.
    @param percentile p05, p25, p50, p75, p95, or run<i> for a single run.
    """
    if percentile.startswith("run"):
        path = scenario.file(quantity, f"results_{percentile}.h5")
    else:
        path = scenario.file(quantity, percentile, "Results.h5")
    if not os.path.exists(path):
        raise FileNotFoundError(f"Missing results file: {path}")
    with h5py.File(path, "r") as handle:
        return {key: np.array(value) for key, value in handle["0"].items()}


def num_runs(scenario: Scenario) -> int:
    return len(glob.glob(scenario.file("total_infections", "results_run*.h5")))


def load_total_infections(scenario: Scenario, percentile: str = "p50"):
    """@return (time in days, cumulative number of infected persons)."""
    data = load_h5(scenario, "total_infections", percentile)
    return data["Time"], data["Total"][:, 0]


def load_infection_events(scenario: Scenario) -> pd.DataFrame:
    """All infection events of all runs, with a ``run`` column added.

    Columns: run, Timestep, Person_ID, Location_ID, Location_Type, Age_Group.
    The rows with Location_Type == "Initial" are the seeded infections; they have no location.
    """
    if scenario._events is not None:
        return scenario._events

    files = sorted(glob.glob(scenario.file("runs", "run_*_infection_events.csv")))
    if not files:
        raise FileNotFoundError(
            f"No infection events in {scenario.file('runs')}. "
            "Run the simulation without --no_infection_events.")

    frames = []
    for path in files:
        run = int(re.search(r"run_(\d+)_", os.path.basename(path)).group(1))
        frame = pd.read_csv(path)
        frame["run"] = run
        frames.append(frame)

    events = pd.concat(frames, ignore_index=True)
    scenario._events = events
    return events


def load_person_locations(scenario: Scenario, run: int) -> pd.DataFrame:
    """Where every Person was at every hour of one run. Needs --write_person_locations."""
    path = scenario.file("runs", f"run_{run}_person_locations.csv")
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"Missing {path}. Rerun the simulation with --write_person_locations.")
    return pd.read_csv(path)


def load_households(scenario: Scenario) -> pd.DataFrame:
    """Person -> home. Columns: Person_ID, Household_ID."""
    return pd.read_csv(scenario.file("household_id.csv"))


def load_workplaces(scenario: Scenario) -> pd.DataFrame:
    """Person -> workplace, only the persons that have one."""
    frame = pd.read_csv(scenario.file("work_id.csv"))
    return frame.dropna(subset=["Work_ID"]).astype({"Work_ID": int})


def load_locations(scenario: Scenario) -> pd.DataFrame:
    """Location -> type. Columns: Location_ID, Location_Type."""
    return pd.read_csv(scenario.file("location_id_and_type.csv"))


def load_city_config(scenario: Scenario) -> Dict[str, str]:
    """The derived city statistics as a plain dict."""
    frame = pd.read_csv(scenario.file("city_config.csv"))
    return dict(zip(frame["key"], frame["value"]))


def population_of(scenario: Scenario) -> int:
    return int(load_city_config(scenario)["total_population"])


def simulated_days(scenario: Scenario) -> float:
    time, _ = load_total_infections(scenario)
    return float(time[-1])


# --------------------------------------------------------------------------------------------
# Small helpers shared by the plots
# --------------------------------------------------------------------------------------------

def per_run_curves(scenario: Scenario) -> np.ndarray:
    """Cumulative infections of every single run, shaped (runs, time)."""
    curves = []
    for run in range(num_runs(scenario)):
        _, total = load_total_infections(scenario, f"run{run}")
        curves.append(total)
    return np.array(curves)


def events_per_run_by(scenario: Scenario, column: str, values: List[str]) -> pd.DataFrame:
    """Count events per run for each value of a column. Rows: run, columns: values."""
    events = load_infection_events(scenario)
    counts = (events.groupby(["run", column]).size().unstack(fill_value=0))
    for value in values:
        if value not in counts.columns:
            counts[value] = 0
    return counts[values]


def new_infections_grid(scenario: Scenario, locations: List[str], num_hours: int) -> np.ndarray:
    """Mean new infections per hour and location type, shaped (locations, hours)."""
    events = load_infection_events(scenario)
    events = events[events["Location_Type"].isin(locations)]
    runs = max(1, num_runs(scenario))

    grid = np.zeros((len(locations), num_hours))
    index = {name: i for i, name in enumerate(locations)}
    for (location, hour), count in events.groupby(["Location_Type", "Timestep"]).size().items():
        if 0 <= hour < num_hours:
            grid[index[location], int(hour)] += count
    return grid / runs


def savefig(figure, args: argparse.Namespace, name: str) -> str:
    os.makedirs(args.output, exist_ok=True)
    path = os.path.join(args.output, f"{name}.{args.format}")
    figure.savefig(path, dpi=args.dpi, bbox_inches="tight")
    print(f"  wrote {path}")
    return path
