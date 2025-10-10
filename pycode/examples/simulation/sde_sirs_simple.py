#############################################################################
# Copyright (C) 2020-2025 MEmilio
#
# Authors: Maximilian Betz
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
import argparse

import numpy as np

from memilio.simulation import AgeGroup, Damping
from memilio.simulation.ssirs import InfectionState as State
from memilio.simulation.ssirs import (
    Model, simulate_stochastic, interpolate_simulation_result)


def run_sde_sirs_simulation():
    """Runs SDE SIRS model"""

    tmax = 5.  # simulation time frame
    dt = 0.001

    # Initialize Model
    model = Model()

    # Mean time in Infected compartment
    model.parameters.TimeInfected.value = 10.
    model.parameters.TimeImmune.value = 100.

    model.parameters.TransmissionProbabilityOnContact.value = 1.

    # Initial number of people per compartment
    total_population = 10000
    model.populations[State.Infected] = 100
    model.populations[State.Recovered] = 1000
    model.populations.set_difference_from_total(
        (State.Susceptible), total_population)

    model.parameters.ContactPatterns.baseline = np.ones(
        (1, 1)) * 20.7
    model.parameters.ContactPatterns.minimum = np.zeros(
        (1, 1))
    model.parameters.ContactPatterns.add_damping(
        Damping(coeffs=np.r_[0.6], t=2, level=0, type=0))

    # Check parameter constraints
    model.check_constraints()

    # Run Simulation
    result = simulate_stochastic(0., days, dt, model)

    # Interpolate results
    result = interpolate_simulation_result(result)

    result.print_table(False, ["Susceptible", "Infected", "Recovered"], 16, 5)


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(
        'sde_sirs_simple',
        description='Simple example demonstrating the setup and simulation of the SDE SIRS model.')
    args = arg_parser.parse_args()
    run_sde_sirs_simulation()
