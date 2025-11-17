#############################################################################
# Copyright (C) 2020-2025 MEmilio
#
# Authors: Henrik Zunker
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
from memilio.simulation.oseirdb import (Model, InfectionState as State,
                                        interpolate_simulation_result,
                                        simulate)


def run_ode_seirdb_simulation():
    """Runs the C++ ODE SEIRDB model"""

    days = 60.
    dt = 1.0

    # Initialize Parameters
    num_groups = 1
    model = Model(num_groups)
    A0 = AgeGroup(0)

    total_population = 10000

    # Initial number of people in each compartment
    model.populations[A0, State.Exposed] = 100
    model.populations[A0, State.Infected] = 100
    model.populations[A0, State.Recovered] = 100
    model.populations[A0, State.Dead] = 100
    model.populations[A0, State.Buried] = 100
    model.populations.set_difference_from_total(
        (A0, State.Susceptible), total_population)

    # Compartment transition durations
    model.parameters.TimeExposed[A0] = 5.2
    model.parameters.TimeInfected[A0] = 6.0
    model.parameters.TimeToBurial[A0] = 5.5

    # Transition probabilities
    model.parameters.TransmissionProbabilityOnContact[A0] = 0.1
    model.parameters.TransmissionProbabilityFromDead[A0] = 0.01
    model.parameters.ProbabilityToRecover[A0] = 0.75

    model.parameters.ContactPatterns.cont_freq_mat[0].baseline = np.full(
        (num_groups, num_groups), 2.7)
    model.parameters.ContactPatterns.cont_freq_mat[0].minimum = np.zeros(
        (num_groups, num_groups))

    model.check_constraints()

    result = simulate(0.0, days, dt, model)
    result = interpolate_simulation_result(result)

    result.print_table(
        False, ["S", "E", "I", "R", "D", "B"], 16, 5)


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(
        'ode_seirdb',
        description='Simple example demonstrating the setup and simulation of the ODE SEIRDB model.')
    args = arg_parser.parse_args()
    run_ode_seirdb_simulation()
