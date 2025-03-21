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

from memilio.simulation import Damping, AgeGroup
from memilio.simulation.oseir import Index_InfectionState
from memilio.simulation.oseir import InfectionState as State
from memilio.simulation.oseir import (Model, interpolate_simulation_result,
                                      simulate_flows)


def run_ode_seir_flows_simulation():
    """Runs the c++ ODE SEIR model using a flow simulation"""

    # Define population of age groups
    populations = [83000]

    days = 100  # number of days to simulate
    dt = 0.1

    # Initialize Parameters
    model = Model(1)

    A0 = AgeGroup(0)

    # Compartment transition duration
    model.parameters.TimeExposed[A0] = 5.2
    model.parameters.TimeInfected[A0] = 6.

    # Compartment transition propabilities
    model.parameters.TransmissionProbabilityOnContact[A0] = 1.

    # Initial number of people in each compartment
    model.populations[A0, State.Exposed] = 100
    model.populations[A0, State.Infected] = 50
    model.populations[A0, State.Recovered] = 10
    model.populations.set_difference_from_total(
        (A0, State.Susceptible), populations[0])

    # model.parameters.ContactPatterns = ContactMatrix(np.r_[0.5])
    model.parameters.ContactPatterns.cont_freq_mat[0].baseline = np.ones(
        (1, 1))
    model.parameters.ContactPatterns.cont_freq_mat[0].minimum = np.zeros(
        (1, 1))
    model.parameters.ContactPatterns.cont_freq_mat.add_damping(
        Damping(coeffs=np.r_[0.9], t=30.0, level=0, type=0))

    # Check logical constraints to parameters
    model.check_constraints()

    # Run flow simulation
    (result, flows) = simulate_flows(0, days, dt, model)

    print(result.print_table(["S", "E", "I", "R"], 16, 5))
    print(flows.print_table(["S->E", "E->I", "I->R"], 16, 5))


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(
        'ode_seir_flows',
        description='Simple example demonstrating the setup and flow simulation of the ODE SEIR model.')
    args = arg_parser.parse_args()
    run_ode_seir_flows_simulation()
