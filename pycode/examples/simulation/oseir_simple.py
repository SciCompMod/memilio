#############################################################################
# Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
#
# Authors: Martin J. Kuehn, Wadim Koslow, Daniel Abele
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
from memilio.simulation import Damping
from memilio.simulation.oseir import Model, simulate, Index_InfectionState
from memilio.simulation.oseir import InfectionState as State
import numpy as np
import argparse


def run_secir_simulation(show_plot = True):
    """
    Runs the c++ secir model using one age group 
    and plots the results
    """

    # Define Comartment names
    compartments = ['Susceptible', 'Exposed', 'Infected', 'Recovered']
    # Define population of age groups
    populations = [83000]

    days = 100  # number of days to simulate
    dt = 0.1
    num_compartments = len(compartments)

    # Initialize Parameters
    model = Model()

    # Compartment transition duration
    model.parameters.LatentTime.value = 5.2
    model.parameters.InfectiousTime.value = 6.

    # Compartment transition propabilities
    model.parameters.InfectionProbabilityFromContact.value = 1.

    # Initial number of people in each compartment
    model.populations[Index_InfectionState(State.Exposed)] = 100
    model.populations[Index_InfectionState(State.Infected)] = 50
    model.populations[Index_InfectionState(State.Recovered)] = 10
    model.populations.set_difference_from_total(
        (Index_InfectionState(State.Susceptible)), populations[0])

    # model.parameters.ContactPatterns = ContactMatrix(np.r_[0.5])
    model.parameters.ContactPatterns.baseline = np.ones((1, 1)) * 1
    model.parameters.ContactPatterns.minimum = np.ones((1, 1)) * 0
    model.parameters.ContactPatterns.add_damping(
        Damping(coeffs=np.r_[0.9], t=30.0, level=0, type=0))

    # Apply mathematical constraints to parameters
    model.apply_constraints()

    # Run Simulation
    result = simulate(0, days, dt, model)
    print(result.get_last_value())


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(
        'secir_simple', 
        description = 'Simple example demonstrating the setup and simulation of the OSEIR model.')
    arg_parser.add_argument('-p', '--show_plot', action='store_const', const=True, default=False)
    args = arg_parser.parse_args()
    run_secir_simulation(**args.__dict__)
