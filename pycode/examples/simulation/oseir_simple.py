#############################################################################
# Copyright (C) 2020-2024 MEmilio
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
import argparse

import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigvals
from scipy.sparse.linalg import eigs

from memilio.simulation import Damping
from memilio.simulation.oseir import Index_InfectionState
from memilio.simulation.oseir import InfectionState as State
from memilio.simulation.oseir import (Model, interpolate_simulation_result,
                                      simulate)


def stiffness(params, pop_total, t_idx, y):
    if not (t_idx < y.get_num_time_points()):
        raise ValueError("t_idx is not a valid index for the TimeSeries")

    coeffStoE = params.ContactPatterns.get_matrix_at(
        y.get_time(t_idx))[0, 0] * params.TransmissionProbabilityOnContact.value / pop_total
    Jacobian = np.zeros((4,4))

    S = y.get_value(t_idx)[0]
    I = y.get_value(t_idx)[2]

    Jacobian[0,0] = - coeffStoE * I
    Jacobian[0,2] = - coeffStoE * S
    Jacobian[1,0] = coeffStoE * I
    Jacobian[1,1] = - (1.0 / params.TimeExposed.value)
    Jacobian[1,2] = coeffStoE * S
    Jacobian[2,1] = (1.0 / params.TimeExposed.value)
    Jacobian[2,2] = -(1.0 / params.TimeInfected.value)
    Jacobian[3,2] = (1.0 / params.TimeInfected.value)

    eigen_vals = eigvals(Jacobian)
    eigen_vals_abs = np.abs(eigen_vals)

    eigen_vals_abs2 = [val for val in eigen_vals_abs if val >= 1e-12]

    max_val = max(eigen_vals_abs2)
    min_val = min(eigen_vals_abs2)

    xi = max_val / min_val
    # print("Max eigenvalue: ", max_val)
    # print("Min eigenvalue: ", min_val)
    return xi


def run_oseir_simulation():
    """
    Runs the c++ oseir model
    """

    # Define population of age groups
    populations = [83000]

    days = 100  # number of days to simulate
    dt = 0.1

    # Initialize Parameters
    model = Model()

    # Compartment transition duration
    model.parameters.TimeExposed.value = 5.2
    model.parameters.TimeInfected.value = 6.

    # Compartment transition propabilities
    model.parameters.TransmissionProbabilityOnContact.value = 1.

    # Initial number of people in each compartment
    model.populations[Index_InfectionState(State.Exposed)] = 100
    model.populations[Index_InfectionState(State.Infected)] = 50
    model.populations[Index_InfectionState(State.Recovered)] = 10
    model.populations.set_difference_from_total(
        (Index_InfectionState(State.Susceptible)), populations[0])

    # model.parameters.ContactPatterns = ContactMatrix(np.r_[0.5])
    model.parameters.ContactPatterns.baseline = np.ones((1, 1))
    model.parameters.ContactPatterns.minimum = np.zeros((1, 1))
    model.parameters.ContactPatterns.add_damping(
        Damping(coeffs=np.r_[0.9], t=30.0, level=0, type=0))

    # Check logical constraints to parameters
    model.check_constraints()

    # Run Simulation
    result = simulate(0, days, dt, model)
    # interpolate results
    result = interpolate_simulation_result(result)

    print(result.get_last_value())

    time = []
    stiff_val = []
    for time_idx in range(result.get_num_time_points()):
        time.append(result.get_time(time_idx))
        stiff_val.append(stiffness(model.parameters, model.populations.get_total(), time_idx, result))

    # print time and stiffness values
    plt.plot(time, stiff_val)
    plt.xlabel("Time")
    plt.ylabel("\u03BE") # xi
    plt.grid()
    plt.yscale("log")
    plt.show()

    print("Min value of xi: ", min(stiff_val))
    


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(
        'oseir_simple',
        description='Simple example demonstrating the setup and simulation of the OSEIR model.')
    args = arg_parser.parse_args()
    run_oseir_simulation()
