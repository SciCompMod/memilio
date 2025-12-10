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

from memilio.simulation.sseirvv import InfectionState as State
from memilio.simulation.sseirvv import (Model, simulate_stochastic)


def run_sde_seirvv_simulation(params, N0=180000, T1=400):
    """
    Perform a forward SEIR simulation with mutation given a list of parameters.

    Parameters
    ----------
    params : list
        [gamma_inv, kappa_inv, beta, t1, I0]
    N0 : int, optional
        Initial population size (default: 180000)
    T1 : int, optional
        Total simulation time (default: 400)

    Returns
    -------
    tuple of np.ndarray
        Time series for compartments:
        (E_wt, S, R_var, E_both, I_wt, I_both, E_var, R_both, R_wt, I_var, N_t)
    """
    var_infc = 100
    dt = 1
    dt_inv = 1

    gamma_inv, kappa_inv, beta, t1, I0 = params
    params.append(var_infc)  

    # Initialize compartment arrays
    compartments = [np.zeros((T1 + 1, 1)) for _ in range(10)]
    E_wt, S, R_var, E_both, I_wt, I_both, E_var, R_both, R_wt, I_var = compartments

    # Initialize model
    model = Model()
    model.parameters.TimeExposedV1.value = kappa_inv
    model.parameters.TimeExposedV2.value = kappa_inv
    model.parameters.TimeInfectedV1.value = gamma_inv
    model.parameters.TimeInfectedV2.value = gamma_inv * 1.35
    model.parameters.TransmissionProbabilityOnContactV1.value = beta
    model.parameters.TransmissionProbabilityOnContactV2.value = beta

    # Initial populations
    model.populations[State.InfectedV1] = I0
    model.populations.set_difference_from_total(State.Susceptible, N0)
    model.parameters.ContactPatterns.baseline = np.ones((1, 1))

    # First phase simulation
    result = simulate_stochastic(0.0, t1, dt, model)
    for i in range(t1):
        vals = result.get_value(dt_inv * i)
        S[i], E_wt[i], I_wt[i], R_wt[i], E_var[i], I_var[i], R_var[i], E_both[i], I_both[i], R_both[i] = vals

    # Update populations at t1
    last_values = result.get_value(result.get_num_time_points() - 1)
    states = [
        State.Susceptible, State.ExposedV1, State.InfectedV1, State.RecoveredV1,
        State.ExposedV2, State.InfectedV2, State.RecoveredV2,
        State.ExposedV1V2, State.InfectedV1V2, State.RecoveredV1V2,
    ]
    for state, val in zip(states, last_values):
        model.populations[state] = val

    # Introduce variant infections
    model.populations[State.InfectedV2] = var_infc

    # Second phase simulation
    result = simulate_stochastic(t1, T1, dt, model)
    for i in range(T1 - t1 + 1):
        vals = result.get_value(dt_inv * i)
        idx = i + t1
        S[idx], E_wt[idx], I_wt[idx], R_wt[idx], E_var[idx], I_var[idx], R_var[idx], E_both[idx], I_both[idx], R_both[idx] = vals

    # Total population over time
    arr = np.array([E_wt, S, R_var, E_both, I_wt, I_both, E_var, R_both, R_wt, I_var])
    N_t = arr.sum(axis=0)

    return E_wt, S, R_var, E_both, I_wt, I_both, E_var, R_both, R_wt, I_var, N_t


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(
        'sde_sirs_simple',
        description='Simple example demonstrating the setup and simulation of the SDE SEIRVV model.')
    args = arg_parser.parse_args()
    params = [17,5,0.08,150,500]
    run_sde_seirvv_simulation(params)
