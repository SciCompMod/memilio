#############################################################################
# Copyright (C) 2020-2026 MEmilio
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
"""
Simple example demonstrating the ODE meningitis model simulation,
mirroring the C++ example in cpp/examples/ode_meng.cpp.
"""

import argparse

import numpy as np

from memilio.simulation import AgeGroup
from memilio.simulation.omeng import (InfectionState as State, Model,
                                      interpolate_simulation_result,
                                      simulate_flows)


def run_ode_meng_simulation():
    """Runs the ODE meningitis model and prints compartment and flow results."""

    t0 = 0.
    tmax = 10.
    dt = 0.1

    nb_total_t0 = 10000.
    nb_inf_t0 = 100.
    cont_freq = 10.

    # One age group
    model = Model(1)
    A0 = AgeGroup(0)

    # --- Parameters ---
    model.parameters.RateCarrierToInfected[A0] = 0.00022   # sigma
    model.parameters.RateCarrierToRecovered[A0] = 0.8      # eta_2
    model.parameters.RateInfectedToRecovered[A0] = 0.43    # eta_1
    model.parameters.RateInfectedToDead[A0] = 0.495        # d
    model.parameters.RateNaturalDeath[A0] = 0.0152207      # mu
    model.parameters.RateImmunityLoss[A0] = 0.851          # xi
    model.parameters.ProbabilityImmunityLossSusLow[A0] = 0.6997  # theta
    model.parameters.ModificationRate[A0] = 0.23           # a
    model.parameters.IncomeFractionSusLow[A0] = 0.585      # Delta
    model.parameters.IncomeRate[A0] = 19787.               # Pi
    model.parameters.RiskOfInfectionFromFromCarrier[A0] = 0.742   # omega
    model.parameters.RiskOfInfectionFromFromInfected[A0] = 0.425  # omega1
    model.parameters.TransmissionProbabilityOnContact[A0] = 0.05

    # Contact matrix: constant contact frequency
    model.parameters.ContactPatterns.cont_freq_mat[0].baseline = np.full(
        (1, 1), cont_freq)
    model.parameters.ContactPatterns.cont_freq_mat[0].minimum = np.zeros(
        (1, 1))

    # --- Initial conditions ---
    model.populations[A0, State.Incoming] = 0.
    model.populations[A0, State.SusceptibleHigh] = 0.
    model.populations[A0, State.Carrier] = 0.
    model.populations[A0, State.Infected] = nb_inf_t0
    model.populations[A0, State.Recovered] = 0.
    model.populations[A0, State.Dead] = 0.
    model.populations[A0, State.DeadNatural] = 0.
    model.populations.set_difference_from_total(
        (A0, State.SusceptibleLow), nb_total_t0 - nb_inf_t0)

    model.check_constraints()

    # --- Simulation ---
    # simulate_flows returns a pair: [0] = compartments, [1] = flows.
    # simulate_flows uses the model-specific Simulation class internally, which
    # overrides advance() to include dynamic NPI checks (contact dampings).
    (result, flows) = simulate_flows(t0, tmax, dt, model)

    # Interpolate to integer day time points (t = 0, 1, 2, ..., 10)
    result_interp = interpolate_simulation_result(result)
    flows_interp = interpolate_simulation_result(flows)

    # Note: The "Incoming" compartment (Inc) acts as a pure source node for
    # births (rate Pi). Since the FlowModel framework only supports conservative
    # flows between compartments, there is no external inflow into Incoming.
    # Its value therefore becomes increasingly negative over time:
    #   dInc/dt = -(1-Delta)*Pi - Delta*Pi = -Pi
    # This is expected model behaviour: Inc accumulates the total number of
    # births with a negative sign and does NOT represent a real population.
    # The living population is S_L + S_H + C + I + R.

    print("Compartments:")
    print(result_interp.print_table(
        False, ["Inc", "S_L", "S_H", "C", "I", "R", "D_D", "D_N"], 16, 5))

    # print("\nFlows:")
    # flows_interp.print_table(
    #     False,
    #     ["Inc->S_H", "Inc->S_L", "S_H->C", "S_L->C", "C->I", "C->R",
    #      "I->R", "I->D_D", "R->S_L", "R->S_H",
    #      "S_H->D_N", "S_L->D_N", "C->D_N", "I->D_N", "R->D_N"],
    #     16, 5)


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(
        'ode_meng_simple',
        description='Simple example demonstrating the ODE meningitis model.')
    arg_parser.parse_args()
    run_ode_meng_simulation()
