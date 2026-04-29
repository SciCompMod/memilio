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
import argparse

import numpy as np

from memilio.simulation.omseirs4 import InfectionState as State
from memilio.simulation.omseirs4 import Model, simulate

# MSEIRS4 model (maternal immunity + 4 infection levels) as described in https://doi.org/10.1016/S0025-5564(01)00066-9.


def run_ode_omseirs4_simulation():
    # Create model and set parameters (per day)
    model = Model()
    params = model.parameters

    params.BaseTransmissionRate.value = 0.4
    params.SeasonalAmplitude.value = 0.15
    params.SeasonalPhase.value = 0.0
    params.NaturalBirthDeathRate.value = 1.0 / (70.0 * 365.0)
    params.LossMaternalImmunityRate.value = 1.0 / 90.0
    params.ProgressionRate.value = 1.0 / 7.0
    params.RecoveryRate.value = 1.0 / 14.0
    params.ImmunityWaningRate.value = 1.0 / (5.0 * 365.0)
    params.Beta2Factor.value = 0.5
    params.Beta3Factor.value = 0.35
    params.Beta4Factor.value = 0.25

    # Initial populations
    N = 1_000_000.0
    pop = model.populations

    pop[State.MaternalImmune] = 5000.0

    pop[State.E1] = 300.0
    pop[State.E2] = 150.0
    pop[State.E3] = 80.0
    pop[State.E4] = 70.0

    pop[State.I1] = 200.0
    pop[State.I2] = 100.0
    pop[State.I3] = 50.0
    pop[State.I4] = 50.0

    pop[State.R1] = 40000.0
    pop[State.R2] = 30000.0
    pop[State.R3] = 20000.0
    pop[State.R4] = 10000.0

    pop[State.S2] = 100000.0
    pop[State.S3] = 50000.0
    pop[State.S4] = 50000.0

    assigned = 5000.0 + (300.0 + 150.0 + 80.0 + 70.0) + (200.0 + 100.0 + 50.0 + 50.0) + \
        (40000.0 + 30000.0 + 20000.0 + 10000.0) + (100000.0 + 50000.0 + 50000.0)
    s1 = max(0.0, N - assigned)
    pop[State.S1] = s1

    model.check_constraints()

    # simulate 10 days with daily time steps
    t0 = 0.0
    tmax = 10.0
    dt = 1.0

    result = simulate(t0, tmax, dt, model)

    # print last results
    print("Last results:")
    print(f"Time: {result.get_last_time()}")
    print(f"compartments: {result.get_last_value()}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser("ode_omseirs4_simple",
                                     description="Simple example demonstrating the ODE MSEIRS4 RSV model")
    args = parser.parse_args()
    run_ode_omseirs4_simulation()
