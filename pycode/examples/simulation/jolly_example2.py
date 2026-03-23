#############################################################################
# Copyright (C) 2020-2026 MEmilio
#
# Authors: Kilian Volmer, Henrik Zunker
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
from memilio.simulation import jolly


def run_jolly_simulation(farm_file):
    """
    Runs the Jolly simulation (Stochastic Metapopulation Model) via Python bindings.
    """
    # Set parameters
    tmax = 20
    dt = 1.0
    suspicion_threshold = 0.2
    sensitivity = 0.95
    h0 = 0.0002
    r0 = 4000
    alpha = 10
    infection_baseline = 0.0001
    culling_factor = 1.0
    A0_SEI = 0.5
    A0_EI = 0.3
    A0_ID = 0.1
    A0_DeathRate = 0.0001
    A1_SEI = 0.5
    A1_EI = 0.3
    A1_ID = 0.1
    A1_DeathRate = 0.0001
    A2_SEI = 0.5
    A2_EI = 0.3
    A2_ID = 0.1
    A2_DeathRate = 0.0001
    A3_SEI = 0.5
    A3_EI = 0.3
    A3_ID = 0.1
    A3_DeathRate = 0.0001
    A4_SEI = 0.5
    A4_EI = 0.3
    A4_ID = 0.1
    A4_DeathRate = 0.0001
    foi_inner_factor0 = 1.0
    foi_inner_factor1 = 1.0
    foi_inner_factor2 = 1.0
    foi_inner_factor3 = 1.0
    foi_inner_factor4 = 1.0
    foi_outer_factor0 = 1.0
    foi_outer_factor1 = 1.0
    foi_outer_factor2 = 1.0
    foi_outer_factor3 = 1.0
    foi_outer_factor4 = 1.0
    damping0 = 0.5
    damping1 = 1.0
    damping2 = 1.0
    damping3 = 1.0
    damping4 = 1.0
    first_infection_day = 0
    second_infection_day = 2
    third_infection_day = 2
    seed = 434

    print("Starting Jolly Simulation...")

    sim_start = jolly.simulate_with_init(
        farm_file, tmax, dt, suspicion_threshold, sensitivity, h0, r0, alpha, infection_baseline, culling_factor, A0_SEI,
        A0_EI, A0_ID, A0_DeathRate, A1_SEI, A1_EI, A1_ID, A1_DeathRate, A2_SEI, A2_EI, A2_ID, A2_DeathRate, A3_SEI,
        A3_EI, A3_ID, A3_DeathRate, A4_SEI, A4_EI, A4_ID, A4_DeathRate, foi_inner_factor0, foi_outer_factor0,
        foi_inner_factor1, foi_outer_factor1, foi_inner_factor2, foi_outer_factor2, foi_inner_factor3,
        foi_outer_factor3, foi_inner_factor4, foi_outer_factor4, damping0, damping1, damping2, damping3, damping4,
        first_infection_day, second_infection_day, third_infection_day, seed)

    result = jolly.get_result(sim_start, tmax)

    infected_farms = sum(1 for r in result if r >= 0)
    print(f"Result: {infected_farms} farms are infected at {tmax}.")

    continued_tmax = 60

    sim_continued = jolly.simulate_continued(sim_start, tmax, continued_tmax, dt, suspicion_threshold, sensitivity, h0, r0, alpha, infection_baseline, culling_factor, A0_SEI,
                                             A0_EI, A0_ID, A0_DeathRate, A1_SEI, A1_EI, A1_ID, A1_DeathRate, A2_SEI, A2_EI, A2_ID, A2_DeathRate, A3_SEI,
                                             A3_EI, A3_ID, A3_DeathRate, A4_SEI, A4_EI, A4_ID, A4_DeathRate, foi_inner_factor0, foi_outer_factor0,
                                             foi_inner_factor1, foi_outer_factor1, foi_inner_factor2, foi_outer_factor2, foi_inner_factor3,
                                             foi_outer_factor3, foi_inner_factor4, foi_outer_factor4, damping0, damping1, damping2, damping3, damping4,
                                             first_infection_day, second_infection_day, third_infection_day, seed)

    result = jolly.get_result(sim_continued, continued_tmax)

    infected_farms = sum(1 for r in result if r >= 0)
    print(f"Result: {infected_farms} farms are infected at {continued_tmax}.")


if __name__ == "__main__":
    farm_file = "/hpc_data/bick_ju/jolly/farms.csv"
    run_jolly_simulation(farm_file)
