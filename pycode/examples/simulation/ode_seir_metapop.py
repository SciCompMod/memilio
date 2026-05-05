#############################################################################
# Copyright (C) 2020-2026 MEmilio
#
# Authors: Carlotta Gerstein
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

import numpy as np

import memilio.simulation as mio
import memilio.simulation.oseir_metapop as oseir_metapop
from memilio.simulation import AgeGroup, LogLevel, set_log_level
set_log_level(LogLevel.Off)

t0 = 0
tmax = 10
dt = 0.1

model = oseir_metapop.Model(3, 1)

for i in range(model.parameters.num_regions.get()):
    # Set infection state stay times (in days)
    model.populations[mio.Region(i), AgeGroup(
        0), oseir_metapop.InfectionState.Susceptible] = 10000

model.populations[mio.Region(0), AgeGroup(
    0), oseir_metapop.InfectionState.Infected].value += 100
model.populations[mio.Region(0), AgeGroup(
    0), oseir_metapop.InfectionState.Susceptible].value -= 100

mobility_data_commuter = np.array(
    [[0.4, 0.3, 0.3], [0.2, 0.7, 0.1], [0.4, 0.1, 0.5]])
model.set_commuting_strengths(mobility_data_commuter)

# Set contact frequency
model.parameters.ContactPatterns.cont_freq_mat[0].baseline = np.ones(
    (3, 3)) * 2.7

model.parameters.TimeExposed[mio.Region(0), AgeGroup(0)] = 3.
model.parameters.TimeExposed[mio.Region(1), AgeGroup(0)] = 4.
model.parameters.TimeExposed[mio.Region(2), AgeGroup(0)] = 5.
model.parameters.TimeInfected[mio.Region(0), AgeGroup(0)] = 7.
model.parameters.TimeInfected[mio.Region(1), AgeGroup(0)] = 8.
model.parameters.TimeInfected[mio.Region(2), AgeGroup(0)] = 9.
model.parameters.TransmissionProbabilityOnContact[mio.Region(
    0), AgeGroup(0)] = 0.07333
model.parameters.TransmissionProbabilityOnContact[mio.Region(
    1), AgeGroup(0)] = 0.07333
model.parameters.TransmissionProbabilityOnContact[mio.Region(
    2), AgeGroup(0)] = 0.07333

result = oseir_metapop.simulate(t0, tmax, dt, model)
interpolated_result = oseir_metapop.interpolate_simulation_result(result)

print("Infected individuals per Region over time [%]:")

col_width = 12

print("Time".ljust(
    col_width) + "".join(f"Region {i}".ljust(col_width) for i in range(model.parameters.num_regions.get())))

for t_idx in range(interpolated_result.get_num_time_points()):
    row = f"{interpolated_result.get_times()[t_idx]:.1f}".ljust(col_width)
    values = interpolated_result.get_value(t_idx)

    for i in range(model.parameters.num_regions.get()):
        population = model.populations.get_group_total_Region(mio.Region(i))
        percentage = values[model.parameters.num_regions.get()
                            * 3 + i] / population * 100
        row += f"{percentage:.5f}".ljust(col_width)

    print(row)
