#############################################################################
# Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
#
# Authors: 
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
import memilio.simulation as mio
import memilio.simulation.secir as secir
import numpy as np
import argparse

def parameter_study():
    mio.set_log_level(mio.LogLevel.Warning)

    t0 = 0
    tmax = 50

    #setup basic parameters
    model = secir.SecirModel(1)

    model.parameters.IncubationTime[secir.AgeGroup(0)] = 5.2
    model.parameters.SerialInterval[secir.AgeGroup(0)] = 4.2
    model.parameters.TimeInfectedSymptoms[secir.AgeGroup(0)] = 6
    model.parameters.HospitalizedToHomeTime[secir.AgeGroup(0)] = 12
    model.parameters.HospitalizedToICUTime[secir.AgeGroup(0)] = 2
    model.parameters.ICUToHomeTime[secir.AgeGroup(0)] = 8
    model.parameters.ICUToDeathTime[secir.AgeGroup(0)] = 5

    model.parameters.ContactPatterns.cont_freq_mat[0].baseline = np.r_[0.5]
    model.parameters.ContactPatterns.cont_freq_mat[0].add_damping(mio.Damping(np.r_[0.3], t = 0.3))

    model.parameters.InfectionProbabilityFromContact[secir.AgeGroup(0)] = 1.0
    model.parameters.RelativeCarrierInfectability[secir.AgeGroup(0)] = 0.67
    model.parameters.AsymptomaticCasesPerInfectious[secir.AgeGroup(0)] = 0.09
    model.parameters.RiskOfInfectionFromSymptomatic[secir.AgeGroup(0)] = 0.25
    model.parameters.HospitalizedCasesPerInfectious[secir.AgeGroup(0)] = 0.2
    model.parameters.ICUCasesPerHospitalized[secir.AgeGroup(0)] = 0.25
    model.parameters.DeathsPerICU[secir.AgeGroup(0)] = 0.3

    #two regions with different populations and with some migration between them
    graph = secir.MigrationGraph()
    model.populations[secir.AgeGroup(0), secir.InfectionState.Exposed] = 100
    model.populations[secir.AgeGroup(0), secir.InfectionState.Carrier] = 50
    model.populations[secir.AgeGroup(0), secir.InfectionState.Infected] = 50
    model.populations[secir.AgeGroup(0), secir.InfectionState.Hospitalized] = 20
    model.populations[secir.AgeGroup(0), secir.InfectionState.ICU] = 10 
    model.populations[secir.AgeGroup(0), secir.InfectionState.Recovered] = 10
    model.populations[secir.AgeGroup(0), secir.InfectionState.Dead] = 0
    model.populations.set_difference_from_group_total_AgeGroup(secir.Index_Agegroup_InfectionState(
        secir.AgeGroup(0), secir.InfectionState.Susceptible), 10000)
    model.apply_constraints()
    graph.add_node(id = 0, model = model, t0 = t0) #copies the model into the graph
    model.populations[secir.AgeGroup(0), secir.InfectionState.Exposed] = 0
    model.populations[secir.AgeGroup(0), secir.InfectionState.Carrier] = 0
    model.populations[secir.AgeGroup(0), secir.InfectionState.Infected] = 0
    model.populations[secir.AgeGroup(0), secir.InfectionState.Hospitalized] = 0
    model.populations[secir.AgeGroup(0), secir.InfectionState.ICU] = 0
    model.populations[secir.AgeGroup(0), secir.InfectionState.Recovered] = 0
    model.populations[secir.AgeGroup(0), secir.InfectionState.Dead] = 0
    model.populations.set_difference_from_group_total_AgeGroup(secir.Index_Agegroup_InfectionState(
        secir.AgeGroup(0), secir.InfectionState.Susceptible), 2000)
    model.apply_constraints()
    graph.add_node(id = 1, model = model, t0 = t0)
    migration_coefficients = 0.1 * np.ones(8)
    migration_params = mio.MigrationParameters(migration_coefficients)
    graph.add_edge(0, 1, migration_params) #one coefficient per (age group x compartment)
    graph.add_edge(1, 0, migration_params) #directed graph -> add both directions so coefficients can be different
    
    #run simulation
    sim = secir.MigrationSimulation(graph, t0, dt = 0.5)
    sim.advance(tmax)

    #process results
    region0_result = sim.graph.get_node(0).property.result
    region1_result = sim.graph.get_node(1).property.result

if __name__ == "__main__":    
    arg_parser = argparse.ArgumentParser(
        'migration', 
        description = 'Example demonstrating the setup and simulation of a geographically resolved SECIR model with travel.')
    args = arg_parser.parse_args()
    parameter_study()
