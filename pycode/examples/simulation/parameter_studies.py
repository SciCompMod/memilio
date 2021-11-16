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
# import epidemiology.secir as secir
import numpy as np
import memilio.simulation as mio
import memilio.simulation.secir as secir

def parameter_study():    
    #setup basic parameters
    num_groups = 6
    model = secir.SecirModel(num_groups)

    for i in range(num_groups):
        group = secir.AgeGroup(i)

        model.parameters.IncubationTime[group] = 5.2
        model.parameters.InfectiousTimeMild[group] = 6
        model.parameters.SerialInterval[group] = 4.2
        model.parameters.HospitalizedToHomeTime[group] = 12
        model.parameters.HomeToHospitalizedTime[group] = 5
        model.parameters.HospitalizedToICUTime[group] = 2
        model.parameters.ICUToHomeTime[group] = 8
        model.parameters.ICUToDeathTime[group] = 5

        model.populations[group, secir.Index_InfectionState(secir.InfectionState.Exposed)] = 100
        model.populations[group, secir.Index_InfectionState(secir.InfectionState.Carrier)] = 50
        model.populations[group, secir.Index_InfectionState(secir.InfectionState.Infected)] = 20
        model.populations[group, secir.Index_InfectionState(secir.InfectionState.Hospitalized)] = 20
        model.populations[group, secir.Index_InfectionState(secir.InfectionState.ICU)] = 10
        model.populations[group, secir.Index_InfectionState(secir.InfectionState.Recovered)] = 50
        model.populations[group, secir.Index_InfectionState(secir.InfectionState.Dead)] = 10
        model.populations.set_difference_from_group_total_AgeGroup((group, secir.Index_InfectionState(secir.InfectionState.Susceptible)), 10000)


        model.parameters.InfectionProbabilityFromContact[group].set_distribution(mio.ParameterDistributionUniform(0.1, 0.2))
        model.parameters.AsymptoticCasesPerInfectious[group] = 0.09
        model.parameters.RiskOfInfectionFromSympomatic[group] = 0.25
        model.parameters.HospitalizedCasesPerInfectious[group] = 0.2
        model.parameters.ICUCasesPerHospitalized[group] = 0.25
        model.parameters.DeathsPerHospitalized[group] = 0.3


    model.parameters.ContactPatterns.cont_freq_mat = mio.ContactMatrixGroup(4,num_groups)
    model.parameters.ContactPatterns.cont_freq_mat[0] = mio.ContactMatrix(np.ones((num_groups, num_groups))*0.5)
    model.parameters.ContactPatterns.cont_freq_mat[1] = mio.ContactMatrix(np.ones((num_groups, num_groups))*0.5)
    model.parameters.ContactPatterns.cont_freq_mat[2] = mio.ContactMatrix(np.ones((num_groups, num_groups))*0.5)
    model.parameters.ContactPatterns.cont_freq_mat[3] = mio.ContactMatrix(np.ones((num_groups, num_groups))*0.5)
    model.parameters.ContactPatterns.cont_freq_mat.add_damping(mio.Damping(np.ones((num_groups, num_groups))*0.7, 30.0))
    print(model.parameters.ContactPatterns.cont_freq_mat[1].baseline)
    #process the result of one run
    parameter_study.c = 0
    def handle_result(graph):

        group = secir.AgeGroup(0)
        print("run {} with infection rate {:.2G}".format(parameter_study.c, graph.get_node(0).property.model.parameters.InfectionProbabilityFromContact[group].value))
        print("compartments at t = {}:".format(graph.get_node(0).property.result.get_time(0)))
        print(graph.get_node(0).property.result.get_value(0))
        print("compartments at t = {}:".format(graph.get_node(0).property.result.get_last_time()))
        print(graph.get_node(0).property.result.get_last_value())
        parameter_study.c += 1

    #study the effect of different infection rates

    model.apply_constraints()

    graph = secir.SecirModelGraph()
    graph.add_node(0, model)
    graph.add_node(1, model)
    graph.add_edge(0, 1, 0.01 * np.ones(8*num_groups))
    graph.add_edge(1, 0, 0.01 * np.ones(8*num_groups))

    study = secir.ParameterStudy(graph, t0 = 1, tmax = 10, dt = 0.5, num_runs = 3)
    study.run(handle_result)

if __name__ == "__main__":
    parameter_study()