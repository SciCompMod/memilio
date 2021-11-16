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

def parameter_study():    
    t0 = 0
    tmax = 50

    #setup basic parameters
    model = secir.SecirModel1()

    model.parameters.times[0].set_incubation(5.2)
    model.parameters.times[0].set_infectious_mild(6)
    model.parameters.times[0].set_serialinterval(4.2)
    model.parameters.times[0].set_hospitalized_to_home(12)
    model.parameters.times[0].set_home_to_hospitalized(5)
    model.parameters.times[0].set_hospitalized_to_icu(2)
    model.parameters.times[0].set_icu_to_home(8)
    model.parameters.times[0].set_icu_to_death(5)

    model.parameters.get_contact_patterns().get_cont_freq_mat().set_cont_freq(0.5, 0, 0)
    model.parameters.get_contact_patterns().get_cont_freq_mat().add_damping(mio.Damping(30, 0.3), 0, 0)

    model.parameters.probabilities[0].set_infection_from_contact(1.0)
    model.parameters.probabilities[0].set_carrier_infectability(0.67)
    model.parameters.probabilities[0].set_asymp_per_infectious(0.09)
    model.parameters.probabilities[0].set_risk_from_symptomatic(0.25)
    model.parameters.probabilities[0].set_hospitalized_per_infectious(0.2)
    model.parameters.probabilities[0].set_icu_per_hospitalized(0.25)
    model.parameters.probabilities[0].set_dead_per_icu(0.3)

    #two regions with different populations and with some migration between them
    graph = secir.MigrationGraph1()
    model.populations.set(100, secir.AgeGroup1.Group0, secir.InfectionType.E)
    model.populations.set(50, secir.AgeGroup1.Group0, secir.InfectionType.C)
    model.populations.set(50, secir.AgeGroup1.Group0, secir.InfectionType.I)
    model.populations.set(20, secir.AgeGroup1.Group0, secir.InfectionType.H)
    model.populations.set(10, secir.AgeGroup1.Group0, secir.InfectionType.U)
    model.populations.set(10, secir.AgeGroup1.Group0, secir.InfectionType.R)
    model.populations.set(0, secir.AgeGroup1.Group0, secir.InfectionType.D)
    model.populations.set_difference_from_total(10000, secir.AgeGroup1.Group0, secir.InfectionType.S)
    model.apply_constraints()
    graph.add_node(model, t0)
    model.populations.set(0, secir.AgeGroup1.Group0, secir.InfectionType.E)
    model.populations.set(0, secir.AgeGroup1.Group0, secir.InfectionType.C)
    model.populations.set(0, secir.AgeGroup1.Group0, secir.InfectionType.I)
    model.populations.set(0, secir.AgeGroup1.Group0, secir.InfectionType.H)
    model.populations.set(0, secir.AgeGroup1.Group0, secir.InfectionType.U)
    model.populations.set(0, secir.AgeGroup1.Group0, secir.InfectionType.R)
    model.populations.set(0, secir.AgeGroup1.Group0, secir.InfectionType.D)
    model.populations.set_difference_from_total(2000, secir.AgeGroup1.Group0, secir.InfectionType.S)
    model.apply_constraints()
    graph.add_node(model, t0)
    migration_coefficients = 0.1 * np.ones(8)
    migration_params = mio.MigrationParams(migration_coefficients)
    graph.add_edge(0, 1, migration_params) #one coefficient per (age group x compartment)
    graph.add_edge(1, 0, migration_params) #directed graph -> add both directions so coefficients can be different
    
    #run simulation
    sim = secir.MigrationSimulation1(graph, t0, dt = 0.5)
    sim.advance(tmax)

    #process results
    region0_result = sim.graph.get_node(0).result
    region1_result = sim.graph.get_node(1).result

if __name__ == "__main__":
    parameter_study()