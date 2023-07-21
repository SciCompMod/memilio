/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
*
* Authors: Daniel Abele
*
* Contact: Martin J. Kuehn <Martin.Kuehn@DLR.de>
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
*     http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.
*/
#include "ode_secir/model.h"
#include "ode_secir/infection_state.h"
#include "ode_secir/parameters.h"
#include "memilio/compartments/simulation.h"
#include "memilio/mobility/graph.h"
#include "memilio/mobility/metapopulation_mobility_stochastic.h"
#include "memilio/mobility/metapopulation_mobility_instant.h"

#include <iostream>

int main(int argc, char** argv)
{
    const auto t0   = 0.;
    const auto tmax = 10.;
    const auto dt   = 0.1; //initial time step

    //total compartment sizes
    double num_total = 10000, num_exp = 200, num_ins = 50, num_is = 50, num_isev = 10, num_icri = 5, num_rec = 0,
           num_dead = 0;

    //model with 3 age groups
    mio::osecir::Model model(3);

    auto& params = model.parameters;

    auto num_age_groups = params.get_num_groups();
    double fact         = 1.0 / (double)(size_t)num_age_groups;

    //set initialization and model parameters
    for (auto i = mio::AgeGroup(0); i < num_age_groups; i++) {
        params.get<mio::osecir::IncubationTime>()[i]       = 5.2;
        params.get<mio::osecir::TimeInfectedSymptoms>()[i] = 6.;
        params.get<mio::osecir::SerialInterval>()[i]       = 4.2;
        params.get<mio::osecir::TimeInfectedSevere>()[i]   = 12;
        params.get<mio::osecir::TimeInfectedCritical>()[i] = 8;

        //initial populations is equally distributed among age groups
        model.populations[{i, mio::osecir::InfectionState::Exposed}]            = fact * num_exp;
        model.populations[{i, mio::osecir::InfectionState::InfectedNoSymptoms}] = fact * num_ins;
        model.populations[{i, mio::osecir::InfectionState::InfectedSymptoms}]   = fact * num_is;
        model.populations[{i, mio::osecir::InfectionState::InfectedSevere}]     = fact * num_isev;
        model.populations[{i, mio::osecir::InfectionState::InfectedCritical}]   = fact * num_icri;
        model.populations[{i, mio::osecir::InfectionState::Recovered}]          = fact * num_rec;
        model.populations[{i, mio::osecir::InfectionState::Dead}]               = fact * num_dead;
        model.populations.set_difference_from_group_total<mio::AgeGroup>({i, mio::osecir::InfectionState::Susceptible},
                                                                         fact * num_total);

        params.get<mio::osecir::TransmissionProbabilityOnContact>()[i] = 0.05;
        params.get<mio::osecir::RelativeTransmissionNoSymptoms>()[i]   = 0.67;
        params.get<mio::osecir::RecoveredPerInfectedNoSymptoms>()[i]   = 0.09;
        params.get<mio::osecir::RiskOfInfectionFromSymptomatic>()[i]   = 0.25;
        params.get<mio::osecir::SeverePerInfectedSymptoms>()[i]        = 0.2;
        params.get<mio::osecir::CriticalPerSevere>()[i]                = 0.25;
        params.get<mio::osecir::DeathsPerCritical>()[i]                = 0.3;
    }

    //add contact pattern and contact damping
    mio::ContactMatrixGroup& contact_matrix = params.get<mio::osecir::ContactPatterns>();
    contact_matrix[0] =
        mio::ContactMatrix(Eigen::MatrixXd::Constant((size_t)num_age_groups, (size_t)num_age_groups, fact * 10));
    contact_matrix.add_damping(Eigen::MatrixXd::Constant((size_t)num_age_groups, (size_t)num_age_groups, 0.6),
                               mio::SimulationTime(5.));

    model.apply_constraints();

    //modify model for second node
    auto model2 = model;

    for (auto i = mio::AgeGroup(0); i < num_age_groups; i++) {

        model2.populations[{i, mio::osecir::InfectionState::Exposed}]            = fact * 100;
        model2.populations[{i, mio::osecir::InfectionState::InfectedNoSymptoms}] = fact * 100;
        model2.populations[{i, mio::osecir::InfectionState::InfectedCritical}]   = fact * 10;
        model2.populations.set_difference_from_group_total<mio::AgeGroup>({i, mio::osecir::InfectionState::Susceptible},
                                                                          fact * num_total);
    }

    mio::Graph<mio::SimulationNode<mio::Simulation<mio::osecir::Model>>, mio::MigrationEdgeStochastic> graph;
    graph.add_node(1001, model, t0);
    graph.add_node(1002, model2, t0);

    auto transition_rates = mio::MigrationCoefficients(model.populations.numel());
    ScalarType kappa      = 0.01;

    for (auto age = mio::AgeGroup(0); age < num_age_groups; age++) {
        for (auto compartment = mio::Index<mio::osecir::InfectionState>(0);
             compartment < model.populations.size<mio::osecir::InfectionState>(); compartment++) {
            auto coeff_idx = model.populations.get_flat_index({age, compartment});
            //age group 0 commutes less
            if (age == mio::AgeGroup(0)) {
                transition_rates.get_baseline()[coeff_idx] = 0.01;
            }
            else {
                transition_rates.get_baseline()[coeff_idx] = 0.02;
            }
        }
    }
    //infected people commute less, hospitalized and dead do not commute
    for (auto age = mio::AgeGroup(0); age < num_age_groups; age++) {
        transition_rates
            .get_baseline()[model.populations.get_flat_index({age, mio::osecir::InfectionState::InfectedSymptoms})] *=
            0.8;
        transition_rates
            .get_baseline()[model.populations.get_flat_index({age, mio::osecir::InfectionState::InfectedSevere})] = 0.;
        transition_rates
            .get_baseline()[model.populations.get_flat_index({age, mio::osecir::InfectionState::InfectedCritical})] =
            0.;
        transition_rates.get_baseline()[model.populations.get_flat_index({age, mio::osecir::InfectionState::Dead})] =
            0.;
    }

    transition_rates.get_baseline() *= kappa;

    graph.add_edge(0, 1, std::move(transition_rates));
    graph.add_edge(1, 0, std::move(transition_rates));

    auto sim = mio::make_migration_sim(t0, dt, std::move(graph));

    sim.advance(tmax);

    return 0;
}
