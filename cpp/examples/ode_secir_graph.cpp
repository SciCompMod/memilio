/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Daniel Abele, Henrik Zunker
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
#include "memilio/config.h"
#include "ode_secir/model.h"
#include "ode_secir/infection_state.h"
#include "ode_secir/parameters.h"
#include "memilio/mobility/metapopulation_mobility_instant.h"
#include "memilio/compartments/simulation.h"
#include <vector>

#include <iostream>

int main()
{
    mio::set_log_level(mio::LogLevel::off);
    const auto t0   = 0.;
    const auto tmax = 2.;
    const auto dt   = 0.5; //time step of Mobility, daily Mobility every second step

    double cont_freq = 10; // see Polymod study

    double nb_total_t0 = 10000, nb_exp_t0 = 100, nb_inf_t0 = 50, nb_car_t0 = 50, nb_hosp_t0 = 20, nb_icu_t0 = 10,
           nb_rec_t0 = 10, nb_dead_t0 = 0;

    const auto num_groups = 6;
    mio::osecir::Model<double> model(num_groups);
    auto nb_groups = model.parameters.get_num_groups();
    double fact    = 1.0 / (double)(size_t)nb_groups;

    auto& params = model.parameters;

    params.set<mio::osecir::StartDay>(60);
    params.set<mio::osecir::Seasonality<double>>(0.2);
    params.get<mio::osecir::TestAndTraceCapacity<double>>() = 35;

    for (auto i = mio::AgeGroup(0); i < nb_groups; i++) {
        params.get<mio::osecir::TimeExposed<double>>()[i]            = 3.2;
        params.get<mio::osecir::TimeInfectedNoSymptoms<double>>()[i] = 2.;
        params.get<mio::osecir::TimeInfectedSymptoms<double>>()[i]   = 5.8;
        params.get<mio::osecir::TimeInfectedSevere<double>>()[i]     = 9.5;
        params.get<mio::osecir::TimeInfectedCritical<double>>()[i]   = 7.1;

        model.populations[{i, mio::osecir::InfectionState::Exposed}]                     = fact * nb_exp_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedNoSymptoms}]          = fact * nb_car_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedNoSymptomsConfirmed}] = 0;
        model.populations[{i, mio::osecir::InfectionState::InfectedSymptoms}]            = fact * nb_inf_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedSymptomsConfirmed}]   = 0;
        model.populations[{i, mio::osecir::InfectionState::InfectedSevere}]              = fact * nb_hosp_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedCritical}]            = fact * nb_icu_t0;
        model.populations[{i, mio::osecir::InfectionState::Recovered}]                   = fact * nb_rec_t0;
        model.populations[{i, mio::osecir::InfectionState::Dead}]                        = fact * nb_dead_t0;
        model.populations.set_difference_from_group_total<mio::AgeGroup>({i, mio::osecir::InfectionState::Susceptible},
                                                                         fact * nb_total_t0);

        params.get<mio::osecir::TransmissionProbabilityOnContact<double>>()[i]  = 0.05;
        params.get<mio::osecir::RelativeTransmissionNoSymptoms<double>>()[i]    = 0.7;
        params.get<mio::osecir::RecoveredPerInfectedNoSymptoms<double>>()[i]    = 0.09;
        params.get<mio::osecir::RiskOfInfectionFromSymptomatic<double>>()[i]    = 0.25;
        params.get<mio::osecir::MaxRiskOfInfectionFromSymptomatic<double>>()[i] = 0.45;
        params.get<mio::osecir::SeverePerInfectedSymptoms<double>>()[i]         = 0.2;
        params.get<mio::osecir::CriticalPerSevere<double>>()[i]                 = 0.25;
        params.get<mio::osecir::DeathsPerCritical<double>>()[i]                 = 0.3;
    }

    model.apply_constraints();

    mio::ContactMatrixGroup& contact_matrix = params.get<mio::osecir::ContactPatterns<double>>();
    contact_matrix[0] =
        mio::ContactMatrix(Eigen::MatrixXd::Constant((size_t)nb_groups, (size_t)nb_groups, fact * cont_freq));

    //two mostly identical groups
    auto model_group1 = model;
    //some contact restrictions in model_group1
    mio::ContactMatrixGroup& contact_matrix_m1 =
        model_group1.parameters.get<mio::osecir::ContactPatterns<ScalarType>>();
    contact_matrix_m1[0].add_damping(0.1, mio::SimulationTime(15.));

    //infection starts in group 1
    model_group1.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Susceptible}] = 9990;
    model_group1.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Exposed}]     = 100;

    mio::Graph<mio::SimulationNode<mio::osecir::Simulation<>>, mio::MobilityEdge<ScalarType>> g;

    const size_t num_nodes = 400;

    // Create nodes with varying parameters
    for (size_t i = 0; i < num_nodes; ++i) {
        auto model_node = model;

        // some dampings
        if (i < num_nodes / 10) {
            mio::ContactMatrixGroup& contact_matrix_node =
                model_node.parameters.get<mio::osecir::ContactPatterns<ScalarType>>();
            contact_matrix_node[0].add_damping(0.1, mio::SimulationTime(15.));
        }

        // Add initial infection to first node only
        if (i == 0) {
            model_node.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Susceptible}] = 9900;
            model_node.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Exposed}]     = 200;
        }

        // Add node to graph
        g.add_node(i, model_node, t0);
    }

    // Create edges between nodes
    for (size_t i = 0; i < num_nodes; ++i) {
        for (size_t j = 0; j < num_nodes; ++j) {
            if (i != j) {
                g.add_edge(i, j,
                           Eigen::VectorXd::Constant((size_t)mio::osecir::InfectionState::Count * num_groups, 0.001));
            }
        }
    }

    auto sim = mio::make_mobility_sim(t0, dt, std::move(g));

    sim.advance(tmax);

    std::vector<double> a(1000000, 1.0);
    std::vector<double> b(1000000, 1.0);
    std::vector<double> c(1000000);
    #pragma acc parallel loop
    for (size_t i = 0; i < a.size(); i++){
        c[i] = a[i] + b[i];
    }

    double result = 0;
    for (size_t i = 0; i < c.size(); i++){
        result += c[i];
    }
    std::cout << result << std::endl;
    return 0;
}
