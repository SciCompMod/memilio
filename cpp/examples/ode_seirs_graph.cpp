/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Henrik Zunker
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
#include "ode_seirs/model.h"
#include "ode_seirs/infection_state.h"
#include "ode_seirs/parameters.h"
#include "memilio/compartments/simulation.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/time_series.h"

#include "memilio/mobility/metapopulation_mobility_instant.h"
#include "memilio/compartments/simulation.h"

auto setup_model()
{
    mio::oseirs::Model<ScalarType> model(1);

    ScalarType total_population                                                          = 10000;
    model.populations[{mio::AgeGroup(0), mio::oseirs::InfectionState::Exposed}]          = 100;
    model.populations[{mio::AgeGroup(0), mio::oseirs::InfectionState::ExposedIsolated}]  = 10;
    model.populations[{mio::AgeGroup(0), mio::oseirs::InfectionState::Infected}]         = 100;
    model.populations[{mio::AgeGroup(0), mio::oseirs::InfectionState::Infectedisolated}] = 10;
    model.populations[{mio::AgeGroup(0), mio::oseirs::InfectionState::Recovered}]        = 100;
    model.populations[{mio::AgeGroup(0), mio::oseirs::InfectionState::Susceptible}] =
        total_population - model.populations[{mio::AgeGroup(0), mio::oseirs::InfectionState::Exposed}] -
        model.populations[{mio::AgeGroup(0), mio::oseirs::InfectionState::Infected}] -
        model.populations[{mio::AgeGroup(0), mio::oseirs::InfectionState::Recovered}];

    model.parameters.set<mio::oseirs::TimeExposed<ScalarType>>(5.2);
    model.parameters.set<mio::oseirs::TimeInfected<ScalarType>>(6);
    model.parameters.set<mio::oseirs::TimeImmunity<ScalarType>>(30);
    model.parameters.set<mio::oseirs::TransmissionProbabilityOnContact<ScalarType>>(0.1);
    model.parameters.set<mio::oseirs::ShareContagionsIsolated<ScalarType>>(0.1);
    model.parameters.set<mio::oseirs::TestAndTraceCapacity<ScalarType>>(10);
    model.parameters.set<mio::oseirs::DetectionRateExposedMinRisk<ScalarType>>(0.1);
    model.parameters.set<mio::oseirs::DetectionRateExposedMaxRisk<ScalarType>>(0.2);
    model.parameters.set<mio::oseirs::DetectionRateInfectedMinRisk<ScalarType>>(0.1);
    model.parameters.set<mio::oseirs::DetectionRateInfectedMaxRisk<ScalarType>>(0.2);

    // mio::ContactMatrixGroup& contact_matrix = model.parameters.get<mio::oseirs::ContactPatterns<ScalarType>>();
    // contact_matrix[0].get_baseline().setConstant(2.7);
    // contact_matrix[0].add_damping(0.7, mio::SimulationTime(30.));

    model.check_constraints();

    return model;
}

int main()
{
    mio::set_log_level(mio::LogLevel::debug);

    ScalarType t0   = 0;
    ScalarType tmax = 50.;
    ScalarType dt   = 0.5;

    mio::log_info("Simulating ODE SEIR; t={} ... {} with dt = {}.", t0, tmax, dt);

    auto model1 = setup_model();
    auto model2 = setup_model();

    model1.populations[{mio::AgeGroup(0), mio::oseirs::InfectionState::Exposed}]  = 250;
    model1.populations[{mio::AgeGroup(0), mio::oseirs::InfectionState::Infected}] = 500;

    // mio::Graph<mio::Simulation<double, mio::oseirs::Model<double>>, mio::MobilityEdge<ScalarType>> g;
    mio::Graph<mio::SimulationNode<mio::Simulation<double, mio::oseirs::Model<double>>>, mio::MobilityEdge<ScalarType>>
        g;

    g.add_node(1001, model1, t0);
    g.add_node(1002, model2, t0);
    g.add_edge(0, 1, Eigen::VectorXd::Constant((size_t)mio::oseirs::InfectionState::Count, 0.1));
    g.add_edge(1, 0, Eigen::VectorXd::Constant((size_t)mio::oseirs::InfectionState::Count, 0.1));

    auto sim = mio::make_mobility_sim(t0, dt, std::move(g));

    sim.advance(tmax);

    const auto& results_m1 = sim.get_graph().nodes()[0].property.get_result();

    results_m1.print_table({"S", "E", "E_I", "I", "I_I", "R"});
}
