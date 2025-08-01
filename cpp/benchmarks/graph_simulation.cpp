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
#include "benchmarks/graph_simulation.h"
#include "memilio/compartments/simulation.h"
#include "memilio/mobility/metapopulation_mobility_instant.h"
#include "benchmark/benchmark.h"
#include "ode_secirvvs/model.h"
#include "memilio/math/adapt_rk.h"
#include <string>

const std::string config_path = "../../benchmarks/graph_simulation.config";

mio::osecirvvs::Model<double> create_model(size_t num_agegroups, const double tmax)
{
    mio::osecirvvs::Model<double> model(num_agegroups);
    const size_t pop_total = 10000;
    const size_t init_val  = 20;
    for (mio::AgeGroup i = 0; i < (mio::AgeGroup)num_agegroups; i++) {
        model.populations[{i, mio::osecirvvs::InfectionState::ExposedNaive}]                                = init_val;
        model.populations[{i, mio::osecirvvs::InfectionState::ExposedImprovedImmunity}]                     = init_val;
        model.populations[{i, mio::osecirvvs::InfectionState::ExposedPartialImmunity}]                      = init_val;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedNoSymptomsNaive}]                     = init_val;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedNoSymptomsNaiveConfirmed}]            = init_val;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedNoSymptomsPartialImmunity}]           = init_val;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedNoSymptomsPartialImmunityConfirmed}]  = init_val;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedNoSymptomsImprovedImmunity}]          = init_val;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedNoSymptomsImprovedImmunityConfirmed}] = init_val;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedSymptomsNaive}]                       = init_val;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedSymptomsNaiveConfirmed}]              = init_val;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedSymptomsPartialImmunity}]             = init_val;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedSymptomsPartialImmunityConfirmed}]    = init_val;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedSymptomsImprovedImmunity}]            = init_val;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedSymptomsImprovedImmunityConfirmed}]   = init_val;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedSevereNaive}]                         = init_val;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedSevereImprovedImmunity}]              = init_val;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedSeverePartialImmunity}]               = init_val;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedCriticalNaive}]                       = init_val;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedCriticalPartialImmunity}]             = init_val;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedCriticalImprovedImmunity}]            = init_val;
        model.populations[{i, mio::osecirvvs::InfectionState::SusceptibleImprovedImmunity}] =
            static_cast<double>(pop_total) / 3;
        model.populations[{i, mio::osecirvvs::InfectionState::SusceptiblePartialImmunity}] =
            static_cast<double>(pop_total) / 3;
        model.populations[{i, mio::osecirvvs::InfectionState::DeadNaive}]            = 0;
        model.populations[{i, mio::osecirvvs::InfectionState::DeadPartialImmunity}]  = 0;
        model.populations[{i, mio::osecirvvs::InfectionState::DeadImprovedImmunity}] = 0;
        model.populations.set_difference_from_group_total<mio::AgeGroup>(
            {i, mio::osecirvvs::InfectionState::SusceptibleNaive}, pop_total);
    }

    const size_t vacc_first                                              = 5;
    const size_t vacc_full                                               = 5;
    model.parameters.get<mio::osecirvvs::ICUCapacity<double>>()          = 100;
    model.parameters.get<mio::osecirvvs::TestAndTraceCapacity<double>>() = 0.0143;
    model.parameters.get<mio::osecirvvs::DailyPartialVaccinations<double>>().resize(mio::SimulationDay(tmax));
    model.parameters.get<mio::osecirvvs::DailyPartialVaccinations<double>>().array().setConstant(vacc_first);
    model.parameters.get<mio::osecirvvs::DailyFullVaccinations<double>>().resize(mio::SimulationDay(tmax));
    model.parameters.get<mio::osecirvvs::DailyFullVaccinations<double>>().array().setConstant(vacc_full);

    auto& contacts       = model.parameters.get<mio::osecirvvs::ContactPatterns<double>>();
    auto& contact_matrix = contacts.get_cont_freq_mat();
    contact_matrix[0].get_baseline().setConstant(0.5);
    contact_matrix[0].get_baseline().diagonal().setConstant(5.0);
    contact_matrix[0].add_damping(0.3, mio::SimulationTime<double>(5.0));

    for (mio::AgeGroup i = 0; i < (mio::AgeGroup)num_agegroups; i++) {
        //times
        model.parameters.get<mio::osecirvvs::TimeExposed<double>>()[i]            = 3.33;
        model.parameters.get<mio::osecirvvs::TimeInfectedNoSymptoms<double>>()[i] = 1.87;
        model.parameters.get<mio::osecirvvs::TimeInfectedSymptoms<double>>()[i]   = 7;
        model.parameters.get<mio::osecirvvs::TimeInfectedSevere<double>>()[i]     = 6;
        model.parameters.get<mio::osecirvvs::TimeInfectedCritical<double>>()[i]   = 7;

        //probabilities
        model.parameters.get<mio::osecirvvs::TransmissionProbabilityOnContact<double>>()[i]  = 0.15;
        model.parameters.get<mio::osecirvvs::RelativeTransmissionNoSymptoms<double>>()[i]    = 0.5;
        model.parameters.get<mio::osecirvvs::RiskOfInfectionFromSymptomatic<double>>()[i]    = 0.0;
        model.parameters.get<mio::osecirvvs::MaxRiskOfInfectionFromSymptomatic<double>>()[i] = 0.4;
        model.parameters.get<mio::osecirvvs::RecoveredPerInfectedNoSymptoms<double>>()[i]    = 0.2;
        model.parameters.get<mio::osecirvvs::SeverePerInfectedSymptoms<double>>()[i]         = 0.1;
        model.parameters.get<mio::osecirvvs::CriticalPerSevere<double>>()[i]                 = 0.1;
        model.parameters.get<mio::osecirvvs::DeathsPerCritical<double>>()[i]                 = 0.1;

        model.parameters.get<mio::osecirvvs::ReducExposedPartialImmunity<double>>()[i]                     = 0.8;
        model.parameters.get<mio::osecirvvs::ReducExposedImprovedImmunity<double>>()[i]                    = 0.331;
        model.parameters.get<mio::osecirvvs::ReducInfectedSymptomsPartialImmunity<double>>()[i]            = 0.65;
        model.parameters.get<mio::osecirvvs::ReducInfectedSymptomsImprovedImmunity<double>>()[i]           = 0.243;
        model.parameters.get<mio::osecirvvs::ReducInfectedSevereCriticalDeadPartialImmunity<double>>()[i]  = 0.1;
        model.parameters.get<mio::osecirvvs::ReducInfectedSevereCriticalDeadImprovedImmunity<double>>()[i] = 0.091;
        model.parameters.get<mio::osecirvvs::ReducTimeInfectedMild<double>>()[i]                           = 0.9;
    }
    model.parameters.get<mio::osecirvvs::Seasonality<double>>() = 0.2;
    return model;
}

template <class Integrator>
auto create_simulation()
{
    auto cfg = mio::benchmark::GraphConfig::initialize(config_path);

    mio::osecirvvs::Model<double> model = create_model(cfg.num_agegroups, cfg.t_max);

    mio::Graph<mio::SimulationNode<double, mio::Simulation<double, mio::osecirvvs::Model<double>>>,
               mio::MobilityEdge<double>>
        g;
    for (size_t county_id = 0; county_id < cfg.num_regions; county_id++) {
        g.add_node(county_id, model, cfg.t0);
        g.nodes()[county_id].property.get_simulation().set_integrator(std::make_shared<Integrator>());
    }

    // Graph is always complete here
    for (size_t county_idx_i = 0; county_idx_i < g.nodes().size(); ++county_idx_i) {
        for (size_t county_idx_j = 0; county_idx_j < g.nodes().size(); ++county_idx_j) {
            if (county_idx_i == county_idx_j)
                continue;
            g.add_edge(
                county_idx_i, county_idx_j,
                Eigen::VectorXd::Constant((size_t)mio::osecirvvs::InfectionState::Count * cfg.num_agegroups, 0.01));
        }
    }

    return mio::make_mobility_sim(cfg.t0, cfg.dt, std::move(g));
}

template <class Integrator>
void init_benchmark(::benchmark::State& state)
{
    for (auto _ : state) {
        // This code gets timed
        create_simulation<Integrator>();
    }
}

template <class Integrator>
void graph_sim_secirvvs(::benchmark::State& state)
{
    mio::set_log_level(mio::LogLevel::critical);
    auto cfg = mio::benchmark::GraphConfig::initialize(config_path);

    for (auto _ : state) {
        // This code gets timed
        auto sim = create_simulation<Integrator>();
        sim.advance(cfg.t_max);
    }
}

// register functions as a benchmarks and set a name
// mitigate influence of cpu scaling
BENCHMARK_TEMPLATE(graph_sim_secirvvs, mio::RKIntegratorCore<double>)->Name("Dummy 1/3");
BENCHMARK_TEMPLATE(graph_sim_secirvvs, mio::RKIntegratorCore<double>)->Name("Dummy 2/3");
BENCHMARK_TEMPLATE(graph_sim_secirvvs, mio::RKIntegratorCore<double>)->Name("Dummy 3/3");
// register functions as a benchmarks and set a name
BENCHMARK_TEMPLATE(init_benchmark, mio::RKIntegratorCore<double>)->Name("Initialize Graph without simulation");
BENCHMARK_TEMPLATE(graph_sim_secirvvs, mio::EulerIntegratorCore<double>)
    ->Name("Graph Simulation - simple explicit euler");
BENCHMARK_TEMPLATE(graph_sim_secirvvs, mio::RKIntegratorCore<double>)->Name("Graph Simulation - adapt_rk");
BENCHMARK_TEMPLATE(graph_sim_secirvvs,
                   mio::ControlledStepperWrapper<double, boost::numeric::odeint::runge_kutta_cash_karp54>)
    ->Name("Graph Simulation - rk_ck54 (boost)");
// BENCHMARK_TEMPLATE(graph_sim_secirvvs, mio::ControlledStepperWrapper<boost::numeric::odeint::runge_kutta_dopri5>)
// ->Name("Graph Simulation - rk_dopri5 (boost)"); // TODO: reenable once boost bug is fixed
BENCHMARK_TEMPLATE(graph_sim_secirvvs,
                   mio::ControlledStepperWrapper<double, boost::numeric::odeint::runge_kutta_fehlberg78>)
    ->Name("Graph Simulation - rkf78 (boost)");
// run all benchmarks
BENCHMARK_MAIN();
