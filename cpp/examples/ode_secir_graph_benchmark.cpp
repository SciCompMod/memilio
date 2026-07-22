/*
* Copyright (C) 2020-2026 MEmilio
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

/**
 * Benchmark comparing the overhead of the model-specific osecirvvs::Simulation::advance()
 * (which uses substeps and calls apply_vaccination / apply_variant / dynamic NPI checks)
 * versus the generic mio::Simulation::advance() (a single integrator call for the entire range).
 */
#include "memilio/compartments/simulation.h"
#include "memilio/math/stepper_wrapper.h"
#include "memilio/mobility/metapopulation_mobility_instant.h"
#include "ode_secir/model.h"

#include "boost/numeric/odeint/stepper/runge_kutta4.hpp"

#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef __linux__
#include <sys/resource.h>
#endif

namespace
{
using FP    = double;
using State = mio::osecir::InfectionState;
using Sim   = mio::osecir::Simulation<FP>;
using Graph = mio::Graph<mio::SimulationNode<FP, Sim>, mio::MobilityEdge<FP>>;
using RK4   = mio::ExplicitStepperWrapper<FP, boost::numeric::odeint::runge_kutta4>;
using Clock = std::chrono::steady_clock;

double seconds(Clock::time_point begin, Clock::time_point end)
{
    return std::chrono::duration<double>(end - begin).count();
}

size_t count_steps(FP tmax, FP dt)
{
    size_t steps = 0;
    for (FP t = 0.0; t < tmax; ++steps) {
        const FP step = t + dt > tmax ? tmax - t : dt;
        t += step;
    }
    return steps;
}

mio::osecir::Model<FP> make_model(size_t num_groups, FP seed_scale)
{
    mio::osecir::Model<FP> model{int(num_groups)};
    auto& p       = model.parameters;
    const FP fact = 1.0 / FP(num_groups);

    p.set<mio::osecir::StartDay<FP>>(0.0);
    p.set<mio::osecir::Seasonality<FP>>(0.0);
    p.get<mio::osecir::TestAndTraceCapacity<FP>>() = 1e12;
    auto& contacts                                 = p.get<mio::osecir::ContactPatterns<FP>>().get_cont_freq_mat();
    contacts[0]                                    = mio::ContactMatrix<FP>(
        Eigen::MatrixX<FP>::Constant(Eigen::Index(num_groups), Eigen::Index(num_groups), 10.0 * fact));

    for (auto age = mio::AgeGroup(0); age < mio::AgeGroup(num_groups); ++age) {
        p.get<mio::osecir::TimeExposed<FP>>()[age]                       = 3.2;
        p.get<mio::osecir::TimeInfectedNoSymptoms<FP>>()[age]            = 2.0;
        p.get<mio::osecir::TimeInfectedSymptoms<FP>>()[age]              = 5.8;
        p.get<mio::osecir::TimeInfectedSevere<FP>>()[age]                = 9.5;
        p.get<mio::osecir::TimeInfectedCritical<FP>>()[age]              = 7.1;
        p.get<mio::osecir::TransmissionProbabilityOnContact<FP>>()[age]  = 0.1;
        p.get<mio::osecir::RelativeTransmissionNoSymptoms<FP>>()[age]    = 0.7;
        p.get<mio::osecir::RecoveredPerInfectedNoSymptoms<FP>>()[age]    = 0.09;
        p.get<mio::osecir::RiskOfInfectionFromSymptomatic<FP>>()[age]    = 0.25;
        p.get<mio::osecir::MaxRiskOfInfectionFromSymptomatic<FP>>()[age] = 0.45;
        p.get<mio::osecir::SeverePerInfectedSymptoms<FP>>()[age]         = 0.2;
        p.get<mio::osecir::CriticalPerSevere<FP>>()[age]                 = 0.25;
        p.get<mio::osecir::DeathsPerSevere<FP>>()[age]                   = 0.0;
        p.get<mio::osecir::DeathsPerCritical<FP>>()[age]                 = 0.3;
        model.populations[{age, State::Exposed}]                         = seed_scale * 100.0 * fact;
        model.populations[{age, State::InfectedNoSymptoms}]              = seed_scale * 50.0 * fact;
        model.populations[{age, State::InfectedSymptoms}]                = seed_scale * 50.0 * fact;
        model.populations[{age, State::InfectedSevere}]                  = seed_scale * 20.0 * fact;
        model.populations[{age, State::InfectedCritical}]                = seed_scale * 10.0 * fact;
        model.populations[{age, State::Recovered}]                       = seed_scale * 10.0 * fact;
        model.populations.set_difference_from_group_total<mio::AgeGroup>({age, State::Susceptible}, 10000.0 * fact);
    }
    model.apply_constraints();
    return model;
}

mio::MobilityParameters<FP> make_mobility(size_t num_groups, size_t degree)
{
    Eigen::VectorX<FP> coefficients = Eigen::VectorX<FP>::Constant(Eigen::Index(num_groups * size_t(State::Count)),
                                                                   degree > 0 ? 0.01 / FP(degree) : 0.0);
    for (size_t age = 0; age < num_groups; ++age) {
        coefficients[Eigen::Index(age * size_t(State::Count) + size_t(State::Dead))] = 0.0;
    }
    return mio::MobilityParameters<FP>(coefficients);
}

void add_edges(Graph& graph, size_t num_nodes, size_t base_degree, size_t extra,
               const mio::MobilityParameters<FP>& base_mobility, const mio::MobilityParameters<FP>& extra_mobility)
{
    // Cyclic successors give a deterministic directed density; wrap first to keep edges sorted.
    for (size_t source = 0; source < num_nodes; ++source) {
        const size_t degree  = base_degree + size_t(source < extra);
        const auto& mobility = source < extra ? extra_mobility : base_mobility;
        const size_t wrapped = source + degree >= num_nodes ? source + degree - num_nodes + 1 : 0;
        for (size_t target = 0; target < wrapped; ++target) {
            graph.add_edge(source, target, mobility);
        }
        for (size_t target = source + 1; target < num_nodes && target <= source + degree; ++target) {
            graph.add_edge(source, target, mobility);
        }
    }
}
} // namespace

int main(int argc, char** argv)
try {
    if (argc < 5 || argc > 6) {
        throw std::runtime_error("usage: ode_secir_graph_benchmark N AGE_GROUPS TMAX DT [STATE.bin]");
    }
    mio::set_log_level(mio::LogLevel::critical);
    const size_t num_nodes  = std::stoull(argv[1]);
    const size_t num_groups = std::stoull(argv[2]);
    const FP tmax           = std::stod(argv[3]);
    const FP dt             = std::stod(argv[4]);
    if (num_nodes < 5 || num_groups == 0 || tmax <= 0.0 || dt <= 0.0) {
        throw std::runtime_error("require N >= 5, AGE_GROUPS > 0, TMAX > 0 and DT > 0");
    }

    const size_t num_edges   = (2 * num_nodes * (num_nodes - 1) + 2) / 5;
    const size_t base_degree = num_edges / num_nodes;
    const size_t extra       = num_edges % num_nodes;
    const auto all_begin     = Clock::now();
    auto model_a             = make_model(num_groups, 1.0);
    auto model_b             = make_model(num_groups, 0.5);
    const auto model_end     = Clock::now();

    Graph graph;
    for (size_t node = 0; node < num_nodes; ++node) {
        graph.add_node(int(node), node % 2 == 0 ? model_a : model_b, 0.0, dt);
    }
    const auto nodes_end = Clock::now();
    for (auto& node : graph.nodes()) {
        node.property.get_simulation().set_integrator_core(std::make_unique<RK4>());
    }
    const auto integrator_end = Clock::now();
    const auto base_mobility  = make_mobility(num_groups, base_degree);
    const auto extra_mobility = make_mobility(num_groups, base_degree + 1);
    add_edges(graph, num_nodes, base_degree, extra, base_mobility, extra_mobility);
    const auto graph_end = Clock::now();

    auto simulation    = mio::make_mobility_sim<FP>(0.0, dt, std::move(graph));
    const auto sim_end = Clock::now();
    simulation.advance(tmax);
    const auto advance_end = Clock::now();
    std::vector<FP> state;
    state.reserve(num_nodes * num_groups * size_t(State::Count));
    for (const auto& node : simulation.get_graph().nodes()) {
        const auto& value = node.property.get_result().get_last_value();
        state.insert(state.end(), value.data(), value.data() + value.size());
    }
    const auto extract_end = Clock::now();
    const FP state_sum     = std::accumulate(state.begin(), state.end(), FP(0));

    if (argc == 6) {
        std::ofstream output(argv[5], std::ios::binary);
        output.write(reinterpret_cast<const char*>(state.data()), std::streamsize(state.size() * sizeof(FP)));
        if (!output) {
            throw std::runtime_error("failed to write state file");
        }
    }
    long max_rss_kb = 0;
#ifdef __linux__
    rusage usage{};
    getrusage(RUSAGE_SELF, &usage);
    max_rss_kb = usage.ru_maxrss;
#endif
    const FP connectivity = FP(num_edges) / FP(num_nodes * (num_nodes - 1));
    std::cout << std::setprecision(12) << "cpp," << num_nodes << ',' << num_groups << ',' << connectivity << ','
              << num_edges << ',' << tmax << ',' << dt << ',' << count_steps(tmax, dt) << ",rk4,"
              << seconds(all_begin, model_end) << ',' << seconds(model_end, nodes_end) << ','
              << seconds(nodes_end, integrator_end) << ',' << seconds(integrator_end, graph_end) << ','
              << seconds(graph_end, sim_end) << ',' << seconds(all_begin, sim_end) << ','
              << seconds(sim_end, advance_end) << ',' << seconds(advance_end, extract_end) << ','
              << seconds(all_begin, extract_end) << ',' << max_rss_kb << ',' << state_sum << '\n';
    return 0;
}
catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 2;
}
